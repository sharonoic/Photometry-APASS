import numpy as np
import glob
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry

# ================= Settings =================
TARGET_RA  = 166.730713
TARGET_DEC = -4.117950

N_REF = 7                   # number of ensemble reference stars
MAX_SEP_ARCMIN = 6.0        # max distance of ref stars from target
VMAG_MIN, VMAG_MAX = 12.0, 15.5   # avoid saturation / too faint

AP_ARCSEC = 6.0             # aperture radius (target + refs)
SKY_IN = 10.0               # sky annulus inner (arcsec)
SKY_OUT = 20.0              # sky annulus outer (arcsec)

APPLY_VALIDATION_CORR = False   # with ensemble median, validation is a check;
                                # set True only if you want the single-star corr applied
# ============================================


def measure_star(data_sub, wcs, ra, dec, ap_arcsec, pix_scale):
    """Aperture photometry with local sky annulus. Returns sky-subtracted flux."""
    x, y = wcs.world_to_pixel_values(ra, dec)
    ap = CircularAperture((float(x), float(y)), r=ap_arcsec / pix_scale)
    sky = CircularAnnulus((float(x), float(y)),
                          r_in=SKY_IN / pix_scale, r_out=SKY_OUT / pix_scale)
    flux_raw = aperture_photometry(data_sub, ap)['aperture_sum'][0]
    sky_mean = aperture_photometry(data_sub, sky)['aperture_sum'][0] / sky.area
    return flux_raw - sky_mean * ap.area, ap


# --- Load APASS catalog and select candidate stars ---
refs = Table.read('apass_ref.csv')
print(f"Loaded {len(refs)} APASS stars from apass_ref.csv")

complete = (~refs['Bmag'].mask & ~refs['Vmag'].mask & ~refs["r'mag"].mask
            if hasattr(refs['Bmag'], 'mask') else np.ones(len(refs), bool))
cand = refs[complete]

tgt_coord = SkyCoord(TARGET_RA, TARGET_DEC, unit='deg')
cand_coords = SkyCoord(cand['RAJ2000'], cand['DEJ2000'], unit='deg')
seps = tgt_coord.separation(cand_coords).arcmin

sel = (seps < MAX_SEP_ARCMIN) & (seps > 0.1) \
      & (cand['Vmag'] > VMAG_MIN) & (cand['Vmag'] < VMAG_MAX)
cand = cand[sel]
cand_seps = seps[sel]

# Sort by brightness; take N_REF for ensemble, next one as validation
order = np.argsort(cand['Vmag'])
cand = cand[order]
cand_seps = cand_seps[order]

if len(cand) < 2:
    raise SystemExit(f"Only {len(cand)} usable catalog stars — relax cuts or check query!")

n_ens = min(N_REF, len(cand) - 1)
ens = cand[:n_ens]
val = cand[n_ens]          # first star NOT in the ensemble -> independent validation

print(f"\nEnsemble: {n_ens} stars (V = {ens['Vmag'].min():.2f} to {ens['Vmag'].max():.2f})")
for s in ens:
    print(f"  V={s['Vmag']:.2f}  RA={s['RAJ2000']:.6f}  Dec={s['DEJ2000']:.6f}")
def coord_name(ra, dec):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra, dec, unit='deg')
    return "APASS J" + c.ra.to_string(unit=u.hourangle, sep='', precision=1, pad=True) \
           + c.dec.to_string(sep='', precision=0, alwayssign=True, pad=True)

VAL_NAME = coord_name(float(val['RAJ2000']), float(val['DEJ2000']))
print(f"Validation star: {VAL_NAME}  V={val['Vmag']:.2f}  RA={val['RAJ2000']:.6f}  Dec={val['DEJ2000']:.6f}")
print(f"Target: RA={TARGET_RA:.6f}  Dec={TARGET_DEC:.6f}\n")

results = []

for fname in sorted(glob.glob('*.new')):
    hdul = fits.open(fname)
    data = hdul[0].data.astype(float)
    head = hdul[0].header
    filt = fname.replace('.new', '')[-1]
    wcs = WCS(head)

    if filt not in ('B', 'V', 'R'):
        print(f"{fname}: unrecognized filter '{filt}', skipping")
        continue

    print(f"{fname}: filter={filt}")

    _, median, std = sigma_clipped_stats(data, sigma=3.0)
    data_sub = data - median
    pix_scale = np.mean(np.abs(wcs.pixel_scale_matrix.diagonal())) * 3600
    ny, nx = data.shape

    cat_col = {'B': 'Bmag', 'V': 'Vmag', 'R': "r'mag"}[filt]

    # --- Measure each ensemble star -> per-star zero point ---
    zp_list, used = [], []
    for s in ens:
        x, y = wcs.world_to_pixel_values(s['RAJ2000'], s['DEJ2000'])
        margin = 20
        if not (margin < x < nx - margin and margin < y < ny - margin):
            continue
        flux, _ = measure_star(data_sub, wcs, s['RAJ2000'], s['DEJ2000'],
                               AP_ARCSEC, pix_scale)
        if flux <= 0:
            continue
        zp_i = float(s[cat_col]) + 2.5 * np.log10(flux)
        zp_list.append(zp_i)
        used.append(float(s['Vmag']))

    if len(zp_list) == 0:
        print("  ERROR: no ensemble stars measurable on this image, skipping\n")
        continue

    zp_arr = np.array(zp_list)
    zp_med = np.median(zp_arr)

    # Sigma-clip outliers around the median (needs >=3 stars to be meaningful)
    if len(zp_arr) >= 3:
        mad_std = 1.4826 * np.median(np.abs(zp_arr - zp_med))
        keep = np.abs(zp_arr - zp_med) < max(3 * mad_std, 0.05)
        zp_arr = zp_arr[keep]
        zp_med = np.median(zp_arr)

    zp = zp_med
    zp_mean = np.mean(zp_arr)
    zp_scatter = np.std(zp_arr) if len(zp_arr) > 1 else 0.05
    zp_err = zp_scatter / np.sqrt(len(zp_arr)) if len(zp_arr) > 1 else 0.05

    print(f"  ZP: median={zp:.3f}, mean={zp_mean:.3f}, scatter={zp_scatter:.3f}, "
          f"err={zp_err:.3f} ({len(zp_arr)} stars)")
    resid = np.array(zp_list) - zp
    print("  Per-star ZP residuals: " +
          ", ".join(f"{r:+.3f}" for r in resid))

    # --- Independent validation star ---
    corr = 0.0
    vx, vy = wcs.world_to_pixel_values(val['RAJ2000'], val['DEJ2000'])
    if 20 < vx < nx - 20 and 20 < vy < ny - 20:
        vflux, _ = measure_star(data_sub, wcs, val['RAJ2000'], val['DEJ2000'],
                                AP_ARCSEC, pix_scale)
        if vflux > 0:
            v_meas = -2.5 * np.log10(vflux) + zp
            v_cat = float(val[cat_col])
            corr = v_cat - v_meas
            flag = "  <-- LARGE, check!" if abs(corr) > 0.1 else ""
            print(f"  Validation ({VAL_NAME}): catalog={v_cat:.3f}, "
                  f"meas={v_meas:.3f}, corr={corr:+.3f}{flag}")
        else:
            print("  Validation: negative flux, skipped")
    else:
        print("  Validation: off-chip for this image")

    # --- Target ---
    tflux, tgt_ap = measure_star(data_sub, wcs, TARGET_RA, TARGET_DEC,
                                 AP_ARCSEC, pix_scale)
    if tflux > 0:
        tgt_mag = -2.5 * np.log10(tflux) + zp
        if APPLY_VALIDATION_CORR:
            tgt_mag += corr
        tgt_snr = tflux / np.sqrt(tflux + tgt_ap.area * std ** 2)
        tgt_err = np.sqrt((1.086 / tgt_snr) ** 2 + zp_err ** 2)
        print(f"  Target: mag={tgt_mag:.3f} +/- {tgt_err:.3f}\n")
    else:
        tgt_mag, tgt_err = 99.9, 99.9
        print("  Target: negative flux\n")

    results.append({'file': fname, 'filter': filt,
                    'mag': round(float(tgt_mag), 3),
                    'err': round(float(tgt_err), 3),
                    'zp': round(float(zp), 3),
                    'zp_mean': round(float(zp_mean), 3),
                    'zp_scatter': round(float(zp_scatter), 3),
                    'n_stars': len(zp_arr),
                    'val_corr': round(float(corr), 3),
                    'val_star': VAL_NAME})

# --- Summary ---
print("=== RESULTS ===")
print(f"{'file':<32}{'filt':<6}{'mag':>8}{'err':>7}{'ZPmed':>8}{'ZPmean':>8}{'scat':>7}{'N':>3}{'val':>8}  {'val_star'}")
for r in results:
    print(f"{r['file']:<32}{r['filter']:<6}{r['mag']:8.3f}{r['err']:7.3f}"
          f"{r['zp']:8.3f}{r['zp_mean']:8.3f}{r['zp_scatter']:7.3f}{r['n_stars']:3d}{r['val_corr']:8.3f}  {r['val_star']}")

if results:
    Table(rows=results).write('nova_phot_apass.txt',
                              format='ascii.fixed_width', overwrite=True)
    print("\nSaved to nova_phot_apass.txt")
