"""
Query APASS9 catalog centered on target coordinates and save to CSV.
Usage: python3 query_apass.py
"""

from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

# --- Target field center ---
TARGET_RA = 294.6625057
TARGET_DEC = -39.6630850
SEARCH_RADIUS_ARCMIN = 15
OUTFILE = 'apass_ref.csv'

def main():
    Vizier.ROW_LIMIT = -1
    target = SkyCoord(ra=TARGET_RA, dec=TARGET_DEC, unit='deg')

    print(f"Querying APASS9 centered on RA={TARGET_RA}, Dec={TARGET_DEC} "
          f"within {SEARCH_RADIUS_ARCMIN} arcmin...")

    result = Vizier.query_region(
        target,
        radius=SEARCH_RADIUS_ARCMIN * u.arcmin,
        catalog='APASS9'
    )

    if len(result) == 0:
        print("No results returned — check catalog name or network access.")
        return

    table = result[0]
    print(f"Found {len(table)} stars")

    table.write(OUTFILE, format='csv', overwrite=True)
    print(f"Saved to {OUTFILE}")

    # Quick coverage sanity check
    cat_coords = SkyCoord(table['RAJ2000'], table['DEJ2000'], unit='deg')
    sep = target.separation(cat_coords).arcmin
    print(f"\nStars within 10 arcmin of target: {(sep < 10).sum()}")
    print(f"RA range: {table['RAJ2000'].min():.4f} to {table['RAJ2000'].max():.4f}")
    print(f"Dec range: {table['DEJ2000'].min():.4f} to {table['DEJ2000'].max():.4f}")


if __name__ == '__main__':
    main()
