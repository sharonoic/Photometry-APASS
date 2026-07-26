# Photometry-APASS
A source detection and photometry package built on the APASS catalogue. 

## Installation
Instructions for setting up `apass_phot3.py` and `query_apass2.py`. 

1. Clone the repository <br>
   `git clone https://github.com/sharonoic/Photometry-APASS`
3. Python version
   Requires Python 3.9+. Check with: <br>
   `python3 --version`
4. Install dependencies
   From the directory containing `requirements.txt`: <br>
   `pip install -r requirements.txt` <br>
   This installs: numpy, astropy, photutils, astroquery.
5. Verify the install
   `python3 -c "import numpy, astropy, photutils, astroquery; print('OK')"` <br>
   Should print `OK` with no errors.
6. Network access
   `query_apass2.py` queries Vizier over the network to fetch APASS9. Confirm outbound access to Vizier works before running it:
   `python3 -c "from astroquery.vizier import Vizier; print(Vizier.query_constraints(catalog='APASS9', ra='0'))"`
   If this hangs or times out, check firewall/proxy settings.
7. Directory setup
   Each target needs its own working directory containing:
   - `apass_ref.csv` (produced by `query_apass2.py`
   - the calibrated `.new` FITS frame to be photometered (one per filter)

## How to use
1. Copy `apass_phot3.py` and `query_apass2.py` to the directory containing your image data files
2. Solve images with astrometry.net <br>
   `solve-field -p -D $PWD *.fits` <br>
   (You may have to install astrometry.net if not yet by: `sudo apt install astrometry.net)
3. Query APASS catalogue <br>
   Change the ra and dec fields in `query_apass2.py` to your object's ra and dec <br>
   `python3 query_apass2.py`
4. Run photometry <br>
   Change the ra and dec fields in `apass_phot3.py` to your object's ra and dec <br>
   `python3 apass_phot.py`
5. Check results <br>
   `cat nova_phot_apass.txt`
