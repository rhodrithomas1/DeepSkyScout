import os
import sqlite3
from math import degrees, atan
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.coordinates import Angle
from astropy.time import Time
import astropy.units as u
from datetime import datetime, timedelta
import numpy as np
from scipy.interpolate import interp1d
import pytz
from functools import lru_cache
import re

from horizon import get_horizon_interpolator

# --- Helper Functions ---
def calculate_fov(sensor_width: float, sensor_height: float, focal_length: float):
    fov_x = degrees(2 * atan(sensor_width / (2 * focal_length)))
    fov_y = degrees(2 * atan(sensor_height / (2 * focal_length)))
    return fov_x, fov_y



def ensure_db_exists(db_path: str):
    if not os.path.exists(db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS objects (
                id INTEGER PRIMARY KEY,
                name TEXT,
                ra REAL,
                dec REAL,
                magnitude REAL,
                size REAL,
                type TEXT,
                narrowband BOOLEAN
            )
            """
        )
        conn.commit()
        conn.close()




def load_objects_from_db(db_path: str, location_name: str = None):
    ensure_db_exists(db_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    if location_name:
        cursor.execute(
            """
            SELECT o.id, o.name, o.ra, o.dec, o.magnitude, o.size, a.max_altitude, o.type, o.narrowband
              FROM objects o
              LEFT JOIN object_altitudes a
                ON o.id = a.object_id AND a.location_name = ?
            """, (location_name,)
        )
    else:
        cursor.execute(
            "SELECT id, name, ra, dec, magnitude, size, NULL, type, narrowband FROM objects"
        )
    rows = cursor.fetchall()
    conn.close()
    return rows




def object_visibility(ra: float, dec: float, lat: float, lon: float, horizon_data):
    local_tz = pytz.timezone("Europe/London")
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    obj = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    # 1) Build times once
    now_local = datetime.now(local_tz)
    dt_list_local = [now_local + timedelta(minutes=i) for i in range(0, 1440, 5)]
    dt_list_utc   = [dt.astimezone(pytz.utc)       for dt in dt_list_local]
    times         = Time(dt_list_utc)

    # 2) Compute alt/az
    altaz = obj.transform_to(AltAz(obstime=times, location=location))
    az_arr = altaz.az.degree
    alt_arr = altaz.alt.degree

    # 3) Horizon interpolation (cached)
    if horizon_data:
        # sort + wrap
        arr = sorted(horizon_data + [(360.0, horizon_data[0][1])])
        az_hz, alt_hz = zip(*arr)
        interp_func = get_horizon_interpolator(tuple(az_hz), tuple(alt_hz))
        horizon_vals = interp_func(az_arr)
    else:
        horizon_vals = np.zeros_like(alt_arr)

    # 4) Visibility mask
    visible = alt_arr > horizon_vals

    # 5) First rise/set detection
    rise = set_ = None
    for i in range(1, len(visible)):
        if visible[i] and not visible[i-1]:
            rise = times[i].datetime.astimezone(local_tz)
        if not visible[i] and visible[i-1]:
            set_ = times[i].datetime.astimezone(local_tz)
        if rise and set_:
            break

    return rise, set_



# --- Camera Database Helper ---
def load_camera_database() -> str:
    from pathlib import Path
    path = Path(__file__).parent / "camera_specs.db"
    if not path.exists():
        conn = sqlite3.connect(path)
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS cameras (
                id INTEGER PRIMARY KEY,
                make TEXT,
                model TEXT,
                pixel_size REAL,
                sensor_width REAL,
                sensor_height REAL,
                resolution_x INTEGER,
                resolution_y INTEGER
            )
            """
        )
        cursor.execute(
            """
            INSERT INTO cameras (make, model, pixel_size, sensor_width, sensor_height, resolution_x, resolution_y)
            VALUES ('Generic','CMOS',4.8,22.3,14.9,6000,4000)
            """
        )
        conn.commit()
        conn.close()
    return str(path)




def calculate_max_altitude(ra, dec, lat, lon):
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    skycoord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    
    # Use timezone-aware UTC time
    now_utc = datetime.now(pytz.utc)
    now = Time(now_utc)

    times = now + np.linspace(0, 1, 100) * u.day
    altitudes = skycoord.transform_to(AltAz(obstime=times, location=location)).alt
    return round(np.max(altitudes).value, 2)





def parse_ra_dec(ra_input, dec_input):
    def _parse(val, is_ra=False):
        if not isinstance(val, str):
            return float(val)
        s = val.strip()
        if re.search(r"[hHdDmMsS]", s):
            return Angle(s).degree
        return float(s)
    return _parse(ra_input, True), _parse(dec_input, False)


def format_ra_hms(ra_deg: float) -> str:
    return Angle(ra_deg * u.deg).to_string(unit=u.hourangle, sep='hms', precision=2, pad=True)


def format_dec_dms(dec_deg: float) -> str:
    return Angle(dec_deg * u.deg).to_string(unit=u.deg, sep='dms', precision=1, alwayssign=True, pad=True)




import os
from astroquery.skyview import SkyView
from astropy import units as u
from astropy.visualization import PercentileInterval, AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from PIL import Image
import numpy as np

def fetch_survey_image(ra_deg, dec_deg,
                       survey="DSS2 Red",
                       width_deg=0.5,
                       height_deg=0.5,
                       object_name: str = None):
    """
    Download (and cache) a small cutout around (ra_deg,dec_deg) from SkyView.
    If object_name is provided, it’s prepended to the cache‐filename.
    """
    # 1) Ensure cache dir exists
    images_dir = os.path.join(os.path.dirname(__file__), "images")
    os.makedirs(images_dir, exist_ok=True)

    # 2) Build a safe filename
    survey_tag = survey.replace(" ", "_")
    ra_tag     = f"{ra_deg:.6f}"
    dec_tag    = f"{dec_deg:.6f}"
    size_tag   = f"{width_deg:.3f}x{height_deg:.3f}"
    if object_name:
        name_tag = object_name.replace(" ", "_")
        fname = f"{name_tag}_{survey_tag}_{ra_tag}_{dec_tag}_{size_tag}.png"
    else:
        fname = f"{survey_tag}_{ra_tag}_{dec_tag}_{size_tag}.png"
    cache_path = os.path.join(images_dir, fname)

    # 3) Return cached if present
    if os.path.exists(cache_path):
        return Image.open(cache_path)

    # 4) Otherwise fetch from SkyView
    pos = f"{ra_deg} {dec_deg}"
    imgs = SkyView.get_images(
        position=pos,
        survey=[survey],
        coordinates="J2000",
        width=width_deg * u.deg,
        height=height_deg * u.deg
    )
    if not imgs:
        raise RuntimeError(f"No image returned for {pos} / {survey}")

    hdu  = imgs[0][0]
    data = hdu.data.astype(float)

    norm   = ImageNormalize(data,
                            interval=PercentileInterval(99.5),
                            stretch=AsinhStretch())
    scaled = (norm(data) * 255).clip(0,255).astype(np.uint8)

    if scaled.ndim == 2:
        rgb = np.stack([scaled]*3, axis=-1)
    else:
        rgb = scaled

    img = Image.fromarray(rgb)

    # 5) Save to cache
    try:
        img.save(cache_path)
    except Exception as e:
        print(f"[WARNING] Could not cache image to {cache_path}: {e}")

    return img
