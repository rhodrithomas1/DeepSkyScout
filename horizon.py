import csv
import re
import matplotlib.pyplot as plt
import numpy as np
import json
import os
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from math import radians
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from astropy.time import Time
from scipy.interpolate import interp1d
from zoneinfo import ZoneInfo
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from functools import lru_cache



CONFIG_FILE = os.path.join(os.path.dirname(__file__), 'config.json')

def load_config():
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, 'r') as f:
                data = json.load(f)
            if isinstance(data, dict):
                return data
        except Exception as e:
            print(f"[ERROR] Failed to load config.json: {e}")
    return {"last_horizon_file": None}

def save_config(config):
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f)

def load_horizon_file(filepath=None):
    """
    Loads a horizon file (space-separated text) and returns parsed data as a list of (azimuth, altitude) tuples.
    If filepath is None, tries to load from the last used path stored in config.json.
    If loading fails, opens a file dialog to let the user pick a new file.
    """
    horizon_data = []
    config = load_config()
    #print("[DEBUG] Config loaded:", config)
    #print("[DEBUG] Config type:", type(config))

    # Step 1: Load from config if no filepath given
    if filepath is None:
        filepath = config.get("last_horizon_file")
        if filepath:
            print(f"[INFO] Attempting to load saved horizon file: {filepath}")
        else:
            print("[INFO] No saved horizon file found.")

    # Step 2: If file does not exist or wasn't found, open file dialog
    if not filepath:
        print("[INFO] Opening file dialog for horizon file...")
        filepath = filedialog.askopenfilename(title="Select Horizon File", filetypes=[("Text Files", "*.txt *.hrz *.csv"), ("All Files", "*.*")])
        if not filepath:
            print("[WARN] No file selected.")
            return []

    # Step 3: Read the file (assuming space-separated values)
    try:
        with open(filepath, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        az = float(parts[0])
                        alt = float(parts[1])
                        horizon_data.append((az, alt))
                        #print(f'horizondata:{horizon_data}')
                    except ValueError:
                        continue

        # Step 4: Save path to config if successful
        if horizon_data:
            config["last_horizon_path"] = filepath
            config["last_horizon_file"] = filepath
            save_config(config)

            #print(f"[INFO] Horizon file loaded and path saved to config: {filepath}")
        else:
            print("[ERROR] Horizon file contained no valid data.")

    except Exception as e:
        print(f"[ERROR] Failed to load horizon file: {e}")

    return horizon_data



def generate_hourly_times(start, end):
    current = start
    times = []
    while current <= end:
        times.append(current)
        current += timedelta(hours=1)
    return Time(times)  # return astropy Time object



def plot_horizon(horizon_data,object_altaz=None,embed_ax=None,object_name=None,show_equatorial=False,show_live_clock=False,original_skycoord=None,observer_location=None,location_name: str = None  ):

    print("[DEBUG] plot_horizon() called")

    if not horizon_data:
        print("[WARNING] No horizon data to plot.")
        return

    az, alt = zip(*sorted(horizon_data))
    az_rad = [radians(a) for a in az]
    r_alt = [a for a in alt]

    ax = embed_ax
    if ax is None:
        print("[DEBUG] No embedded axis provided, creating new figure.")
        fig = plt.Figure(figsize=(5, 4), facecolor='black')  # 🔲 Dark background
        ax = fig.add_subplot(111, polar=True)
    else:
        print("[DEBUG] Using embedded axis, clearing previous content.")
        ax.clear()
        fig = ax.get_figure()

        fig.set_facecolor('black')  # 🔲 Dark background

        for sub_ax in fig.axes:
            if sub_ax != ax and sub_ax.get_label() == '<colorbar>':
                fig.delaxes(sub_ax)

    # 🟦 Horizon line in bright color
    ax.plot(az_rad, r_alt, label="Horizon", color='deepskyblue')

    print("[DEBUG] Horizon profile plotted.")

    if object_altaz is None or len(object_altaz) == 0:
        print("[WARNING] No valid altaz data to plot.")
        return

    try:
        az_deg = object_altaz.az.degree
        alt_deg = object_altaz.alt.degree
        london_tz = ZoneInfo("Europe/London")
        times = [t.to_datetime(timezone=london_tz) for t in object_altaz.obstime]
    except Exception as e:
        print(f"[ERROR] Failed to process object_altaz data: {e}")
        return

    if np.any(np.isnan(alt_deg)) or np.any(np.isnan(az_deg)):
        print("[WARNING] NaN detected in AZ/ALT values.")
        return

    az_rad = [radians(a) for a in az_deg]
    r_alt = [a for a in alt_deg]

    if not (len(az_rad) == len(r_alt) == len(times)):
        print(f"[ERROR] Mismatch in data lengths: az={len(az_rad)}, r_alt={len(r_alt)}, times={len(times)}")
        return

    label = object_name if object_name else "Object Track"
    print(f"[DEBUG] Plotting {label} with colormap:")

    time_progress = np.linspace(0, 1, len(az_rad))
    norm = Normalize(vmin=0, vmax=1)
    cmap = plt.cm.plasma  # 🔥 Good visibility in dark mode

    sc = ax.scatter(az_rad, r_alt, c=time_progress, cmap=cmap, norm=norm, s=20)

    for i, (r, theta, t) in enumerate(zip(r_alt, az_rad, times)):
        if i % 4 == 0:
            ax.text(theta, r, t.strftime("%H:%M"), fontsize=8, ha='center', va='bottom', rotation=0, color='white')  # ⬅️ white text

    # ⏰ Live Sky Clock marker
    if show_live_clock:
        print("[DEBUG] Drawing Live Sky Clock marker...")
        try:
            location = object_altaz.location
            now_utc = Time(datetime.utcnow())
            now_local = now_utc.to_datetime(timezone=ZoneInfo("Europe/London"))
            current_altaz = original_skycoord.transform_to(AltAz(obstime=now_utc, location=location))

            if not np.isnan(current_altaz.alt.deg) and not np.isnan(current_altaz.az.deg):
                az_rad_now = radians(current_altaz.az.deg)
                alt_now = current_altaz.alt.deg
                ax.plot(az_rad_now, alt_now, marker='o', color='lime', markersize=10, label='Now')
                ax.text(az_rad_now, alt_now + 2, now_local.strftime("%H:%M"), color='lime',
                        fontsize=9, ha='center', va='bottom', rotation=0)
                print(f"[DEBUG] Live marker: az={current_altaz.az.deg:.2f}, alt={alt_now:.2f}, time={now_local.strftime('%H:%M')}")
            else:
                print("[WARNING] Live AltAz contains NaN, skipping marker.")
        except Exception as e:
            print(f"[ERROR] Failed to draw Live Sky Clock: {e}")

    # 🎨 Colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    # 🌐 Equatorial Grid (if enabled)
    if show_equatorial:
        ax.grid(False)
        print("[DEBUG] Drawing equatorial grid...")

        # explicit None check instead of `or`
        if observer_location is not None:
            obs = observer_location
        else:
            obs = EarthLocation(lat=51.5, lon=-0.1)

        now = Time(datetime.now())

        # Meridian line
        theta_meridian = radians(180)
        ax.plot(
            [theta_meridian, theta_meridian],
            [0, 90],
            color='darkgreen',
            linewidth=2,
            linestyle='--',
            label='Meridian'
        )

        # Hide default radial labels
        ax.set_yticklabels([])

        # Declination rings
        custom_dec_lines = [-30, 0, 30, 60]
        for dec in custom_dec_lines:
            ra_vals = np.linspace(0, 360, 180)
            coords = SkyCoord(ra=ra_vals * u.deg, dec=dec * u.deg)
            altaz = coords.transform_to(AltAz(obstime=now, location=obs))
            ax.plot(
                np.radians(altaz.az.deg),
                altaz.alt.deg,
                '-',
                color='white',
                linewidth=0.5
            )

        # Label rings at meridian
        for dec in [d for d in custom_dec_lines if d > 0]:
            ra_vals = np.linspace(0, 360, 360)
            coords = SkyCoord(ra=ra_vals * u.deg, dec=dec * u.deg)
            altaz = coords.transform_to(AltAz(obstime=now, location=obs))
            az_deg = altaz.az.deg
            alt_deg = altaz.alt.deg
            idx = np.argmin(np.abs(az_deg - 180))
            ax.text(
                radians(az_deg[idx]),
                alt_deg[idx],
                f"{dec}°",
                fontsize=10,
                color='white',
                ha='center',
                va='center'
            )

        # RA spokes
        for ra in range(0, 360, 30):
            dec_vals = np.linspace(-30, 85, 60)
            coords = SkyCoord(ra=ra * u.deg, dec=dec_vals * u.deg)
            altaz = coords.transform_to(AltAz(obstime=now, location=obs))
            ax.plot(
                np.radians(altaz.az.deg),
                altaz.alt.deg,
                '-',
                color='white',
                linewidth=0.5
            )

    # …rest of your styling & return…
    ax.set_facecolor('black')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlim(90, 0)
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_title(
    f"Horizon Profile and {object_name} Path" +
    (f" in {location_name}" if location_name else "") ,
    color='white'
)

    ax.tick_params(colors='white')
    for lbl in ax.get_yticklabels(): lbl.set_color('white')
    for lbl in ax.get_xticklabels(): lbl.set_color('white')
    return fig






@lru_cache(maxsize=None)
def get_horizon_interpolator(az_horizon: tuple, alt_horizon: tuple):
    return interp1d(az_horizon, alt_horizon, kind='linear', fill_value='extrapolate')

def is_object_visible_through_horizon(
    altaz_path, horizon_data,
    sample_interval_minutes=60,
    min_visible_minutes=30,
    return_mask=False
):
    if not horizon_data:
        return [] if return_mask else False

    # 1) Prepare cached interpolator
    arr = sorted(horizon_data + [(360.0, horizon_data[0][1])])
    az_hz, alt_hz = zip(*arr)
    interp_func  = get_horizon_interpolator(tuple(az_hz), tuple(alt_hz))

    # 2) Vectorized check
    az_arr  = np.array([p.az.deg % 360 for p in altaz_path])
    alt_arr = np.array([p.alt.deg        for p in altaz_path])
    horizon_vals = interp_func(az_arr)
    mask    = alt_arr > horizon_vals

    if return_mask:
        return mask.tolist()

    visible_minutes = int(mask.sum()) * sample_interval_minutes
    return visible_minutes >= min_visible_minutes
