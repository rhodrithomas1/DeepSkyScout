import tkinter as tk
from tkinter import filedialog
import os
import sqlite3
from datetime import datetime, timedelta, date
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body
from astropy.time import Time
import astropy.units as u
from sunrise import get_sun_times
from horizon import load_horizon_file, plot_horizon, is_object_visible_through_horizon, load_config, save_config
from horizon import get_horizon_interpolator
from helpers import calculate_fov, load_objects_from_db, object_visibility,load_camera_database,ensure_db_exists,calculate_max_altitude,parse_ra_dec,format_ra_hms,format_dec_dms,fetch_survey_image
from zoneinfo import ZoneInfo  

from tkinter import simpledialog

from astral import LocationInfo
from astral.sun import sun
from ttkbootstrap.scrolled import ScrolledText
from scipy.interpolate import interp1d

from ttkbootstrap import ttk
import ttkbootstrap as tb
from ttkbootstrap.constants import *
from concurrent.futures import ThreadPoolExecutor, as_completed


from ttkbootstrap.style import Style
import json

from PIL import Image, ImageDraw, ImageTk
from astral import moon
import tkinter.font as tkfont

class DeepSkyApp:
    def load_camera_options(self):
        conn = sqlite3.connect(self.camera_db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT DISTINCT make FROM cameras ORDER BY make")
        makes = [row[0] for row in cursor.fetchall()]
        self.make_menu['values'] = makes
        conn.close()

    def update_models(self, event=None):
        make = self.make_var.get()
        conn = sqlite3.connect(self.camera_db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT model FROM cameras WHERE make = ? ORDER BY model", (make,))
        models = [row[0] for row in cursor.fetchall()]
        self.model_menu['values'] = models
        conn.close()
        if models:
            self.model_var.set(models[0])
            self.load_camera_specs()

    def load_camera_specs(self, event=None):
        make = self.make_var.get()
        model = self.model_var.get()
        conn = sqlite3.connect(self.camera_db_path)
        cursor = conn.cursor()
        cursor.execute("""
            SELECT sensor_width, sensor_height FROM cameras
            WHERE make = ? AND model = ?
        """, (make, model))
        result = cursor.fetchone()
        conn.close()
        if result:
            self.sensor_w_entry.delete(0, tk.END)
            self.sensor_h_entry.delete(0, tk.END)
            self.sensor_w_entry.insert(0, str(result[0]))
            self.sensor_h_entry.insert(0, str(result[1]))


    def load_saved_locations(self):
        self.locations_file = os.path.join(os.path.dirname(__file__), "saved_locations.json")
        if not os.path.exists(self.locations_file):
            with open(self.locations_file, 'w') as f:
                json.dump({}, f)
        with open(self.locations_file, 'r') as f:
            self.saved_locations = json.load(f)
        self.location_menu['values'] = list(self.saved_locations.keys())

    def save_current_location(self):
        name = tk.simpledialog.askstring("Save Location", "Enter a name for this location:")
        if name:
            try:
                lat = float(self.lat_entry.get())
                lon = float(self.lon_entry.get())
                self.saved_locations[name] = {"lat": lat, "lon": lon}
                with open(self.locations_file, 'w') as f:
                    json.dump(self.saved_locations, f, indent=4)
                self.load_saved_locations()
                self.location_var.set(name)
            except ValueError:
                tk.messagebox.showerror("Invalid Input", "Latitude and Longitude must be valid numbers.")


    def use_saved_location(self, event=None):
        name = self.location_var.get()
        if name not in self.saved_locations:
            tk.messagebox.showerror(
                "Location Not Found",
                f"No saved location named '{name}'."
            )
            return

        loc = self.saved_locations[name]
        if 'lat' not in loc or 'lon' not in loc:
            tk.messagebox.showerror(
                "Missing Data",
                f"Saved location '{name}' missing latitude or longitude."
            )
            return

        # 1) Update entry fields
        lat = float(loc['lat'])
        lon = float(loc['lon'])
        self.lat_entry.delete(0, tk.END)
        self.lat_entry.insert(0, str(lat))
        self.lon_entry.delete(0, tk.END)
        self.lon_entry.insert(0, str(lon))

        # 2) Persist for next launch
        cfg = load_config()
        cfg["last_location"] = {"latitude": lat, "longitude": lon}
        save_config(cfg)

        # 3) Update the EarthLocation used by the plot routines
        self.location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

        # 4) Immediately refresh sunrise/sunset
        self.update_sun_times()

        # 5) If you’ve already loaded a horizon, re-filter & redraw
        if self.horizon_data:
            # Kick off the visibility filter (on the UI thread or via your executor)
            self._filter_visible_objects()



    def update_plot_with_equatorial_grid(self):
        if hasattr(self, 'selected_altaz_path') and self.selected_altaz_path is not None:
            plot_horizon(
            self.horizon_data,
            object_altaz=self.selected_altaz_path,
            embed_ax=self.ax,
            object_name=self.selected_object_name,
            show_equatorial=self.show_equatorial_grid.get(),
            observer_location=self.location,
            location_name=self.location_var.get()

        )
            self.canvas.draw()
        else:
            print("[INFO] No object selected yet to update with equatorial grid toggle.")


    def update_clock(self):
        now = datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        self.clock_label.config(text=f"Current UTC Time: {now}")
        self.clock_job = self.root.after(1000, self.update_clock)  # Store job ID


    def __init__(self, root):
        

        style = tb.Style()  
        
        # 2) Force Helvetica 12 on every class you’ll use
       
        style.configure('.', font=('Helvetica', 12))
        root.option_add('*Menu.font', 'Helvetica 12')


        # 3) Now store style & root
        self.style = style
        self.root  = root



        self.root.title("DeepSkyScout")
        self._executor = ThreadPoolExecutor(max_workers=1)
        
        # ───── Database Setup ─────
        self.db_path        = os.path.join(os.path.dirname(__file__), "deepsky_objects.db")

        

        # keep one connection & cursor for the whole app
        self._db_conn   = sqlite3.connect(self.db_path, check_same_thread=False)
        self._db_cursor = self._db_conn.cursor()



        print(f"[DEBUG] Using database at: {self.db_path}")
        self.camera_db_path = load_camera_database()
        ensure_db_exists(self.db_path)
        self.objects        = load_objects_from_db(self.db_path, None)
        self.all_objects    = list(self.objects)

        # ───── Cache RA/Dec & SkyCoord ─────


        self._radec_map     = {}   # oid -> (ra_deg, dec_deg)
        self._skycoord_map  = {}   # oid -> SkyCoord

        for rec in self.all_objects:
            oid, name, raw_ra, raw_dec, *rest = rec
            try:
                ra, dec = parse_ra_dec(raw_ra, raw_dec)
                self._radec_map[oid]    = (ra, dec)
                self._skycoord_map[oid] = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            except Exception:
                # skip any bad rows
                continue


        self.horizon_data   = []

        # ───── Parameter Panels ─────
        params_container = ttk.Frame(root)
        params_container.pack(fill="x", padx=10, pady=5)

        # ─── Telescope parameters (left) ───
        input_frame = tb.LabelFrame(params_container, text="Telescope Parameters")
        input_frame.pack(side="left", fill="both", expand=True, padx=(0,5))

        # Camera Make/Model
        tb.Label(input_frame, text="Camera Make:").grid(row=0, column=2)
        self.make_var  = tk.StringVar()
        self.make_menu = tb.Combobox(input_frame, textvariable=self.make_var, state="readonly")
        self.make_menu.grid(row=0, column=3)
        self.make_menu.bind("<<ComboboxSelected>>", self.update_models)

        tb.Label(input_frame, text="Camera Model:").grid(row=1, column=2)
        self.model_var  = tk.StringVar()
        self.model_menu = tb.Combobox(input_frame, textvariable=self.model_var, state="readonly")
        self.model_menu.grid(row=1, column=3)
        self.model_menu.bind("<<ComboboxSelected>>", self.load_camera_specs)

        self.load_camera_options()

        # FOV entries
        tb.Label(input_frame, text="Focal Length (mm):").grid(row=0, column=0)
        self.focal_entry = tb.Entry(input_frame, bootstyle="secondary")
        self.focal_entry.grid(row=0, column=1)

        tb.Label(input_frame, text="Sensor Width (mm):").grid(row=1, column=0)
        self.sensor_w_entry = tb.Entry(input_frame, bootstyle="secondary")
        self.sensor_w_entry.grid(row=1, column=1)

        tb.Label(input_frame, text="Sensor Height (mm):").grid(row=2, column=0)
        self.sensor_h_entry = tb.Entry(input_frame, bootstyle="secondary")
        self.sensor_h_entry.grid(row=2, column=1)

        # Equipment profiles
        tb.Label(input_frame, text="Equipment Profile:").grid(row=5, column=0)
        self.equipment_var  = tk.StringVar()
        self.equipment_menu = tb.Combobox(input_frame, textvariable=self.equipment_var, state="readonly")
        self.equipment_menu.grid(row=5, column=1)
        self.equipment_menu.bind("<<ComboboxSelected>>", self.load_selected_equipment)

        tb.Button(input_frame, text="Save As...", command=self.save_equipment_dialog).grid(row=6, column=0)
        tb.Button(input_frame, text="Delete",  command=self.delete_equipment).grid(row=6, column=1)

        self.refresh_equipment_menu()

        # Size Filter Controls
        tb.Label(input_frame, text="Min Size (px):").grid(row=7, column=0)
        self.min_pixels_entry = tb.Entry(input_frame, bootstyle="secondary")
        self.min_pixels_entry.grid(row=7, column=1)

        tb.Label(input_frame, text="Max Size (×FOV):").grid(row=8, column=0)
        self.max_ratio_entry = tb.Entry(input_frame, bootstyle="secondary")
        self.max_ratio_entry.grid(row=8, column=1)

        self.min_size_label = tb.Label(input_frame, text="0′ 0″")
        self.min_size_label.grid(row=7, column=2, columnspan=2, sticky="w")
        self.max_size_label = tb.Label(input_frame, text="0′ 0″")
        self.max_size_label.grid(row=8, column=2, columnspan=2, sticky="w")

        self.min_pixels_entry.bind("<KeyRelease>", self.update_size_display)
        self.max_ratio_entry.bind("<KeyRelease>", self.update_size_display)

        tb.Button(input_frame, text="Apply Size Filter", command=self.apply_size_filter)\
        .grid(row=9, column=0, columnspan=2, pady=(5,10))
        tb.Button(input_frame, text="Clear Filter",      command=self.clear_size_filter)\
        .grid(row=9, column=2, columnspan=2, pady=(5,10))

        tb.Button(input_frame, text="Calculate FOV", command=self.update_fov)\
        .grid(row=3, column=0, columnspan=2, pady=5)
        self.fov_label = tb.Label(input_frame, text="FOV: N/A")
        self.fov_label.grid(row=4, column=0, columnspan=2)

        # ─── Observer Location (right) ───
        loc_frame = tb.LabelFrame(params_container, text="Observer Location")
        loc_frame.pack(side="left", fill="both", expand=True, padx=(5,0))

        tb.Label(loc_frame, text="Saved Locations:").grid(row=0, column=2)
        self.location_var  = tk.StringVar()
        self.location_menu = tb.Combobox(loc_frame, textvariable=self.location_var, state="readonly")
        self.location_menu.grid(row=0, column=3)
        # when you pick a location, update lat/lon, sun & moon
        self.location_menu.bind(
            "<<ComboboxSelected>>",
            lambda e: (self.use_saved_location(),
                    self.update_sun_times(),
                    self._update_moon())
        )

        self.load_saved_locations_button = tb.Button(
            loc_frame, text="Save Current Location", command=self.save_current_location
        )
        self.load_saved_locations_button.grid(row=2, column=2, columnspan=2, pady=5)
        self.load_saved_locations()

        tb.Label(loc_frame, text="Latitude (°):").grid(row=0, column=0)
        self.lat_entry = tb.Entry(loc_frame, bootstyle="secondary")
        self.lat_entry.grid(row=0, column=1)

        tb.Label(loc_frame, text="Longitude (°):").grid(row=1, column=0)
        self.lon_entry = tb.Entry(loc_frame, bootstyle="secondary")
        self.lon_entry.grid(row=1, column=1)

        # Auto‐select last saved location
        config = load_config()
        if "last_location" in config:
            last_lat = config["last_location"]["latitude"]
            last_lon = config["last_location"]["longitude"]
            for name, coords in self.saved_locations.items():
                if coords.get("lat")==last_lat and coords.get("lon")==last_lon:
                    self.location_var.set(name)
                    # defer so lat_entry/lon_entry exist
                    self.root.after(0, self.use_saved_location)
                    self.root.after(0, lambda: (self.update_sun_times(), self._update_moon()))
                    break

        # Sun label & button
        self.sun_label = tb.Label(loc_frame, text="Sunrise/Sunset: N/A")
        self.sun_label.grid(row=3, column=0, columnspan=2)
        self.sun_button = tb.Button(
            loc_frame,
            text="Get Sun & Moon",
            command=lambda: (self.update_sun_times(), self._update_moon())
        )
        self.sun_button.grid(row=2, column=0, columnspan=2, pady=5)

        # Horizon status + spinner
        self.horizon_label = tb.Label(loc_frame, text="No horizon file loaded.")
        self.horizon_label.grid(row=4, column=0, columnspan=4, sticky="w", padx=5, pady=(5,0))
        self.progress = ttk.Progressbar(loc_frame, mode="indeterminate")
        self.progress.grid(row=5, column=0, columnspan=2, sticky="ew", pady=(5,10))

        # Moon Phase widget
        moon_frame = tb.LabelFrame(params_container, text="Moon Phase")
        moon_frame.pack(side="left", fill="both", expand=False, padx=(5,0))
        self.moon_canvas = tk.Canvas(
            moon_frame, width=200, height=200,
            bg='#333333', highlightthickness=0
        )
        self.moon_canvas.pack(padx=5, pady=(5,0))
        self.moon_label = tk.Label(
            moon_frame, fg='white', bg='#333333', font=("Segoe UI", 10)
        )
        self.moon_label.pack(pady=(0,5))

        # Kick off initial sun & moon
        self.root.after(0, lambda: (self.update_sun_times(), self._update_moon()))

        # Auto-load last horizon file
        # Auto-load last horizon file
        if "last_horizon_file" in config:
            try:
                path = config["last_horizon_file"]
                self.horizon_data = load_horizon_file(path)
                self.horizon_label.config(
                    text=f"Loaded {len(self.horizon_data)} points from last used horizon file."
                )
                # Kick off the same background filter → UI refresh
                def _start_initial_filter():
                    fut = self._executor.submit(self._filter_visible_objects)
                    fut.add_done_callback(
                        lambda f: self.root.after(0, self._after_filter_done)
                    )
                self.root.after(0, _start_initial_filter)
            except Exception as e:
                print(f"[ERROR] Could not load saved horizon file: {e}")

        # Auto-load last equipment profile
        last_eq = config.get("last_equipment_profile")
        if last_eq and last_eq in config.get("equipment_profiles", {}):
            self.equipment_var.set(last_eq)
            self.load_selected_equipment()




        # ───── Middle Section ─────
        middle_frame = ttk.Frame(root)
        middle_frame.pack(fill="both", expand=False, padx=10, pady=5)

        # Info box (left)
        self.info_text = ScrolledText(middle_frame, height=8, wrap="word", bootstyle="dark")
        self.info_text.pack(side="left", fill="y", padx=(0,10), pady=5)

        # Sky plot (right)
        self.plot_frame = tb.LabelFrame(middle_frame, text="Sky Plot")
        self.plot_frame.pack(side="left", fill="both", expand=True)
        self.fig = plt.Figure(figsize=(5,4))
        self.ax  = self.fig.add_subplot(111, polar=True)
        self.ax.set_facecolor('black')
        self.fig.patch.set_facecolor('black')
        self.ax.tick_params(colors='white')
        for lbl in (*self.ax.get_xticklabels(), *self.ax.get_yticklabels()):
            lbl.set_color('white')
        self.ax.set_theta_zero_location('N')
        self.ax.set_theta_direction(-1)
        self.ax.set_yticks([80,60,40,20,0])
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

     

        
        # ─── After you create self.canvas ───
        # build a horizontal row of controls
        control_frame = ttk.Frame(root)
        control_frame.pack(fill='x', padx=10, pady=(0,10))

        self.show_equatorial_grid = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            control_frame,
            text="Show Equatorial Grid",
            variable=self.show_equatorial_grid,
            command=self.update_plot_with_equatorial_grid
        ).pack(side='left', padx=5)

        self.show_live_sky_clock = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            control_frame,
            text="Show Live Position",
            variable=self.show_live_sky_clock,
            command=self.refresh_plot_with_clock
        ).pack(side='left', padx=5)

        # … after setting up Equatorial & Live checkbuttons …
        self.show_all_objects = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            control_frame,
            text="Show All Objects",
            variable=self.show_all_objects,
            command=self._toggle_show_all_objects
        ).pack(side='left', padx=5)


        # ─── after your existing Show All Objects checkbox ───
        self.show_image_popup = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            control_frame,
            text="Show Object Image",
            variable=self.show_image_popup
        ).pack(side='left', padx=5)



        # ─── Then build your Treeview as before ───
        self.tree = ttk.Treeview(
            root,
            columns=("Name","RA (°)","Dec (°)","Mag","Size (arcmin)",
                     "Max Alt (°)","Type","Narrowband","Visible Time (min)"),
            show="headings"
        )
        # … set up headings, bind, pack, etc. …

        # configure a “greyed” tag
        self.tree.tag_configure('greyed', foreground='red')




        for col in self.tree["columns"]:
            self.tree.heading(col, text=col,
                command=lambda _col=col: self._sort_by_column(_col, False))
        self.tree.pack(fill="both", expand=True, padx=10, pady=5)
        self.tree.bind("<<TreeviewSelect>>", self.on_object_select)
        
        

        # in your __init__, right after you create self.tree:
        self._base_headings = { col: col for col in self.tree["columns"] }


        # ───── Menu Bar + Clock ─────
        menubar_frame = ttk.Frame(root)
        menubar_frame.pack(fill="x", padx=10, pady=(5,0))
        self.clock_label = tb.Label(
            menubar_frame,
            text=f"Current UTC Time: {datetime.utcnow():%Y-%m-%d %H:%M:%S}"
        )
        self.clock_label.pack(side="right")
        self.update_clock()

        self.load_horizon_button = tb.Button(
            menubar_frame,
            text="Load Horizon File",
            command=self.load_horizon
        )
        self.load_horizon_button.pack(side="right", padx=(0,10))

        menubar = tk.Menu(root)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Load Database",     command=self.load_db)
        file_menu.add_command(label="Load Horizon File", command=self.load_horizon)
        # ─── after your existing File → Load Horizon DB entries ───
        file_menu.add_separator()
        file_menu.add_command(
            label="Download Images",
            command=self.download_images
        )

        
        
        menubar.add_cascade(label="File", menu=file_menu)

        theme_menu = tk.Menu(menubar, tearoff=0)
        for theme in Style().theme_names():
            theme_menu.add_command(label=theme, command=lambda t=theme: self.set_theme(t))
        menubar.add_cascade(label="Theme", menu=theme_menu)
        root.config(menu=menubar)

        # Finalize
        self.refresh_table()
        self.clock_job = None
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)



    def load_db(self):
        self.db_path = filedialog.askopenfilename(filetypes=[("SQLite DB", "*.db")])
        if self.db_path:
            ensure_db_exists(self.db_path)
            self.objects = load_objects_from_db(self.db_path)
            self.refresh_table()


    from helpers import parse_ra_dec, format_ra_hms, format_dec_dms

 
    def refresh_table(self):
        self.tree.delete(*self.tree.get_children())
        for obj in self.objects:
            # obj is a 10-tuple: (id, name, ra, dec, mag, size, max_alt, type, narrowband, vis_min)
            try:
                raw_ra, raw_dec = obj[2], obj[3]
                ra_deg, dec_deg = parse_ra_dec(raw_ra, raw_dec)
                ra_str  = format_ra_hms(ra_deg)
                dec_str = format_dec_dms(dec_deg)
            except Exception:
                continue

            name    = obj[1]
            mag     = obj[4]
            size    = obj[5]
            max_alt = obj[6]
            typ     = obj[7]
            nb      = obj[8]

            # Convert visible minutes into "Hh Mm"
            vis_str = ""
            if len(obj) > 9 and obj[9] is not None:
                total_mins = int(obj[9])
                hours = total_mins // 60
                mins  = total_mins % 60
                if hours > 0:
                    vis_str = f"{hours}h {mins}m"
                else:
                    vis_str = f"{mins}m"

            values = (name, ra_str, dec_str, mag, size, max_alt, typ, nb, vis_str)
            self.tree.insert("", "end", values=values)





    def update_fov(self):
        try:
            # 1) Read inputs
            f = float(self.focal_entry.get())
            w = float(self.sensor_w_entry.get())
            h = float(self.sensor_h_entry.get())
            fov_x_deg, fov_y_deg = calculate_fov(w, h, f)

            # 2) Helper to convert a decimal‐degree value to (°,' ,")
            def deg_to_dms(deg: float):
                d = int(deg)
                m_float = (deg - d) * 60
                m = int(m_float)
                s = (m_float - m) * 60
                return d, m, s

            x_d, x_m, x_s = deg_to_dms(fov_x_deg)
            y_d, y_m, y_s = deg_to_dms(fov_y_deg)

            # 3) Format both as D°M′S″
            fov_x_str = f"{x_d}°{x_m}′{x_s:.0f}″"
            fov_y_str = f"{y_d}°{y_m}′{y_s:.0f}″"

            # 4) Update the label
            self.fov_label.config(
                text=f"FOV: {fov_x_str} × {fov_y_str}"
            )
        except ValueError:
            self.fov_label.config(text="Invalid input")




    def update_sun_times(self):
        try:
            lat = float(self.lat_entry.get())
            lon = float(self.lon_entry.get())
            from datetime import datetime, timedelta

            now = datetime.utcnow()
            today = date.today()
            if now.time() > datetime.combine(today, datetime.min.time()).replace(hour=12).time():
                obs_date = today
            else:
                obs_date = today - timedelta(days=1)

            city = LocationInfo(name="Custom", region="Earth", timezone="UTC", latitude=lat, longitude=lon)
            s = sun(city.observer, date=obs_date, tzinfo="UTC")
            sunrise = s['sunrise'].time()
            sunset = s['sunset'].time()
            self.sun_label.config(text=f"Sunrise: {sunrise}, Sunset: {sunset}")
        except Exception as e:
            self.sun_label.config(text=f"Error: {str(e)}")

   
    def _select_horizon_file(self) -> str | None:
        """Prompt the user to pick a horizon file and return its path."""
        path = filedialog.askopenfilename(
            initialdir=os.path.dirname(__file__),
            filetypes=[("Horizon Files", "*.hrz *.csv *.txt")]
        )
        return path or None

    def _load_and_save_horizon(self, file_path: str) -> None:
        """Load horizon_data, update config, and refresh the horizon label."""
        # Load the raw horizon points
        self.horizon_data = load_horizon_file(file_path)

        # Persist to config.json
        cfg = load_config()
        cfg["last_horizon_file"] = file_path
        save_config(cfg)

        # Update the UI
        self.horizon_label.config(
            text=f"Loaded {len(self.horizon_data)} points from horizon file."
        )


    from concurrent.futures import ThreadPoolExecutor, as_completed
    from astropy.coordinates import SkyCoord, AltAz, EarthLocation
    from astropy.time import Time
    import astropy.units as u
    import numpy as np
    from datetime import datetime, timedelta
    from sunrise import get_sun_times
    from helpers import (
        get_horizon_interpolator,
        parse_ra_dec,
        calculate_max_altitude,
        load_objects_from_db
    )




    def _filter_visible_objects(self) -> None:
        import numpy as np
        import astropy.units as u
        from astropy.coordinates import SkyCoord, AltAz, EarthLocation
        from astropy.time import Time
        from datetime import datetime, timedelta
        from sunrise import get_sun_times
        from helpers import (
            get_horizon_interpolator,
            parse_ra_dec,
            calculate_max_altitude,
            load_objects_from_db
        )

        # 1) Observer position & loc_name
        lat = float(self.lat_entry.get().strip())
        lon = float(self.lon_entry.get().strip())
        loc_name = (self.location_var.get().strip() or f"{lat}_{lon}") \
                   .replace(",", "_").replace(" ", "_")

        # 2) Prepare DB table & preload existing max‐alts
        self._db_cursor.execute("""
          CREATE TABLE IF NOT EXISTS object_altitudes (
            object_id     INTEGER,
            location_name TEXT,
            max_altitude  REAL,
            PRIMARY KEY(object_id, location_name)
          )
        """)
        self._db_cursor.execute(
          "SELECT object_id, max_altitude FROM object_altitudes WHERE location_name=?",
          (loc_name,)
        )
        existing_alt = dict(self._db_cursor.fetchall())

        # 3) Build night‐time Time array
        sunrise, sunset = get_sun_times(lat, lon)
        # 1a) Determine next sunrise
        # If sunrise comes *before* sunset today, roll it to tomorrow
        if sunrise <= sunset:
            next_sunrise = sunrise + timedelta(days=1)
        else:
            next_sunrise = sunrise

        night_start = datetime.combine(sunset.date(), sunset.time())
        night_end   = datetime.combine(next_sunrise.date(), next_sunrise.time())

        # 1b) Compute actual night length in minutes
        total_night_minutes = int((night_end - night_start).total_seconds() // 60)

        # You can *coarsen* this resolution to speed up:
        resolution = 30  # minutes per sample; try 10 or 15 if still slow
        total_min = int((night_end - night_start).total_seconds() // 60)
        times = Time([
          night_start + timedelta(minutes=resolution * i)
          for i in range(0, total_min // resolution + 1)
        ])

        # 4) One‐time horizon interpolator
        observer_loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
        az_hz, alt_hz = zip(*sorted(self.horizon_data + [(360.0, self.horizon_data[0][1])]))
        interp_hz = get_horizon_interpolator(tuple(az_hz), tuple(alt_hz))

        # 5) Load all objects
        objs = load_objects_from_db(self.db_path, loc_name)

        # 6) Per‐object work
        def process_one(entry):
            oid, name, raw_ra, raw_dec, mag, size, db_alt, otype, narrow = entry
            # new:
            ra, dec = self._radec_map[oid]
            coord    = self._skycoord_map[oid]

            # vector transform over times only
            altaz = coord.transform_to(AltAz(obstime=times, location=observer_loc))
            az_vals = altaz.az.deg % 360
            alt_vals = altaz.alt.deg
            # never exceed the true night length
            vis_min = int((alt_vals > interp_hz(az_vals)).sum()) * resolution
            vis_min = min(vis_min, total_night_minutes)
  
            if vis_min < 30:
                return None
            # determine max_alt
            if oid in existing_alt:
                max_alt = existing_alt[oid]
                newrow = None
            else:
                max_alt = calculate_max_altitude(ra, dec, lat, lon)
                newrow = (oid, loc_name, max_alt)
            # return 10‐tuple plus new row
            return (
              (oid, name, ra, dec, mag, size, max_alt, otype, narrow, vis_min),
              newrow
            )

        # 7) Run sequentially (or in threads if you prefer)
        filtered = []
        new_rows = []
        # Single‐threaded loop is often fastest for Astropy transforms:
        for entry in objs:
            res = process_one(entry)
            if res:
                tpl, nr = res
                filtered.append(tpl)
                if nr:
                    new_rows.append(nr)

        # 8) Bulk‐insert new altitudes
        if new_rows:
            self._db_cursor.executemany(
              "INSERT OR REPLACE INTO object_altitudes (object_id, location_name, max_altitude) VALUES (?,?,?)",
              new_rows
            )
            self._db_conn.commit()

        # 9) Sort & cache
        self.objects = sorted(filtered, key=lambda x: x[6], reverse=True)
        self._visible_objects = list(self.objects)














    def load_horizon(self):
        file_path = self._select_horizon_file()
        if not file_path:
            return

        # Load & cache the horizon immediately (fast)
        self._load_and_save_horizon(file_path)

        # Disable the button so user can’t click again
        self.progress.start()
        self.load_horizon_button.config(state="disabled")
        # Optionally start a spinner/progressbar here…

        # Submit the heavy work
        future = self._executor.submit(self._filter_visible_objects)

        # When it’s done, re-enable the button and stop spinner on the UI thread
        def _on_done(fut):
            # any exception will be re-raised here if you want to catch it
            try:
                fut.result()
            except Exception as e:
                print("[ERROR] Background filtering failed:", e)
            # Schedule our UI‐updates back on the main Tk thread
            self.root.after(0, self._after_filter_done)

        future.add_done_callback(_on_done)










    def _after_filter_done(self):
        # 1) Re‐enable UI
        self.progress.stop()
        self.load_horizon_button.config(state="normal")

        # 2) Clear any existing arrows in all headers
        for col in self.tree["columns"]:
            self.tree.heading(col,
                text=col,
                command=lambda _col=col: self._sort_by_column(_col, False)
            )

        # 3) Refresh & then default‐sort by Visible Time descending
        self.refresh_table()
        self._sort_by_column("Visible Time (min)", reverse=True)

        # 4) Re‐write that header with a ▼ glyph
        self.tree.heading(
            "Visible Time (min)",
            text="Visible Time (min) ▼",
            command=lambda: self._sort_by_column("Visible Time (min)", False)
        )

        print("[INFO] Visibility filtering complete.")













    def plot_horizon(self):
        if not hasattr(self, 'selected_altaz_path') or self.selected_altaz_path is None:
            print("[ERROR] No selected AltAz path to plot.")
            return

        print("[DEBUG] Plotting selected object horizon path.")
        plot_horizon(self.horizon_data, self.selected_altaz_path, embed_ax=self.ax, object_name=self.selected_object_name)
        self.canvas.draw()


    from imaging import fetch_survey_image  # make sure this is at the top of gui.py

    def on_object_select(self, event):
        print("[DEBUG] on_object_select() called")

        # Ensure horizon data is loaded
        if not self.horizon_data:
            self.info_text.delete("1.0", tk.END)
            self.info_text.insert(tk.END, "Load horizon data first.")
            return

        # Get selected tree item
        selected = self.tree.focus()
        if not selected:
            print("[ERROR] No tree item selected.")
            return

        # ─── Parse the selection ───
        try:
            values = self.tree.item(selected, 'values')
            print(f"[DEBUG] Selected values: {values}")

            if len(values) < 3:
                raise ValueError("Expected at least RA and Dec in selection")

            name    = str(values[0])
            raw_ra  = values[1]
            raw_dec = values[2]
            ra, dec = parse_ra_dec(raw_ra, raw_dec)
            self.selected_object_name = name
            print(f"[DEBUG] Parsed object: {name}, RA: {ra}, Dec: {dec}")
        except Exception as e:
            self.info_text.delete("1.0", tk.END)
            self.info_text.insert(tk.END, f"Error parsing selection: {e}")
            print(f"[ERROR] Error parsing selection: {e}")
            return

        # ─── Compute visibility & plot horizon ───
        try:
            lat = float(self.lat_entry.get())
            lon = float(self.lon_entry.get())
            local_tz = ZoneInfo("Europe/London")

            # Compute rise/set times above custom horizon
            rise, set_ = object_visibility(ra, dec, lat, lon, self.horizon_data)
            if rise and set_:
                rise = rise.astimezone(local_tz)
                set_ = set_.astimezone(local_tz)

            # Build AltAz track for plotting
            location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
            obj      = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

            sunrise, sunset = get_sun_times(lat, lon)
            sunrise = sunrise.astimezone(local_tz)
            sunset  = sunset.astimezone(local_tz)

            night_start = datetime.combine(sunset.date(), sunset.time(), tzinfo=local_tz)
            night_end   = datetime.combine(sunset.date() + timedelta(days=1),
                                        sunrise.time(), tzinfo=local_tz)

            # Determine sample times
            if rise and set_ and rise < night_end and set_ > night_start:
                dt_list = [
                    t for t in
                    (rise + timedelta(hours=i) for i in range(int((set_ - rise).total_seconds() / 3600) + 1))
                    if night_start <= t <= night_end
                ]
            else:
                print("[WARNING] No valid rise/set, falling back to 30 min sampling")
                total = int((night_end - night_start).total_seconds() // 60)
                dt_list = [
                    night_start + timedelta(minutes=i)
                    for i in range(0, total, 30)
                ]

            if not dt_list:
                self.info_text.delete("1.0", tk.END)
                self.info_text.insert(tk.END, "No valid time points for plotting.")
                return

            times = Time(dt_list)
            altaz = obj.transform_to(AltAz(obstime=times, location=location))
            self.selected_altaz_path = altaz

            # Compute meridian (highest) time
            now      = datetime.now(tz=local_tz)
            full_day = Time([now + timedelta(minutes=15*i) for i in range(96)])
            alts     = obj.transform_to(AltAz(obstime=full_day, location=location)).alt
            mid_idx  = int(np.argmax(alts))
            meridian = full_day[mid_idx].to_datetime(timezone=local_tz)

            # Update info box
            self.info_text.delete("1.0", tk.END)
            self.info_text.insert(tk.END, f"Object: {name}\n")
            self.info_text.insert(tk.END, f"Rises above custom horizon: {rise}\n")
            self.info_text.insert(tk.END, f"Sets below custom horizon:  {set_}\n")
            self.info_text.insert(tk.END, f"Meridian (highest point):  {meridian:%Y-%m-%d %H:%M:%S}\n")

            # Plot on the embedded canvas
            fig = plot_horizon(
                self.horizon_data,
                object_altaz=altaz,
                embed_ax=self.ax,
                object_name=name,
                show_equatorial=self.show_equatorial_grid.get(),
                show_live_clock=self.show_live_sky_clock.get(),
                original_skycoord=obj,
                observer_location=self.location,
                location_name=self.location_var.get()
            )
            if fig:
                self.canvas.draw()
                print("[DEBUG] Canvas redrawn")
            else:
                print("[ERROR] No figure returned from plot_horizon.")

        except Exception as e:
            self.info_text.delete("1.0", tk.END)
            self.info_text.insert(tk.END, f"Error computing visibility: {e}")
            print(f"[ERROR] Exception in on_object_select: {e}")
            return

        # ─── Show SkyView popup if enabled ───
        if self.show_image_popup.get():
            try:
                # Size in arc-minutes is in column index 4 of values
                size_arcmin = float(values[4] or 0)
            except ValueError:
                size_arcmin = 0.0
            # Pad by 20% and ensure a minimum of 0.1°
            size_deg = max(size_arcmin / 60.0 * 1.2, 0.1)

            try:
                img = fetch_survey_image(
                    ra_deg=ra,
                    dec_deg=dec,
                    survey="DSS2 Red",
                    width_deg=size_deg,
                    height_deg=size_deg,
                    object_name=name
                )

                popup = tk.Toplevel(self.root)
                popup.title(f"{name} ({size_arcmin:.1f}′) from SkyView")
                photo = ImageTk.PhotoImage(img)
                lbl   = tk.Label(popup, image=photo)
                lbl.image = photo
                lbl.pack(padx=5, pady=5)

            except Exception as ie:
                tk.messagebox.showerror(
                    "Image Load Error",
                    f"Could not fetch image for {name}:\n{ie}"
                )












    

    def refresh_plot_with_clock(self):
        self.on_object_select(None)


    def save_equipment_dialog(self):
        name = simpledialog.askstring("Save Equipment", "Enter a name for this setup:")
        if not name:
            return
        config = load_config()
        profiles = config.setdefault("equipment_profiles", {})
        profiles[name] = {
            "focal_length": self.focal_entry.get().strip(),
            "sensor_width": self.sensor_w_entry.get().strip(),
            "sensor_height": self.sensor_h_entry.get().strip(),
            "camera_make": self.make_var.get(),
            "camera_model": self.model_var.get()
        }
        save_config(config)
        self.refresh_equipment_menu()
        self.equipment_var.set(name)
        print(f"[INFO] Saved equipment profile: {name}")


    def delete_equipment(self):
        name = self.equipment_var.get()
        if not name:
            return
        config = load_config()
        if "equipment_profiles" in config and name in config["equipment_profiles"]:
            del config["equipment_profiles"][name]
            save_config(config)
            self.refresh_equipment_menu()
            self.equipment_var.set("")
            print(f"[INFO] Deleted equipment profile: {name}")


    def load_selected_equipment(self, event=None):
        name = self.equipment_var.get()
        config = load_config()
        profile = config.get("equipment_profiles", {}).get(name, {})

        # 1) Set the camera dropdowns (this will auto-load DB specs)
        self.make_var.set(profile.get("camera_make", ""))
        self.update_models()              # repopulates model_menu & calls load_camera_specs()
        self.model_var.set(profile.get("camera_model", ""))

        # 2) Now restore exactly what you saved in the profile
        self.focal_entry.delete(0, tk.END)
        self.focal_entry.insert(0, profile.get("focal_length", ""))

        self.sensor_w_entry.delete(0, tk.END)
        self.sensor_w_entry.insert(0, profile.get("sensor_width", ""))

        self.sensor_h_entry.delete(0, tk.END)
        self.sensor_h_entry.insert(0, profile.get("sensor_height", ""))

        # 3) Remember this as the last‐used profile
        config["last_equipment_profile"] = name
        save_config(config)

        print(f"[INFO] Loaded equipment profile: {name}")




    def refresh_equipment_menu(self):
        config = load_config()
        profile_names = list(config.get("equipment_profiles", {}).keys())
        self.equipment_menu["values"] = profile_names


    def on_closing(self):
        if self.clock_job:
            self.root.after_cancel(self.clock_job)
        self.root.destroy()
    

    def set_theme(self, theme_name):
        try:
            style = tb.Style()
            style.theme_use(theme_name)
            # re‐apply our Helvetica‐12 override
            for cls in (
                'TLabel','TButton','TEntry','TCombobox',
                'Treeview','Treeview.Heading','Menu'
            ):
                style.configure(cls, font=('Helvetica', 12))
            save_config({ "theme": theme_name })
        except Exception as e:
            print("Theme switch failed:", e)



    def _sort_by_column(self, col: str, reverse: bool) -> None:
        # 1) Sort your data exactly as before...

        numeric_cols = {
            "RA (°)", "Dec (°)", "Mag",
            "Size (arcmin)", "Max Alt (°)",
            "Visible Time (min)"
        }
        full_index = self.tree["columns"].index(col) + 1
        if col in numeric_cols:
            keyfunc = lambda obj: obj[full_index] or 0.0
        else:
            keyfunc = lambda obj: str(obj[full_index]).lower()
        self.objects.sort(key=keyfunc, reverse=reverse)

        # 2) Refresh the table rows
        self.refresh_table()

        # 3) Now update all the headings to show the little arrow
        for c in self.tree["columns"]:
            base = self._base_headings[c]
            if c == col:
                arrow = "▼" if reverse else "▲"
                self.tree.heading(c, text=f"{base} {arrow}",
                                command=lambda _c=c, _r=not reverse: self._sort_by_column(_c, _r))
            else:
                # strip any arrow on the others
                self.tree.heading(c, text=base,
                                command=lambda _c=c: self._sort_by_column(_c, False))







    def _get_camera_resolution(self) -> tuple[int,int]:
        """Return (res_x, res_y) for the selected camera."""
        conn   = sqlite3.connect(self.camera_db_path)
        cursor = conn.cursor()
        make   = self.make_var.get()
        model  = self.model_var.get()
        cursor.execute("""
            SELECT resolution_x, resolution_y
              FROM cameras
             WHERE make=? AND model=?
        """, (make, model))
        row = cursor.fetchone()
        conn.close()
        if row and row[0] and row[1]:
            return row[0], row[1]
        raise ValueError("Camera resolution not found.")
    


    def apply_size_filter(self) -> None:
        """
        Filter self.objects by angular size in ARC-MINUTES:
        • lower_arcmin = (min_pixels * FOVx / res_x) * 60
        • upper_arcmin = (max_ratio   * FOVx)       * 60
        """
        try:
            # 1) Read FOV_x (in degrees)
            f = float(self.focal_entry.get())
            w, h     = float(self.sensor_w_entry.get()), float(self.sensor_h_entry.get())
            fov_x_deg, _ = calculate_fov(w, h, f)

            # 2) Read camera resolution
            res_x, _ = self._get_camera_resolution()

            # 3) Compute pixel scale (deg per pixel) and convert to arcmin/pixel
            scale_deg_per_px  = fov_x_deg / res_x
            scale_arcmin_per_px = scale_deg_per_px * 60.0

            # 4) Read filter inputs
            min_px     = int(self.min_pixels_entry.get() or 0)
            max_ratio  = float(self.max_ratio_entry.get() or 0)

            # 5) Compute bounds in ARC-MINUTES
            lower_arcmin = min_px * scale_arcmin_per_px
            upper_arcmin = max_ratio * fov_x_deg * 60.0

            # 6) Filter the master list (size is stored in arcminutes)
            base_list = getattr(self, "_visible_objects", self.objects)
            filtered = [
                obj for obj in base_list
                if lower_arcmin <= obj[5] <= upper_arcmin
            ]

            # 7) Update and redraw
            self.objects = filtered
            self.refresh_table()

        except Exception as e:
            tk.messagebox.showerror(
                "Filter Error",
                f"Could not apply size filter:\n{e}"
            )





    def clear_size_filter(self):
        # Restore the post-visibility list, not the original DB load
        self.objects = list(getattr(self, "_visible_objects", self.objects))
        self.refresh_table()



    def update_size_display(self, event=None):
        """
        Compute and display the min/max angular sizes next to the pixel/ratio entries,
        with arc-second resolution to 2 decimal places.
        """
        try:
            # 1) Compute horizontal FOV in degrees
            f = float(self.focal_entry.get())
            w = float(self.sensor_w_entry.get())
            fov_x_deg, _ = calculate_fov(w, float(self.sensor_h_entry.get()), f)

            # 2) Get horizontal resolution for pixel scale
            res_x, _ = self._get_camera_resolution()
            scale_deg_per_px     = fov_x_deg / res_x
            scale_arcmin_per_px  = scale_deg_per_px * 60.0

            # — Lower bound in arc-minutes & arc-seconds —
            min_px = int(self.min_pixels_entry.get() or 0)
            lower_arcmin = min_px * scale_arcmin_per_px               # arc-min total
            lower_arcmin_int = int(lower_arcmin)                      # whole arc-min
            lower_arcsec     = (lower_arcmin - lower_arcmin_int) * 60 # remaining arc-sec
            # format to 2 decimal places for arc-seconds
            self.min_size_label.config(
                text=f"{lower_arcmin_int}′ {lower_arcsec:.2f}″"
            )

            # — Upper bound in arc-minutes & arc-seconds —
            max_ratio = float(self.max_ratio_entry.get() or 0)
            upper_arcmin = max_ratio * fov_x_deg * 60.0
            upper_arcmin_int = int(upper_arcmin)
            upper_arcsec     = (upper_arcmin - upper_arcmin_int) * 60
            self.max_size_label.config(
                text=f"{upper_arcmin_int}′ {upper_arcsec:.2f}″"
            )

        except Exception:
            # on any parse error, reset to 0′ 0.00″
            self.min_size_label.config(text="0′ 0.00″")
            self.max_size_label.config(text="0′ 0.00″")



    LUNATION = 29.53

    def _get_true_illumination(self, lat, lon, when=None):
        """Compute true illuminated fraction from Sun–Moon separation."""
        when = when or datetime.utcnow()
        t = Time(when)
        loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
        m = get_body("moon", t, location=loc)
        s = get_body("sun",  t, location=loc)
        sep = m.separation(s).rad
        return (1 -np.cos(sep)) / 2

    def _get_terminator_angle(self, lat, lon, when=None):
        """Position angle of the dark/light boundary."""
        when = when or datetime.utcnow()
        t = Time(when)
        loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
        m_az = get_body("moon", t, location=loc).transform_to(AltAz(obstime=t, location=loc))
        s_az = get_body("sun",  t, location=loc).transform_to(AltAz(obstime=t, location=loc))
        pa = m_az.position_angle(s_az).to(u.deg).value
        return (pa + 90) % 360

    def _create_moon_image(self, r, frac):
        """Draw a 0°‐oriented crescent for illuminated fraction frac."""
        d = int(2*r)
        img = Image.new('RGB', (d, d), 'black')
        draw = ImageDraw.Draw(img)
        draw.ellipse((0,0,d,d), fill='white', outline='white')
        mask_w = int((1 - frac)*d)
        if frac < 0.5:
            draw.ellipse((0,0,mask_w,d), fill='black')
        else:
            draw.ellipse((d-mask_w,0,d,d), fill='black')
        return img

    def _update_moon(self):
        """Fetch current lat/lon, redraw the moon on its canvas."""
        try:
            lat = float(self.lat_entry.get())
            lon = float(self.lon_entry.get())
        except ValueError:
            return
        frac  = self._get_true_illumination(lat, lon, datetime.utcnow())
        angle = self._get_terminator_angle(lat, lon, datetime.utcnow())
        img   = self._create_moon_image(r=80, frac=frac).rotate(angle, resample=Image.BICUBIC)
        self._tk_moon = ImageTk.PhotoImage(img)
        self.moon_canvas.delete("all")
        self.moon_canvas.create_image(100, 100, image=self._tk_moon)
        self.moon_label.config(text=f"{frac*100:.1f}% illuminated")




    def _on_show_all_objects(self):
            # 1) collect the IDs that passed the last filter
            visible_ids = {obj[0] for obj in self._visible_objects}

            # 2) clear the table
            self.tree.delete(*self.tree.get_children())

            # 3) re-insert every object, greying out the filtered-out ones
            for oid, name, ra, dec, mag, size, max_alt, obj_type, narrowband in self.all_objects:
                values = (
                    name,
                    ra,
                    dec,
                    mag,
                    size,
                    max_alt,
                    obj_type,
                    narrowband,
                    None  # no “visible minutes” for full-list view
                )
                tags = ('greyed',) if oid not in visible_ids else ()
                self.tree.insert("", "end", values=values, tags=tags)

    def _toggle_show_all_objects(self):
        # If checked: show all objects, greying out the ones not in the last filter
        if self.show_all_objects.get():
            visible_ids = {o[0] for o in self._visible_objects}
            self.tree.delete(*self.tree.get_children())
            for oid, name, raw_ra, raw_dec, mag, size, max_alt, typ, nb in self.all_objects:
                try:
                    ra_deg, dec_deg = parse_ra_dec(raw_ra, raw_dec)
                    ra_str = format_ra_hms(ra_deg)
                    dec_str = format_dec_dms(dec_deg)
                except:
                    continue
                tags = ('greyed',) if oid not in visible_ids else ()
                # No “Visible Time” column when showing all
                self.tree.insert(
                    "", "end",
                    values=(name, ra_str, dec_str, mag, size, max_alt, typ, nb, ""),
                    tags=tags
                )
        else:
            # If unchecked: re‐apply the horizon/night filter
            self.objects = list(self._visible_objects)
            self.refresh_table()


    def download_images(self):
        """
        Download (and cache) SkyView cutouts for every object in self.all_objects,
        showing progress in a popup.
        """
        # 1) Build the popup with a progressbar and label
        popup = tk.Toplevel(self.root)
        popup.title("Downloading All Images")
        popup.transient(self.root)
        ttk.Label(popup, text="Downloading images for all deep-sky objects...")\
            .pack(padx=10, pady=(10,5))

        progress = ttk.Progressbar(
            popup, orient="horizontal", length=400, mode="determinate"
        )
        progress.pack(padx=10, pady=5)

        total = len(self.all_objects)
        progress["maximum"] = total

        status_lbl = ttk.Label(popup, text=f"0 / {total}")
        status_lbl.pack(padx=10, pady=(0,10))

        # 2) Helper to update UI from background thread
        def _update(i):
            progress["value"] = i
            status_lbl.config(text=f"{i} / {total}")

        # 3) Worker that runs in your existing ThreadPoolExecutor
        def _worker():
            for idx, rec in enumerate(self.all_objects, start=1):
                # unpack record
                oid, obj_name, raw_ra, raw_dec, *_ = rec
                try:
                    ra, dec = parse_ra_dec(raw_ra, raw_dec)
                    size_arcmin = float(rec[5] or 0)
                    size_deg = max(size_arcmin/60.0 * 1.2, 0.1)

                    # fetch and cache (passes object_name now) :contentReference[oaicite:0]{index=0}&#8203;:contentReference[oaicite:1]{index=1}
                    fetch_survey_image(
                        ra_deg=ra,
                        dec_deg=dec,
                        survey="DSS2 Red",
                        width_deg=size_deg,
                        height_deg=size_deg,
                        object_name=obj_name
                    )
                except Exception as e:
                    print(f"[ERROR] Downloading {obj_name} failed: {e}")

                # schedule UI update (must use after() to cross threads)
                popup.after(0, _update, idx)

            # when done, close the popup
            popup.after(500, popup.destroy)

        # 4) Kick off the background job
        self._executor.submit(_worker)




class ImagePopup(tk.Toplevel):
    def __init__(self, parent, pil_img, title, ra, dec, size_arcmin):
        super().__init__(parent)
        self.parent = parent
        self.original_img = pil_img
        self.zoom = 1.0
        self.ra = ra
        self.dec = dec
        self.size_arcmin = size_arcmin

        self.title(title)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        # 1) Top toolbar: Zoom, Save
        toolbar = ttk.Frame(self)
        toolbar.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        ttk.Button(toolbar, text="–", width=3, command=self._zoom_out).pack(side="left", padx=2)
        ttk.Button(toolbar, text="+", width=3, command=self._zoom_in).pack(side="left", padx=2)
        ttk.Button(toolbar, text="Save As…", command=self._save_image).pack(side="right")

        # 2) Canvas with scrollbars
        canvas = tk.Canvas(self, bg="black")
        canvas.grid(row=1, column=0, sticky="nsew")
        self._canvas = canvas

        vsb = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal", command=canvas.xview)
        vsb.grid(row=1, column=1, sticky="ns")
        hsb.grid(row=2, column=0, sticky="ew")
        canvas.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        # 3) Put the image in a Label inside a Frame on the Canvas
        self._img_frame = ttk.Frame(canvas)
        self._img_id = canvas.create_window((0,0), window=self._img_frame, anchor="nw")
        self._img_label = tk.Label(self._img_frame)
        self._img_label.pack()

        # Render initial
        self._update_display()

        # Make scrolling region follow image size
        self._img_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

    def _update_display(self):
        # 1) Resize original by self.zoom
        w,h = self.original_img.size
        dw, dh = int(w*self.zoom), int(h*self.zoom)
        resized = self.original_img.resize((dw,dh), Image.LANCZOS)

        # 2) Draw crosshair + scale bar
        draw = ImageDraw.Draw(resized)
        # crosshair
        cx, cy = dw//2, dh//2
        draw.line((cx-10,cy, cx+10,cy), fill="yellow")
        draw.line((cx,cy-10, cx,cy+10), fill="yellow")
        # scale bar: assume size_arcmin corresponds to full width
        # e.g. draw 1/5 of width = size_arcmin/5 arcmin
        bar_len = dw//5
        y0 = dh - 20
        x0 = dw//10
        draw.line((x0,y0, x0+bar_len,y0), fill="white", width=2)
        # label
        arcmin_label = f"{self.size_arcmin/5:.1f}′"
        draw.text((x0, y0-15), arcmin_label, fill="white")

        self._tkimg = ImageTk.PhotoImage(resized)
        self._img_label.config(image=self._tkimg)

    def _zoom_in(self):
        self.zoom = min(self.zoom * 1.25, 10)
        self._update_display()

    def _zoom_out(self):
        self.zoom = max(self.zoom / 1.25, 0.1)
        self._update_display()

    def _save_image(self):
        fn = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG","*.png"),("JPEG","*.jpg")]
        )
        if not fn:
            return
        # save the current resized image
        self._tkimg._PhotoImage__photo.write(fn)  # or keep PIL reference and save that

