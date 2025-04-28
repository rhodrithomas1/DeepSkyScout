#!/usr/bin/env python3
"""
Fetch and display a SkyView cutout of M31 (Andromeda Galaxy)
using Astroquery and Pillow.
"""

from astroquery.skyview import SkyView
from astropy import units as u
from astropy.visualization import PercentileInterval, AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from PIL import Image, ImageTk, ImageDraw
import numpy as np
import tkinter as tk

def fetch_survey_image(ra_deg, dec_deg,
                       survey="DSS2 Red",
                       width_deg=0.5,
                       height_deg=0.5):
    """
    Download a small cutout around (ra_deg,dec_deg) from SkyView.
    Returns a PIL Image.
    """
    pos = f"{ra_deg} {dec_deg}"
    imgs = SkyView.get_images(position=pos,
                              survey=[survey],
                              coordinates="J2000",
                              width=width_deg*u.deg,
                              height=height_deg*u.deg)
    if not imgs:
        raise RuntimeError(f"No image returned for {pos} / {survey}")

    hdu = imgs[0][0]  # primary HDU of first survey
    data = hdu.data.astype(float)

    # simple stretch & normalization
    norm = ImageNormalize(data,
                          interval=PercentileInterval(99.5),
                          stretch=AsinhStretch())

    # map to 8-bit
    scaled = (norm(data) * 255).clip(0,255).astype(np.uint8)
    # if grayscale, stack to RGB
    if scaled.ndim == 2:
        rgb = np.stack([scaled]*3, axis=-1)
    else:
        rgb = scaled

    return Image.fromarray(rgb)

def main():
    # M31 J2000 coordinates
    # RA: 00h42m44.3s → 10.684625°
    # Dec: +41°16′9″ → 41.269167°
    ra_m31  = 10.684625
    dec_m31 = 41.269167

    # Fetch a 0.5° × 0.5° cutout around M31
    img = fetch_survey_image(ra_m31, dec_m31,
                             survey="DSS2 Red",
                             width_deg=0.5,
                             height_deg=0.5)

    # Display using Tkinter
    root = tk.Tk()
    root.title("M31 (Andromeda) from DSS2 Red via SkyView")

    photo = ImageTk.PhotoImage(img)
    lbl = tk.Label(root, image=photo)
    lbl.image = photo
    lbl.pack()

    root.mainloop()

if __name__ == "__main__":
    main()
