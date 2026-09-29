"""
FITS star detection and aperture photometry pipeline.
Converted from fits.ipynb -- every notebook step is now a function.
Run:  python fits_pipeline.py
"""

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.aperture import (
    CircularAperture,
    CircularAnnulus,
    aperture_photometry,
    ApertureStats,
)
from acstools import acszpt


# ---------------------------------------------------------------------------
# Cell 0 : Open the FITS file
# ---------------------------------------------------------------------------
def open_fits(path):
    """Open a FITS file, print its info and return the HDUList."""
    hdul = fits.open(path)
    print(hdul.info())
    return hdul


# ---------------------------------------------------------------------------
# Cell 1 : Extract the channels
# ---------------------------------------------------------------------------
def extract_image(hdul):
    """Return the image data from the primary HDU."""
    image = hdul[0].data
    print(image)
    return image


# ---------------------------------------------------------------------------
# Cell 2 : Image plot
# ---------------------------------------------------------------------------
def plot_image(image):
    """Show the four plots from the notebook in a single frame (2x2):
    linear, log, log+magma, log+Greys."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    im = axes[0, 0].imshow(image, origin='lower')
    fig.colorbar(im, ax=axes[0, 0])

    im = axes[0, 1].imshow(image, origin='lower', norm=LogNorm())
    fig.colorbar(im, ax=axes[0, 1])

    im = axes[1, 0].imshow(image, origin='lower', norm=LogNorm(), cmap='magma')
    fig.colorbar(im, ax=axes[1, 0])

    im = axes[1, 1].imshow(image, origin='lower', norm=LogNorm(), cmap='Greys')
    fig.colorbar(im, ax=axes[1, 1])

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Cell 3 : Stars detection
# ---------------------------------------------------------------------------
def detect_sources(image, fwhm=3.0, sigma=3.0, threshold_sigma=5.0):
    """
    1. Estimate background mean, median and std with sigma clipping.
    2. Find stars at least `threshold_sigma` * std above background (DAOStarFinder).
    Returns the sources table.
    """
    mean, median, std = sigma_clipped_stats(image, sigma=sigma)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
    sources = daofind(image - median)
    return sources


# ---------------------------------------------------------------------------
# Cell 4 : Plot detected stars
# ---------------------------------------------------------------------------
def plot_detections(image, sources):
    """Overplot detected star centroids on the image and print the count."""
    x = sources['x_centroid']
    y = sources['y_centroid']

    plt.figure()
    plt.imshow(image, origin='lower', norm=LogNorm(), cmap='Greys')
    plt.scatter(x, y)
    plt.colorbar()
    plt.show()

    print(f"{len(sources)} radiation emmitting bodies found")
    return x, y


# ---------------------------------------------------------------------------
# Cell 5 : Aperture circling
# ---------------------------------------------------------------------------
def make_apertures(x, y, radius=6.0):
    """Create circular apertures at the detected star positions."""
    positions = np.transpose((x, y))
    apertures = CircularAperture(positions, r=radius)
    return positions, apertures


def plot_apertures(image, apertures):
    """Plot the apertures as hollow cyan circles."""
    plt.figure()
    plt.imshow(image, origin='lower', norm=LogNorm(), cmap='Greys')
    apertures.plot(color='cyan', lw=1.5, alpha=0.9)
    plt.colorbar()
    plt.show()


# ---------------------------------------------------------------------------
# Cell 6 : Aperture photometry
# ---------------------------------------------------------------------------
def make_annulus(positions, r_in=10, r_out=15):
    """Create the background annulus around each star."""
    return CircularAnnulus(positions, r_in=r_in, r_out=r_out)


def plot_aperture_and_annulus(image, apertures, annulus_aperture):
    """Plot apertures (blue) and annuli (green) on separate figures."""
    plt.figure()
    plt.imshow(image, cmap='Greys', origin='lower', norm=LogNorm())
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()

    plt.figure()
    plt.imshow(image, cmap='Greys', origin='lower', norm=LogNorm())
    annulus_aperture.plot(color='green', lw=1.5, alpha=0.5)
    plt.show()


def run_photometry(image, apertures, annulus_aperture):
    """Compute aperture sums and total background; return the star_data table."""
    aperture_stats = ApertureStats(image, annulus_aperture)
    bkg_mean = aperture_stats.mean
    aperture_area = apertures.area_overlap(image)

    total_bkg = bkg_mean * aperture_area
    star_data = aperture_photometry(image, apertures)
    star_data['total_bkg'] = total_bkg

    for col in star_data.colnames:
        star_data[col].info.format = '%.8g'

    star_data.pprint()
    return star_data


# ---------------------------------------------------------------------------
# Cell 7 : Fetch the zero-point
# ---------------------------------------------------------------------------
def fetch_zero_point(date, detector, filt):
    """
    Query the ACS zero-point calculator.
    Returns (zpt_table for all filters, filter_zpt for the given filter).
    """
    q = acszpt.Query(date=date, detector=detector)
    zpt_table = q.fetch()

    q_filter = acszpt.Query(date=date, detector=detector, filt=filt)
    filter_zpt = q_filter.fetch()

    print(filter_zpt)
    return zpt_table, filter_zpt


# ---------------------------------------------------------------------------
# Cell 8 : Magnitude of each detected star
# ---------------------------------------------------------------------------
def compute_magnitudes(star_data, zero_point, exptime):
    """
    mag = zero_point - 2.5 * log10(|aperture_sum - total_bkg| / exptime)
    Adds a 'Magnitudes' column to star_data and returns the list.
    """
    magnitudes = []
    for line in star_data:
        numerator = abs(line['aperture_sum'] - line['total_bkg'])
        mag = zero_point - 2.5 * (math.log10(numerator / exptime))
        magnitudes.append(mag)

    star_data['Magnitudes'] = magnitudes
    star_data.pprint()
    return magnitudes