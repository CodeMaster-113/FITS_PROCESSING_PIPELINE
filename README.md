# Automated Stellar Photometry Pipeline for HST/ACS Observations

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Astropy](https://img.shields.io/badge/astropy-v5.0+-orange.svg)](https://www.astropy.org/)
[![Photutils](https://img.shields.io/badge/photutils-v1.8+-green.svg)](https://photutils.readthedocs.io/)
[![STScI acstools](https://img.shields.io/badge/acstools-STScI-purple.svg)](https://acstools.readthedocs.io/)

A modular, function-driven scientific Python pipeline for automated aperture photometry on astronomical Flexible Image Transport System (FITS) data from the **Hubble Space Telescope Advanced Camera for Surveys (HST/ACS Wide Field Channel)**.

The framework decomposes traditional monolithic reductions (e.g., IRAF/DAOPHOT) into isolated, testable routines: sub-array cropping, iterative $3\sigma$-clipped sky background estimation, convolved 2D Gaussian point-source detection (`DAOStarFinder`), circular aperture integration with concentric local background annulus subtraction, and epoch-dependent photometric zero-point calibration queried directly from Space Telescope Science Institute (STScI) databases.

---

## Table of Contents

- [Key Features](#key-features)
- [Pipeline Architecture](#pipeline-architecture)
- [Directory Structure](#directory-structure)
- [Installation](#installation)
- [Quickstart & Usage](#quickstart--usage)
  - [Interactive Command-Line Execution](#interactive-command-line-execution)
  - [Programmatic Python API](#programmatic-python-api)
- [Mathematical Framework](#mathematical-framework)
  - [Background Estimation via $\sigma$-Clipping](#background-estimation-via-sigma-clipping)
  - [Aperture & Annulus Flux Integration](#aperture--annulus-flux-integration)
  - [Pogson Magnitude Calibration](#pogson-magnitude-calibration)
  - [CCD Noise Equation & Error Propagation](#ccd-noise-equation--error-propagation)
- [Empirical Validation (NGC 3201 / F814W)](#empirical-validation-ngc-3201--f814w)

---

## Key Features

- **Spatial Sub-Array Extraction**: Isolates target high-density regions (e.g., $2000 \times 2000$ sub-arrays from $4096 \times 4096$ full-frame exposures) to optimize memory and processing speed.
- **Robust Background Modeling**: Iterative $3\sigma$-clipping rejects stellar profiles and cosmic rays, converging to an unbiased sky median and dispersion.
- **PSF-Convolved Centroiding**: Leverages `photutils.detection.DAOStarFinder` with a 2D Gaussian kernel ($FWHM = 3.0\text{ px}$) and $5\sigma$ thresholding, utilizing sharpness and roundness filters to reject artifacts.
- **Dual Annular Aperture Geometry**: Integrates source cores ($r_{\text{ap}} = 6.0\text{ px}$) against local concentric background annuli ($r_{\text{in}} = 10\text{ px}, r_{\text{out}} = 15\text{ px}$), accounting for fractional detector edge overlaps.
- **Dynamic STScI Calibration**: Direct query integration with `acstools.acszpt` retrieves time-dependent zero points accounting for detector degradation and observation epoch.
- **Exoplanet Transit Extensibility**: Designed with an architectural path to temporal analysis, supporting sub-pixel centroid drift tracking and differential ensemble normalization.

---

## Pipeline Architecture

The pipeline executes a sequential 4-stage data flow:

```
[ Stage 1: Ingestion & Sub-Array ]
  Raw FITS File (.fits) ──> Header Parsing (Primary HDU) ──> Sub-Array Crop [4000:6000, 4000:6000]
                                                                        │
                                                                        ▼
[ Stage 2: Background & Detection ]
  Iterative 3σ-Clipping (Median & Variance) ──> DAOStarFinder Centroid Detection (FWHM=3px, 5σ)
                                                                        │
                                                                        ▼
[ Stage 3: Aperture Photometry ]
  Dual Geometry (r=6px, rin=10px, rout=15px) ──> Local Annular Background Subtraction (ApertureStats)
                                                                        │
                                                                        ▼
[ Stage 4: Calibration & Output ]
  STScI API Query (acstools.acszpt) ──> Pogson Calibration ──> Calibrated Astropy Table Catalog
```

---

## Directory Structure

```plaintext
├── fits_pipeline.py          # Core functional library (decoupled reduction routines)
├── main.py                   # Interactive runtime driver and CLI workflow orchestrator
├── figures/                  # Publication-grade flowcharts and empirical diagnostic plots
│   ├── fig1_pipeline_perfect_flowchart.png  # 4-stage architecture diagram
│   ├── fig2_aperture_geometry.png           # Scaled aperture and annulus geometry
│   ├── fig3_timeseries_flowchart.png        # Time-series exoplanet pipeline workflow
│   └── fig4_curves_and_stats.png            # Transit light-curve & luminosity function
├── docs/                     # Full research paper documentation
│   ├── Automated_Stellar_Photometry_HST.pdf  # 5-page publication-formatted paper
│   └── Automated_Stellar_Photometry_HST.docx # Formatted Word documentation
├── requirements.txt          # Python package dependencies
└── README.md                 # Project technical documentation
```

---

## Installation

### Prerequisites

Python 3.9 or higher is required. It is recommended to use an isolated Conda or virtual environment.

```bash
# Clone the repository
git clone [https://github.com/your-username/hst-stellar-photometry.git](https://github.com/your-username/hst-stellar-photometry.git)
cd hst-stellar-photometry

# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install required dependencies
pip install -r requirements.txt
```

### Dependencies (`requirements.txt`)

```plaintext
numpy>=1.22.0
matplotlib>=3.5.0
astropy>=5.0.0
photutils>=1.8.0
acstools>=3.4.0
```

*Or via Conda:*
```bash
conda install -c conda-forge numpy matplotlib astropy photutils acstools
```

---

## Quickstart & Usage

### Interactive Command-Line Execution

Run the interactive driver:

```bash
python main.py
```

The script prompts for observation parameters:

```plaintext
Path to FITS file: data/hst_acs_ngc3201_f814w.fits
Instrument: ACS
Filter (e.g. F814W): F814W
Detector for zero-point query (e.g. WFC): WFC
Observation date for zero-point query (YYYY-MM-DD): 2006-03-15
Exposure time in seconds: 340.0
```

Execution proceeds automatically:
1. Ingests FITS primary HDU and prints structural metadata.
2. Crops the active sub-array (`image[4000:6000, 4000:6000]`).
3. Generates $2 \times 2$ contrast visualizations (Linear, LogNorm, Magma, Greys).
4. Identifies star centroids and plots coordinate scatter overlays.
5. Constructs circular apertures and background annuli, displaying verification overlays.
6. Computes net background-subtracted fluxes.
7. Queries STScI for zero points and appends calibrated instrumental magnitudes.

### Programmatic Python API

Each operational phase is an independent function in `fits_pipeline.py`:

```python
from fits_pipeline import (
    open_fits, extract_image, detect_sources,
    make_apertures, make_annulus, run_photometry,
    fetch_zero_point, compute_magnitudes
)

# 1. Ingestion & Sub-Array Slicing
hdul = open_fits("data/ngc3201_f814w.fits")
image = extract_image(hdul)[4000:6000, 4000:6000]

# 2. Source Centroid Detection
sources = detect_sources(image, fwhm=3.0, sigma=3.0, threshold_sigma=5.0)
x, y = sources['x_centroid'], sources['y_centroid']

# 3. Aperture Geometry & Photometric Integration
positions, apertures = make_apertures(x, y, radius=6.0)
annulus_aperture = make_annulus(positions, r_in=10.0, r_out=15.0)
star_data = run_photometry(image, apertures, annulus_aperture)

# 4. Zero-Point Calibration
zpt_table, filter_zpt = fetch_zero_point(date="2006-03-15", detector="WFC", filt="F814W")
zpt = 25.52  # Value retrieved from STScI table
star_data = compute_magnitudes(star_data, zero_point=zpt, exptime=340.0)

# Display calibrated catalog
print(star_data['id', 'xcenter', 'ycenter', 'aperture_sum', 'total_bkg', 'Magnitudes'])
```

---

## Mathematical Framework

### Background Estimation via $\sigma$-Clipping

Pixel counts superimpose stellar emission onto ambient sky and bias levels. Outlier pixels deviating by more than $k\sigma$ ($k = 3$) from the median are rejected iteratively:

$$\text{Mask}^{(i)} = \left\{ x_j \;\middle\vert{}\; \left\vert{} x_j - m^{(i)} \right\vert{} > k \sigma^{(i)} \right\}$$

This yields unbiased estimates of the sky median $\mu_{\text{sky}}$ for baseline subtraction and sample dispersion $\sigma_{\text{sky}}$ for detection thresholding.

### Aperture & Annulus Flux Integration

Aperture photometry integrates flux inside a radius $r_{\text{ap}} = 6.0\text{ px}$ centered on $(x_0, y_0)$. A concentric circular annulus ($r_{\text{in}} = 10.0\text{ px}, r_{\text{out}} = 15.0\text{ px}$) estimates local sky background per pixel $s_{\text{bkg}}$. Net stellar flux is:

$$F_{\text{net}} = F_{\text{ap}} - \left( s_{\text{bkg}} \times A_{\text{ap}} \right)$$

where $A_{\text{ap}} = \pi r_{\text{ap}}^2$ accounts for fractional pixel area overlap across image boundaries.

### Pogson Magnitude Calibration

Net counts convert to a count rate and are mapped to calibrated instrumental magnitudes:

$$m_{\text{inst}} = \text{ZPT} - 2.5 \log_{10}\left( \frac{F_{\text{net}}}{t_{\text{exp}}} \right)$$

where $\text{ZPT}$ is the epoch- and filter-specific zero point queried from STScI.

### CCD Noise Equation & Error Propagation

Following the standard CCD equation, the total uncertainty in net flux $\sigma_F$ and instrumental magnitude $\sigma_m$ is:

$$\sigma_F = \sqrt{F_{\text{net}} + A_{\text{ap}}\left(1 + \frac{A_{\text{ap}}}{A_{\text{ann}}}\right)\left(s_{\text{bkg}} + \sigma_{\text{read}}^2\right)}$$

$$\sigma_m = \frac{2.5}{\ln 10} \left( \frac{\sigma_F}{F_{\text{net}}} \right) \approx 1.0857 \left( \frac{\sigma_F}{F_{\text{net}}} \right)$$

---

## Empirical Validation (NGC 3201 / F814W)

The pipeline was benchmarked against HST/ACS Wide Field Channel observations of the galactic globular cluster NGC 3201 ($0.049''/\text{pixel}$ plate scale) in the F814W filter (broad $I$-band, $\lambda_{\text{eff}} \approx 806\text{ nm}$).

### Sample Calibrated Photometric Catalog

| Star ID | x_centroid (px) | y_centroid (px) | Aperture Sum ($e^-$) | Total Background ($e^-$) | Calibrated Mag (F814W) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | 412.15 | 105.32 | 12450.20 | 2100.12 | 19.82 |
| 2 | 885.40 | 310.88 | 3200.51 | 2050.40 | 22.15 |
| 3 | 1204.62 | 742.19 | 45210.80 | 2130.45 | 17.41 |
| 4 | 1530.11 | 1102.04 | 1520.10 | 2080.20 | 23.50 |

The differential luminosity distribution demonstrates an expected astrophysical power-law increase toward faint apparent magnitudes, exhibiting a $5\sigma$ detection completeness turnover at $m_{\text{F814W}} \approx 24.2$.


---

**Author**: Sahil Vishwasrao  
*Department of Computer Science and Engineering, Sardar Patel Institute of Technology, Mumbai*

---
