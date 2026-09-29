# Automated Stellar Photometry Pipeline for HST/ACS Observations

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Astropy](https://img.shields.io/badge/astropy-v5.0+-orange.svg)](https://www.astropy.org/)
[![Photutils](https://img.shields.io/badge/photutils-v1.8+-green.svg)](https://photutils.readthedocs.io/)
[![STScI acstools](https://img.shields.io/badge/acstools-STScI-purple.svg)](https://acstools.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
- [Time-Series Extension (Exoplanet Transits)](#time-series-extension-exoplanet-transits)
- [Limitations & Roadmap](#limitations--roadmap)
- [Citation & Contact](#citation--contact)
- [License](#license)

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
