"""
Main script: runs the FITS photometry pipeline in notebook order.
Requires fits_pipeline.py in the same folder.
Run:  python main.py
"""

from fits_pipeline import (
    open_fits,
    extract_image,
    plot_image,
    detect_sources,
    plot_detections,
    make_apertures,
    plot_apertures,
    make_annulus,
    plot_aperture_and_annulus,
    run_photometry,
    fetch_zero_point,
    compute_magnitudes,
)


def main():
    # ---- user inputs ----
    path = input("Path to FITS file: ").strip()
    instrument = input("Instrument: ").strip()
    fil = input("Filter (e.g. F814W): ").strip()
    detector = input("Detector for zero-point query (e.g. WFC): ").strip()
    date = input("Observation date for zero-point query (YYYY-MM-DD): ").strip()
    exptime = float(input("Exposure time in seconds: ").strip())

    print(f"\nInstrument: {instrument} | Filter: {fil} | Detector: {detector} "
          f"| Date: {date} | Exposure: {exptime}s\n")

    # Cell 0-1 : open file and extract image
    hdul = open_fits(path)
    image = extract_image(hdul)
    image = image[4000:6000,4000:6000]

    # Cell 2 : image plots
    plot_image(image)

    # Cell 3-4 : detect sources and plot them
    sources = detect_sources(image)
    x, y = plot_detections(image, sources)

    # Cell 5 : apertures
    positions, apertures = make_apertures(x, y)
    plot_apertures(image, apertures)

    # Cell 6 : background annulus and photometry
    annulus_aperture = make_annulus(positions)
    plot_aperture_and_annulus(image, apertures, annulus_aperture)
    star_data = run_photometry(image, apertures, annulus_aperture)

    # Cell 7 : zero point
    fetch_zero_point(date, detector, fil)
    zero_point = float(input("Zero point (read it from the table printed above): ").strip())

    # Cell 8 : magnitudes
    print(f"\nApparent magnitude and source data in {fil} band:")
    compute_magnitudes(star_data, zero_point, exptime)

    hdul.close()


if __name__ == "__main__":
    main()