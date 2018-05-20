# ME398-Magnetic-Mirror

Python functions for an optical setup to scan and create a deformation profile of a 2.5x2.5cm surface. This setup uses a class II laser, a ThorLabs CMOS detector, and two ThorLabs motorized linear stages.

data_fitting.py contains functions to fit the CMOS pixel intensity data to a 2D Gaussian function to find the centroid. The centroid location and optical bench geometry are used to ind the angle and slope of the surface at a given point.

profile_scan.py contains functions to create a 3D surface profile array from position and slope data using Simpson's Rule.

2D Profile Scan.ipynb performs the profile scan and graphs the output.
