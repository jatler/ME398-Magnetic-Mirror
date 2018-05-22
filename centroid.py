import scipy.optimize as opt
import numpy as np
import pylab as plt
import collections

def twoD_Gaussian(X, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """" Computes the value of a 2D Gaussian function at the point X = (x,y)
    
    Args:
        X: Position array (x,y) used to compute the value of the 2D gaussian function
        amplitude: Amplitude of general 2D gaussian function
        x0,y0: X- and y- coordinates of the center of general 2D gaussian function
        sigma_x, sigma_y: x and y standard deviations
        theta:  rotation angle of the 2D gaussian function
        offset: asymptote distance from y = 0
        
    Returns:
        f.ravel(): 1D array of general 2D gaussian function output at (x,y).
    """"

    x,y = X
    (xo,yo) = (float(xo),float(yo))

    # a,b,c: parameters for the general 2D gaussian function
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)

    #f: equation for the 2D gaussian function
    f = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                                      + c*((y-yo)**2)))
    return f.ravel()

def centroid(data,plot=0,A=200,x0=640,y0=512,sx=180,sy=150,theta=0,h=5,run=0):
    """ Fits a 2D gaussian function to the pixel intensity data received from the CMOS. The fit uses scipy.optimize.curve_fit, which uses non-linear least squares to ft the 2D
    gaussian function to the data.

    Args:
        data: A 1280 x 1024 array of CMOS pixel intensity values ranging from 0 to 255
        plot: default to 0 - does not plot, 1 - plots the pixel data with contours of the 2D gaussian fit
        A: Initial guess for the 2D gaussian amplitude.
        x0,y0: Initial guess for the 2D gaussian center coordinates
        sx, sy: Initial guess for the 2D guassian x and y standard deviations
        theta: Initial guess for the 2D gaussian rotation angle
        h: Initial guess for the 2D gaussian asymptote offset from y = 0
        run: run number used in the file name of saved figures

    Returns:
        popt: Fitted parameters of the 2D gaussian function
        perr: Fitted 2D gaussian standard deviation values
    """

    initial_guess = (A,x0,y0,sx,sy,theta,h)
    data = data.ravel()
    x = np.linspace(0, 1279, 1280)
    y = np.linspace(0, 1023, 1024)
    x, y = np.meshgrid(x, y)

    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data, p0=initial_guess)
    perr = np.sqrt(np.diag(pcov))
    
    data_fitted = twoD_Gaussian((x, y), *popt)
    if plot == 1:
        fig, ax = plt.subplots(1, 1)
        ax.imshow(data.reshape(1024, 1280), cmap='gist_gray', origin='bottom',
          extent=(x.min(), x.max(), y.min(), y.max()))
        ax.contour(x, y, data_fitted.reshape(1024,1280), 8, colors='w', linewidths=0.3)
        """ Update directory location when saving figures """
        fig.savefig('C:/Users/Johnny Atler/Desktop/centroid_%s.pdf'%(run)); plt.close()

    return [popt,perr]

def angles(data,d,w):
    """ Computes the x and y angles of the surface from a 1280 x 1024 array of CMOS light intensity and two parameters
    of the optical bench layout. A 2D gaussian function is fit to the data to find the centroid. The centroid location
    and optical bench geometry are used to find the angle of the surface.

    Args:
        data: 1280 x 1024 array of CMOS pixel intensity values ranging from 0 to 255
        d: Parallel distance from the CMOS sensor face (0,0) pixel to the surface
        w: Perpandicular distance from the CMOS sensor face (0,0) pixel to the surface

    Returns:
        angle_xa: Angle in arcseconds of the surface relative to the xz plane
        angle_ya: Angle in arcseconds of the surface relative to the yz plane
    """

    # Centroid image (i.e. Fit Them with 2D Gaussians)
    popt, perr = centroid(data)

    size_pixel = 5.2*10**-6 # size of CMOS pixel in m
    # Calculate angle
    x_cmos = (1280 - popt[1])*size_pixel  # facing the sensor, x pixel = 0 on left, 1280 on right
    y_cmos = (1024 - popt[2])*size_pixel  # facing the sensor, y pixel = 0 on bottom, 1024 on top
    angle_x = 0.5*(np.arctan((w+x_cmos)/d)
    angle_y = 0.5*(np.arctan(y_cmos/d))

    angle_xa = angle_x*180/np.pi*3600  #angle in arcseconds
    angle_ya = angle_y*180/np.pi*3600  #angle in arcseconds
    return angle_xa, angle_ya

def slopes(data, d = 2.235, w = 0.914):
    """ Computes the slope of the surface relative to the xz plane from a 1280x1024 array of CMOS light intensity and
    two parameters of the optical bench layout. A 2D gaussian function is fit to the data to find the centroid. The
    centroid location and optical bench geometry are used to find the angle of the surface, which is used to find the
    slope.

    Args:
        data: 1280 x 1024 array of CMOS pixel intensity values ranging from 0 to 255
        d: Parallel distance from the CMOS sensor face (0,0) pixel to the surface
        w: Perpandicular distance from the CMOS sensor face (0,0) pixel to the surface

    Returns:
        slope_x: Slope of the surface relative to the xz plane
        slope_y: Slope of the surface relative to the yz plane
    """

    angle_xa, angle_ya = angles(data, d, w)
    angle_x = angle_xa*np.pi*3600/180
    angle_y = angle_ya * np.pi * 3600 / 180
    slope_x = np.tan(angle_x)
    slope_y = np.tan(angle_y)

    return slope_x, slope_y