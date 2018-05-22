import numpy as np
from scipy.integrate import simps       #Simpsons rule

# Test functions
def h(x,y): return (x/1.)**2+(y/2.)**2
def dhx(x,y): return (2*x**1.)
def dhy(x,y): return (0.5*y**1.)

def test2d(xlim=[0,5],ylim=[0,6],num=[6,7]):
    """" Tests the Simpson's rule numerical integration method using the test functions h(x,y), dhx(x,y) (partial
    derivative with respect to x), and dhy(x,y) (partial derivative with respect to y). Prints the slopes in the x
    direction, the slopes in the y direction.

    Args:
        xlim: [lower, uppper] numpy array with x upper and lower bounds
        ylim: [lower, upper] numpy array with y lower and upper bounds
        num: [x_elements, y_elements] numpy array with the total number of elements in x and y arrays

    Prints:
        x,y: x and y arrays
        hx2d, hy2d: x and y slope arrays
        Ixy1: Surface profile array from simps
        Ixy==Ixy1: 1: If the surface profile computed from simps is euqal to the surface profile computed in profile2d
    """
    x, y = [np.linspace(xlim[0],xlim[1],num[0]), np.linspace(ylim[0],ylim[1],num[1])]
    print (x)
    print (y)
    hx, hy = dhx(x,y), dhy(x,y)
    hx2d = np.array([hx for j in range(len(hy))])
    hy2d = np.transpose([hy for i in range(len(hx))])
    print (hx2d)
    print (hy2d)
    Ix0 = [simps(hx[0:i],x[0:i]) for i in range(1,num[0]+1)]
    Iy0 = [simps(hy[0:j],y[0:j]) for j in range(1,num[1]+1)]
    Ixy = np.transpose([[simps(hy[0:j],y[0:j]) for j in range(1,num[1]+1)]+xi for xi in Ix0])
    Ixy = Ixy #- np.min(Ixy)
    Ixy1 = profile2d(hx2d,hy2d,x,y)
    print (Ixy1)
    print (Ixy==Ixy1)


def profile2d(hx2d,hy2d,x,y):
    """"" Creates a 3D surface profile array from position and slope data using Simpson's rule.

    Args:
        hx2d: Slopes array relative to the xy plane.
        hy2d: Slopes array relative to the xz plane

    Prints:
        Ix: Surface profile array with respect to x. Indices in the array represent xy positions. Values in the array
            represent heights in z.
        Iy: Surface profile array with respect to y. Indices in the array represent xy positions. Values in the array
            represent heights in z.

    Returns:
        Ixy: 2D surface profile array. Indices in the array represent xy positions. Values in the array represent
            heights in z.
    """
    if hx2d.shape == hy2d.shape: num = [len(hx2d[0]),len(hx2d)]
    Ix = np.array([[simps(hx2d[j][0:i],x[0:i]) for i in range(1,num[0]+1)] for j in range(0,num[1])])
    print (Ix)
    Iy = np.transpose([[simps(np.transpose(hy2d)[i][0:j],y[0:j]) for j in range(1,num[1]+1)] for i in range(0,num[0])])
    print (Iy)
    Ixy = Ix+Iy #- np.min(Ix+Iy)
    return Ixy