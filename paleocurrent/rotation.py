import numpy as np

## Define functions to apply a tilt correction to paleocurrents

# These functions are from the python library ```pmagpy``` and are defined here so that 
# library does not need to be imported.

def dir2cart(d):
    """
    Converts a list or array of vector directions in degrees (declination,
    inclination) to an array of the direction in cartesian coordinates (x,y,z)

    Parameters
    ----------
    d : list or array of [dec,inc] or [dec,inc,intensity]

    Returns
    -------
    cart : array of [x,y,z]

    Examples
    --------
    >>> pmag.dir2cart([200,40,1])
    array([-0.71984631, -0.26200263,  0.64278761])
    """
    ints = np.ones(len(d)).transpose(
    )  # get an array of ones to plug into dec,inc pairs
    d = np.array(d).astype('float')
    rad = np.pi/180.
    if len(d.shape) > 1:  # array of vectors
        decs, incs = d[:, 0] * rad, d[:, 1] * rad
        if d.shape[1] == 3:
            ints = d[:, 2]  # take the given lengths
    else:  # single vector
        decs, incs = np.array(float(d[0])) * rad, np.array(float(d[1])) * rad
        if len(d) == 3:
            ints = np.array(d[2])
        else:
            ints = np.array([1.])
    cart = np.array([ints * np.cos(decs) * np.cos(incs), ints *
                     np.sin(decs) * np.cos(incs), ints * np.sin(incs)])
    cart = np.array([ints * np.cos(decs) * np.cos(incs), ints *
                     np.sin(decs) * np.cos(incs), ints * np.sin(incs)]).transpose()
    return cart

def dotilt(dec, inc, bed_az, bed_dip):
    """
    Does a tilt correction on a direction (dec,inc) using bedding dip direction
    and bedding dip.

    Parameters
    ----------
    dec : declination directions in degrees
    inc : inclination direction in degrees
    bed_az : bedding dip direction
    bed_dip : bedding dip

    Returns
    -------
    dec,inc : a tuple of rotated dec, inc values

    Examples
    -------
    >>> dotilt(91.2,43.1,90.0,20.0)
    (90.952568837153436, 23.103411670066617)
    """
    rad = np.pi / 180.  # converts from degrees to radians
    X = dir2cart([dec, inc, 1.])  # get cartesian coordinates of dec,inc
# get some sines and cosines of new coordinate system
    sa, ca = -np.sin(bed_az * rad), np.cos(bed_az * rad)
    cdp, sdp = np.cos(bed_dip * rad), np.sin(bed_dip * rad)
# do the rotation
    xc = X[0] * (sa * sa + ca * ca * cdp) + X[1] * \
        (ca * sa * (1. - cdp)) + X[2] * sdp * ca
    yc = X[0] * ca * sa * (1. - cdp) + X[1] * \
        (ca * ca + sa * sa * cdp) - X[2] * sa * sdp
    zc = X[0] * ca * sdp - X[1] * sdp * sa - X[2] * cdp
# convert back to direction:
    Dir = cart2dir([xc, yc, -zc])
    # return declination, inclination of rotated direction
    return Dir[0], Dir[1]

def cart2dir(cart):
    """
    Converts a direction in cartesian coordinates into declination, inclinations

    Parameters
    ----------
    cart : input list of [x,y,z] or list of lists [[x1,y1,z1],[x2,y2,z2]...]

    Returns
    -------
    direction_array : returns an array of [declination, inclination, intensity]

    Examples
    --------
    >>> pmag.cart2dir([0,1,0])
    array([ 90.,   0.,   1.])
    """
    cart = np.array(cart)
    rad = np.pi/180.  # constant to convert degrees to radians
    if len(cart.shape) > 1:
        Xs, Ys, Zs = cart[:, 0], cart[:, 1], cart[:, 2]
    else:  # single vector
        Xs, Ys, Zs = cart[0], cart[1], cart[2]
    if np.iscomplexobj(Xs):
        Xs = Xs.real
    if np.iscomplexobj(Ys):
        Ys = Ys.real
    if np.iscomplexobj(Zs):
        Zs = Zs.real
    Rs = np.sqrt(Xs**2 + Ys**2 + Zs**2)  # calculate resultant vector length
    # calculate declination taking care of correct quadrants (arctan2) and
    # making modulo 360.
    Decs = (np.arctan2(Ys, Xs)/rad) % 360.
    try:
        # calculate inclination (converting to degrees) #
        Incs = np.arcsin(Zs/Rs)/rad
    except:
        print('trouble in cart2dir')  # most likely division by zero somewhere
        return np.zeros(3)

    return np.array([Decs, Incs, Rs]).transpose()  # return the directions list

def doflip(dec, inc):
    """
    flips lower hemisphere data to upper hemisphere
    """
    if inc < 0:
        inc = -inc
        dec = (dec + 180.) % 360.
    return dec, inc