"""
Pulls map value at given ra & dec.

Parameters
----------
file : string
path name of map file
ra   : float
right ascension of target in decimal degrees
dec  : float
declination of target in decimal degrees

Returns
-------
value : float
value at given location within map. nan = not in map
"""

def getValue(file, ext, ra, dec):
    
    from IsCoordInMap import IsInMap
    import astropy.io.fits as pyfits
    import os

    if(os.path.isfile(file) == True):

        hdulist = pyfits.open(file)
        map = hdulist[ext].data

        isInMap, xVal, yVal = IsInMap(file, ext, ra, dec)

        if isInMap == True:
            value = map[yVal, xVal]
        else:
            value = float("nan")
    else:
        value = float("nan")

    return(value)
