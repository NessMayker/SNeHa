def IsInMap(file, ext, ra, dec):
    """
    Determines if an ra and dec is contained within the coverage of a map.    

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
    isInMap : bool
        boolean if target is in map
    xVal    : int 
        x coordinate of target in map ("nan" = not in map)
    yVal    : int
        y coordinate of target in map ("nan" = not in map)
    """
    import numpy as np
    from astropy.wcs import WCS
    import astropy.io.fits as pyfits

    hdulist = pyfits.open(file)
    map = hdulist[ext].data
    wcs = WCS(hdulist[ext].header, naxis=2)
    
    pix_x, pix_y = wcs.wcs_world2pix(ra, dec, 0, ra_dec_order=True)

    naxis = wcs._naxis #size of image
    naxis1 = hdulist[ext].header['NAXIS1']
    naxis2 = hdulist[ext].header['NAXIS2']
    is_in_x = (pix_x >= 0) & (pix_x <= naxis[0]-1) #because of 0-indexing
    is_in_y = (pix_y >= 0) & (pix_y <= naxis[1]-1)
    
    x_coord = np.array(pix_x)[is_in_x & is_in_y]
    y_coord = np.array(pix_y)[is_in_x & is_in_y]
    
    isInMap = False
    
    if(is_in_x == True & is_in_y == True):
        isInMap = True
        xVal = int(round(float(x_coord)))
        yVal = int(round(float(y_coord)))
        
    else:
        xVal = float("nan")
        yVal = float("nan") 
    
    return(isInMap, xVal, yVal)