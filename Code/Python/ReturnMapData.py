

import os
import numpy as np
import astropy.io.fits as pyfits
from astropy.io import ascii
from astropy.table import Table
from astropy.wcs import WCS
from reproject import reproject_interp

import sys
sys.path.append('/home/mayker.1/Desktop/PythonFunctions')
from deprojectGalaxy import deproject


def returnMapData(image, aco, centerCoord, incl, pa):

    # read in fits files
    hdu_int  = pyfits.open(image)
    intMap   = hdu_int[0].data

    #Convert x & y pixels to ra and dec
    wcs      = WCS(hdu_int[0].header, naxis=2)
    naxis    = wcs._naxis # size of image naxis[0] = x and [1] = y
    grid     = np.indices((naxis[1],naxis[0]))
    ra, dec  = wcs.wcs_pix2world(grid[1],grid[0],0)

    #deproject ra and dec to dx and dy
    radius, projang, dx, dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=ra, dec=dec, return_offset=True)

    #flatten data structures 
    f_int  = intMap.flatten()
    f_ra   = ra.flatten()
    f_dec  = dec.flatten()    
    f_dx   = dx.flatten()
    f_dy   = dy.flatten()

    #remove nans
    keep  = np.where(np.isfinite(f_int))
    ra    = f_ra[keep]
    dec   = f_dec[keep]
    inten = f_int[keep]
    dx    = f_dx[keep]
    dy    = f_dy[keep]

    if os.path.isfile(aco):
        hdu_aco = pyfits.open(aco)
        acomap, footprint = reproject_interp(hdu_aco, hdu_int[0].header)
        f_aco = acomap.flatten()
        alphaCO = f_aco[keep]
        alphaCO = np.nan_to_num(alphaCO, nan=6.7)
    else: alphaCO = np.full(len(inten),6.7,dtype=float)
    return(inten, alphaCO, ra, dec, dx, dy)


