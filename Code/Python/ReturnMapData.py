
# import packages
import os
import numpy as np
import astropy.io.fits as pyfits
from astropy.io import ascii
from astropy.table import Table
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.coordinates import SkyCoord, Angle
from astropy.nddata import Cutout2D
import astropy.units as u

# import local packages
import sys
sys.path.append('/users/nessmaykerchen/Desktop/Research/SNeHa/Code/Python')
from deprojectGalaxy import deproject
from AngularSize import findAngSize


def returnMapData(image, ext, centerCoord, incl, pa, dist, SNra = None, SNdec = None, local = None):

    # read in fits files
    hdu_int  = pyfits.open(image)
    intMap   = hdu_int[ext].data

    #Convert x & y pixels to ra and dec
    wcs      = WCS(hdu_int[ext].header, naxis=2)
    naxis    = wcs._naxis # size of image naxis[0] = x and [1] = y
    grid     = np.indices((naxis[1],naxis[0]))
    ra, dec  = wcs.wcs_pix2world(grid[1],grid[0],0)

    #deproject ra and dec to dx and dy
    radius, projang, dx, dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=ra, dec=dec, return_offset=True)

    # isolate local environment if specified
    if local != None:
        angBoxSize = findAngSize(local, dist) # 500 pc in decimal degrees

        #cutout2D needs skycoord position to carry units
        ra_cut, dec_cut = Angle(SNra * u.degree), Angle(SNdec * u.degree)
        raRad, decRad  = ra_cut.radian * u.rad, dec_cut.radian * u.rad    
        position = SkyCoord(raRad, decRad) #position is center, use ra & dec of SN location
        size = u.Quantity((angBoxSize,angBoxSize), u.degree) #size is size of box in arcsec 
    
        # make 2D cutout, will assign a new wcs to cutout to keep track of coords
        cutout_int = Cutout2D(intMap, position, size, wcs) 
        f_int   = cutout_int.data.flatten()
        cutout_ra = Cutout2D(ra, position, size, wcs) 
        f_ra   = cutout_ra.data.flatten()
        cutout_dec = Cutout2D(dec, position, size, wcs) 
        f_dec   = cutout_dec.data.flatten()
        cutout_dx = Cutout2D(dx, position, size, wcs) 
        f_dx   = cutout_dx.data.flatten()
        cutout_dy = Cutout2D(dy, position, size, wcs) 
        f_dy   = cutout_dy.data.flatten()

    else: 
        f_int  = intMap.flatten()
        f_ra   = ra.flatten()
        f_dec  = dec.flatten()    
        f_dx   = dx.flatten()
        f_dy   = dy.flatten()

    keep  = np.where(np.isfinite(f_int))
    inten = f_int[keep]
    ra    = f_ra[keep]
    dec   = f_dec[keep]
    dx    = f_dx[keep]
    dy    = f_dy[keep]

    return(inten, ra, dec, dx, dy)


