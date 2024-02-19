
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS


#Checks if Supernovae are in an image and, if in the map, reports back their name, type, ra, dec, map, and, x & y coords
def check_in_image(SNras, SNdecs, SNnames, SNtypes, image, ext):
    '''Checks if a SN is within the footprint of a map and returns the name, type, ra, dec, map, and x & y coords if so'''
    intensity = []
    hdulist  = pyfits.open(image)
    
    map = hdulist[ext].data
    wcs = WCS(hdulist[ext].header, naxis=2)
    naxis = wcs._naxis #size of image
    naxis1 = hdulist[ext].header['NAXIS1']
    naxis2 = hdulist[ext].header['NAXIS2']

    coords_arr = np.column_stack((SNras, SNdecs)) # ras and decs now [ra,dec]
    pix_x, pix_y = wcs.wcs_world2pix(SNras,SNdecs,0)

    #use world coordinates of all SNe to see if any fall in image (our version of footprint_contains)
    is_in_x = (pix_x >= 0) & (pix_x <= naxis[0]-1) #because of 0-indexing
    is_in_y = (pix_y >= 0) & (pix_y <= naxis[1]-1)
             
    #get the name, ra, and dec of the SNe that fall in image
    #boolean array indexing (gives back array of Trues and Falses)
    #we are pulling out the SNe that are True and assigning them to own arrays
    name_in_image = np.array(SNnames)[is_in_x & is_in_y]
    type_in_image = np.array(SNtypes)[is_in_x & is_in_y]
    ra_in_image = np.array(SNras)[is_in_x & is_in_y]
    dec_in_image = np.array(SNdecs)[is_in_x & is_in_y]
    
    x_coord = np.array(pix_x)[is_in_x & is_in_y]
    y_coord = np.array(pix_y)[is_in_x & is_in_y]
    
    return(name_in_image, type_in_image, ra_in_image, dec_in_image, map, x_coord, y_coord)
