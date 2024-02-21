

# import global packages
import os
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from reproject import reproject_interp

# import local packages
import sys
sys.path.append('/users/nessmaykerchen/Desktop/Research/SNeHa/Code/Python')
from deprojectGalaxy import deproject
from AngularSize import findAngSize
from AngDistToPc import angDistToPc
from ReturnMapData import returnMapData
from nearest import findNearest
from getMapValue import getValue


# runs through each galaxy file and finds nearest HII region to each SN site 
# (global, for local set local = resolution).
def nearestHII(galaxy, HII_file, ext, centerCoord, pa, incl, galDist, SNras, SNdecs, local=None):
    '''takes galaxy information and returns nearest HII region (in pc) and its coordinates'''
    inHII = []
    if os.path.isfile(HII_file):
        hdu  = pyfits.open(HII_file)
        wcs  = WCS(hdu[ext].header, naxis=2)

        HII_val, nearestHII, HIIras, HIIdecs = [],[],[],[]

        # iterate through each SN environment and find nearest HII region

            
        for i in range(len(SNras)):
            snra = float(SNras[i])
            sndec = float(SNdecs[i])

            # returns local map (with cutout centered on SN if local!=None)
            HII, ra, dec, dx, dy = returnMapData(HII_file, ext, centerCoord, incl, pa, galDist, snra, sndec, local)

            # Get the value from the HII map at the SN location
            HIIval = getValue(HII_file,ext,snra,sndec)
            HII_val.append(HIIval)

            #remove nans
            keep  = np.where(HII >= 0)
            map_ra    = ra[keep]
            map_dec   = dec[keep]
            map_HII   = HII[keep]
            map_dx    = dx[keep]
            map_dy    = dy[keep]

            # get deprojected SN coordinates
            SN_Rad, SN_PA, SN_dx, SN_dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=snra, dec=sndec, return_offset=True) 

            if len(map_dx) != 0:

                # find nearest HII region and convert distance to pc
                nearestHII_ang, idx = findNearest(map_dx, SN_dx, map_dy, SN_dy)
                HIIra, HIIdec = map_ra[idx], map_dec[idx]
                nearestHII_pc = angDistToPc(nearestHII_ang, galDist)

                nearestHII.append(nearestHII_pc)
                HIIras.append(HIIra) 
                HIIdecs.append(HIIdec)
            else:pass

        # Check if the SN is in an HII region and adjust nearest HII region value accordingly
        # (want a distance value = 0 if the SN is IN the region)
        for i in range(len(SNras)):
            if HII_val[i] > 0:
                inHII.append("Yes")
                nearestHII[i] = 0.
            else:
                inHII.append("No")
    
    else:
        print("No file for ", galaxy)

        nearestHII = float("nan")
        HIIras = float("nan")
        HIIdecs = float("nan")
        
    return(inHII, nearestHII, HIIras, HIIdecs, HII_val)

