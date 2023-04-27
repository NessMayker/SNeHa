
"""For generating values for visualization set return distVals = False
For model 1: set modelType = 1"""

import os
import numpy as np
import astropy.io.fits as pyfits
from astropy.io import ascii
from astropy.table import Table
from astropy.wcs import WCS
from reproject import reproject_interp

import sys
sys.path.append('./Python')
from deprojectGalaxy import deproject
from normalize import norm
from FindNearestMC import angDistToPc, findNearest
from ReturnMapData import returnMapData

def runModels(galaxy, image, centerCoord, pa, incl, galDist, modelType = 1, starLight = None, starRa = None, starDec = None, expSize = 100):
    r_ra, r_dec, r_dx, r_dy, r_sm = [],[],[],[],[]
   
    if os.path.isfile(image):
        inten, ra, dec, dx, dy = returnMapData(image, centerCoord=centerCoord, incl=incl, pa=pa)
      
        #if model is random
        if modelType == 1:
            numrows = len(inten)    

            rand = np.random.randint(low=0, size = expSize, high=numrows)

            r_ra  = ra[rand]
            r_dec = dec[rand]
            r_dx  = dx[rand]
            r_dy  = dy[rand]
            r_int = inten[rand] * np.cos(incl*np.pi/180.)
            

        #if model is gas density weighted
        elif modelType == 2:
            intArr = np.clip(inten, 0.0, None)
            total = sum(intArr)
            prob  = intArr/total 
            prob  = norm(prob)
            n = len(intArr)

            indicies = np.arange(n, dtype=int)
            rand = np.random.choice(indicies, size = expSize, p=prob)

            r_ra  = ra[rand]
            r_dec = dec[rand]
            r_dx  = dx[rand]
            r_dy  = dy[rand]
            #r_int = inten[rand]
            r_int = inten[rand] * np.cos(incl*np.pi/180.)

        #if model is starlight density weighted
        elif modelType == 3:
            intArr = np.clip(starLight, 0.0, None)
            total = sum(intArr)
            prob  = intArr/total 
            prob  = norm(prob)
            n = len(intArr)
            indicies = np.arange(n, dtype=int)
            rand = np.random.choice(indicies, size = expSize, p=prob)

            starRad, starPA, starDx, starDy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=starRa, dec=starDec, return_offset=True)

            r_ra  = starRa[rand]
            r_dec = starDec[rand]
            r_dx  = starDx[rand]
            r_dy  = starDy[rand]
            #r_int = inten[rand]
            r_int = inten[rand] * np.cos(incl*np.pi/180.)

        else: print("Wrong model choice, should be 1, 2, or 3.")
        print("r_int:",r_int)
        return(r_ra, r_dec, r_dx, r_dy, r_int)


    else:
        print("Something wrong with file or missing")


