"""
Holds functions useful for finding the nearest MC to a location in a galaxy.
"""

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
from getMapValue import getValue

def int2mass(x, aco, res=150.0):
    area = (res/2.0)**2*np.pi/np.log(2.0)
    y = x * aco * area
    return(y)
    
def mass2int(mass, aco, res=150.0):
    area = (res/2.0)**2*np.pi/np.log(2.0)
    inten = mass / (aco * area)
    return(inten)

def angDistToPc(x,galDist):
    return(galDist*10**6*np.tan(x*np.pi/180))


def findNearest(x1,x2,y1,y2):
    # Where x1 & y1 are 1d arrays of map coordinates
    # x2,y2 are the coords of the SNe
    n = len(x1)
    m = len(x2)
    if n != 0:
        x1Vec = np.tile(x1, (m,1)) #constant x along column
        y1Vec = np.tile(y1, (m,1)) #constant y along column
        x2Vec = np.tile(x2, (n,1)) #constant x along column
        y2Vec = np.tile(y2, (n,1)) #constant y along column
        x2Vec = np.transpose(x2Vec) #constant x along rows
        y2Vec = np.transpose(y2Vec) #constant y along rows

        dist = np.sqrt((x1Vec-x2Vec)**2 + (y1Vec-y2Vec)**2)    

        mindist = np.nanmin(dist, axis = 1)
    else: 
        mindist = np.full_like(len(x2),float("nan"))
    
    return(mindist)
'''
def nearestMCMethod(galaxy, co_file, err_file, mass_file, centerCoord, pa, incl, galDist, otherra, otherdec, othername=None):
    
    if os.path.isfile(co_file):

        # read in fits files
        area = (150.0/2.0)**2*np.pi/np.log(2.0)
        hdu_int = pyfits.open(co_file)
        intMap  = hdu_int[0].data

        hdu_err = pyfits.open(err_file)
        errMap  = hdu_err[0].data

        hdu_mass = pyfits.open(mass_file)
        massMap  = hdu_mass[0].data 

        #Convert x & y pixels to ra and dec
        wcs      = WCS(hdu_int[0].header, naxis=2)
        naxis    = wcs._naxis # size of image naxis[0] = x and [1] = y
        grid     = np.indices((naxis[1],naxis[0]))
        ra, dec  = wcs.wcs_pix2world(grid[1],grid[0],0)

        #deproject ra and dec to dx and dy
        radius, projang, dx, dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=ra, dec=dec, return_offset=True)

        #flatten data structures 
        f_int  = intMap.flatten()
        f_err  = errMap.flatten()
        f_mass = massMap.flatten()
        f_ra   = ra.flatten()
        f_dec  = dec.flatten()    
        f_dx   = dx.flatten()
        f_dy   = dy.flatten()

        #remove nans
        keep  = np.where(np.isfinite(f_int))
        ra    = f_ra[keep]
        dec   = f_dec[keep]
        inten = f_int[keep]
        err   = f_err[keep]
        mass  = f_mass[keep]
        dx    = f_dx[keep]
        dy    = f_dy[keep]

        SNR = []
        for i in range(len(inten)):
            if err[i] == 0.0:
                SNR.append(0.0)
            elif inten[i] < 0.0:
                SNR.append(0.0)           
            else:
                SNR.append(inten[i]/err[i])
        SNR = np.array(SNR)      
        
        print("at Mass cutoff A for", galaxy)
        
        idx  = (mass > 10**5.5) * (SNR > 3.0)
        mass = mass[idx]
        dx   = dx[idx]
        dy   = dy[idx]

        print("Pixels with Mass > 10**5.5", len(mass))

        otherRad, otherPA, otherdx, otherdy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=otherra, dec=otherdec, return_offset=True) 
        
        nearestMCx55 = findNearest(dx, otherdx, dy, otherdy)
        nearestMC55 = angDistToPc(nearestMCx55,galDist)
        
        print("Nearest 55", nearestMC55, galaxy)        
        
        idx  = (mass > 10**6.5) 
        mass = mass[idx]
        dx   = dx[idx]
        dy   = dy[idx]

        print("Pixels with Mass > 10**6.5", len(mass))            

        nearestMCx65 = findNearest(dx, otherdx, dy, otherdy)
        nearestMC65 = angDistToPc(nearestMCx65,galDist)
                   

        print("Nearest 65", nearestMC65, galaxy)        
        
        print("done with", galaxy, " ")
	return(nearestMC55, nearestMC65)

    else:
        print("No file for ", galaxy)

        n55 = float("nan")
        n66 = float("nan")
	return(n55, n65)'''


def nearestMCMethodSigMol(galaxy, co_file, err_file, aco_file, centerCoord, pa, incl, galDist, otherra, otherdec, othername=None):
    
    if os.path.isfile(co_file):

        # read in fits files
        hdu_int = pyfits.open(co_file)
        intMap  = hdu_int[0].data

        hdu_err = pyfits.open(err_file)
        errMap  = hdu_err[0].data
        
        if os.path.isfile(aco_file):
            hdu_aco = pyfits.open(aco_file)
            acoMap  = hdu_aco[0].data
        else:
            acoMap = intMap * 6.7

        # Accounting for inclination affects in NGC4945 aco map
        if galaxy == "NGC4945":
            acoMap = intMap * 6.7

        #Convert x & y pixels to ra and dec
        wcs      = WCS(hdu_int[0].header, naxis=2)
        naxis    = wcs._naxis # size of image naxis[0] = x and [1] = y
        grid     = np.indices((naxis[1],naxis[0]))
        ra, dec  = wcs.wcs_pix2world(grid[1],grid[0],0)
        
        if incl == 90:
            incl = 80
         
        #deproject ra and dec to dx and dy
        radius, projang, dx, dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=ra, dec=dec, return_offset=True)

        #flatten data structures 
        f_int  = intMap.flatten()
        f_err  = errMap.flatten()
        f_aco  = acoMap.flatten()
        f_ra   = ra.flatten()
        f_dec  = dec.flatten()    
        f_dx   = dx.flatten()
        f_dy   = dy.flatten()

        #remove nans
        keep   = np.where(np.isfinite(f_int))
        ra     = f_ra[keep]
        dec    = f_dec[keep]
        inten  = f_int[keep]
        aco    = f_aco[keep]
        sigmol = aco * inten * np.cos(incl * np.pi/180.0)
        err    = f_err[keep]
        dx     = f_dx[keep]
        dy     = f_dy[keep]


        threeSigma = aco * 3.0 * err * np.cos(incl * np.pi/180.0)
        threesigCutOff = np.nanmedian(threeSigma)
        print("Three Sigma Cutoff: ", threesigCutOff)

        SNR = []
        for i in range(len(inten)):
            if err[i] == 0.0:
                SNR.append(0.0)
            elif inten[i] < 0.0:
                SNR.append(0.0)           
            else:
                SNR.append(inten[i]/err[i])
        SNR = np.array(SNR)      
        
        #print("at SigMol cutoff A for", galaxy)
        
        idx  = (sigmol > threesigCutOff) * (SNR > 3.0)
        #sm = sigmol[idx]
        dx   = dx[idx]
        dy   = dy[idx]

        #print("Pixels with SigMol > threeSigmaValue", len(sm))

        otherRad, otherPA, otherdx, otherdy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=otherra, dec=otherdec, return_offset=True) 

        sm = []
        for l in  range(len(otherra)):
            sm.append(getValue(co_file, otherra[l], otherdec[l]))

        nearestSigMol = findNearest(dx, otherdx, dy, otherdy)
        nearestSM = angDistToPc(nearestSigMol,galDist)
        
        #print("Nearest sigmol", nearestSM, galaxy)        
               
        print("done with", galaxy, " ")  
        return(nearestSM,sm, threesigCutOff)

    else:
        print("No file for ", galaxy)

        n55 = float("nan")
        sm =float("nan")
        threesigCutOff = float("nan")
        return(n55,sm, threesigCutOff)


