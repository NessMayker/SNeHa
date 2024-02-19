import numpy as np

def gc(gal_ra, gal_dec, incl, pa, sn_ra, sn_dec):
    """
    Calculate deprojected radii and projected angles in a disk.    

    Parameters
    ----------
    gal_ra  : float
        ra coord of the galactic center in degrees
    gal_dec : float
        dec coord of the galactic center in degrees
    incl    : float
        galaxy inclination angle in degrees
    pa      : float
        galaxy position angle in degrees
    sn_ra   : float
        ra coord of the supernova in degrees
    sn_dec  : float
        dec coord of the supernova in degrees
    
    Returns
    -------
    radius_deg : float
        galactocentric radius in degrees
    radius_arcsec : float
        galactocentric radius in arcseconds
    proj_ang : float
        projection angle in degrees (theta in (r,theta))
    """

    # recast the ra and dec arrays in term of the center coordinates
    # arrays are now in degrees from the center
    # offsets in ra and dec on celestial sphere (might misbehave at poles)
    dx_deg = (sn_ra - gal_ra) * np.cos(np.deg2rad(gal_dec))
    dy_deg = sn_dec - gal_dec

    # rotation angle (rotate x-axis up to the major axis)
    rotangle = np.pi/2.0 - np.deg2rad(pa)

    # create deprojected coordinate grids
    # offsets after deprojection with the coordinates we are interested in
    deprojdx_deg = (dx_deg * np.cos(rotangle) +
                    dy_deg * np.sin(rotangle))
    deprojdy_deg = (dy_deg * np.cos(rotangle) -
                    dx_deg * np.sin(rotangle))
    deprojdy_deg /= np.cos(np.deg2rad(incl))

    # make map of deprojected distance from the center
    radius_deg = np.sqrt(deprojdx_deg**2 + deprojdy_deg**2)
    radius_arcsec = radius_deg * 3600
    
    # make map of angle w.r.t. position angle
    projang_deg = np.rad2deg(np.arctan2(deprojdy_deg, deprojdx_deg))

    return radius_deg, radius_arcsec, projang_deg


def gcr(radius_deg, dist):
    """
    Calculates the galactocentric radius of a galaxy
    
    Parameters
    ----------
    radius_deg : float
        radius of the galaxy in degrees
        
    dist       : float
        distance to galaxy in kpc
    """   
            
    value = dist * 1000 *  np.tan(radius_deg)
    return (value)
