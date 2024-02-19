
import numpy as np

def angDistToPc(x,galDist):
    '''Converts an angular distance to parsecs when given galaxy distance in Mpc'''
    pc = galDist * 10**6 * np.tan( x * np.pi / 180 )
    return(pc)