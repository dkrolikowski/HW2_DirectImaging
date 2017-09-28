import glob, pdb, mpyfit

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits

def Gauss2D( X, p ):
    x, y = X

    xterm = np.exp( - ( x - p[0] ) ** 2.0 / ( 2.0 * p[1] ** 2.0 ) )
    yterm = np.exp( - ( y - p[2] ) ** 2.0 / ( 2.0 * p[3] ** 2.0 ) )

    gauss = p[4] * xterm * yterm[:,np.newaxis] + p[5]

    return gauss

def least( p, args ):
    X, vals, func = args

    dif = vals - func( X, p )

    return dif.ravel()

datapath  = '../data/ROXs12/'
datafiles = glob.glob( datapath + '*' )

test   = fits.open( datafiles[0] )[0].data
frames = np.zeros( ( len(datafiles), test.shape[0], test.shape[1] ) )
heads  = []

for i in range( frames.shape[0] ):
    f = fits.open( datafiles[i] )[0]

    frames[i] = f.data
    heads.append( f.header )

# Finding star behind coronagraph
roughx  = 470
roughy  = 612
centers = np.zeros( ( frames.shape[0], 2 ) )

for i in range( frames.shape[0] ):
    tofit = frames[i,roughx-15:roughx+15,roughy-15:roughy+15]
    x = np.arange( tofit.shape[1] )
    y = np.arange( tofit.shape[0] )

    p0 = [ 15, 2, 15, 2, tofit[15,15], np.median(tofit) ]

    fit, res = mpyfit.fit( least, p0, ( ( x, y ), tofit, Gauss2D ) )
    
    centers[i] = np.array( [ fit[2] + roughx - 15, fit[0] + roughy - 15 ] )

print centers

