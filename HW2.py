import glob, pdb, mpyfit

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits

datapath  = 'data/ROXs12/'
datafiles = glob.glob( datapath + '*' )

test   = fits.open( datafiles[0] )[0].data
frames = np.zeros( ( len(datafiles), test.shape[0], test.shape[1] ) )
heads  = []

for i in range( frames.shape[0] ):
    f = fits.open( datafiles[i] )[0]

    frames[i] = f.data
    heads.append( f.header )

psf = np.median( frames, axis = 0 )

#plt.clf()
#plt.imshow( frames[8] )
#plt.show()

def Gauss2D( X, p ):
    x, y = X

    xterm = np.exp( - ( x - p[0] ) ** 2.0 / ( 2.0 * p[1] ** 2.0 ) )
    yterm = np.exp( - ( y - p[2] ) ** 2.0 / ( 2.0 * p[3] ** 2.0 ) )

    gauss = p[4] * xterm * yterm[:,np.newaxis]

    return gauss

def least( p, args ):
    X, vals, func = args

    dif = vals - func( X, p )

    return dif.ravel()

x = np.arange( 1024 )
y = np.arange( 1024 )

p = [ 300, 8, 500, 8, 1000 ]
plt.clf()
plt.imshow( Gauss2D( ( x, y ), p ) + frames[8] )
plt.show()
