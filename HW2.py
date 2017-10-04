import glob, pdb, mpyfit, scipy.interpolate, scipy.ndimage

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
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

datapath  = 'data/ROXs12/'
datafiles = glob.glob( datapath + '*' )

## Read in images
test  = fits.open( datafiles[0] )[0].data
frms  = np.zeros( ( len(datafiles), test.shape[0], test.shape[1] ) )
heads = []

for i in range( frms.shape[0] ):
    f = fits.open( datafiles[i] )[0]

    frms[i] = f.data
    heads.append( f.header )

## Part 2 -- Finding star behind coronagraph
starpos = [ 470, 611 ]
centers = np.zeros( ( frms.shape[0], 2 ) )

for i in range( frms.shape[0] ):
    boxw  = 10
    sigma = 2
    tofit = frms[i,starpos[0]-boxw:starpos[0]+boxw,starpos[1]-boxw:starpos[1]+boxw]
    x     = np.arange( tofit.shape[1] )
    y     = np.arange( tofit.shape[0] )

    p0 = [ boxw, sigma, boxw, sigma, tofit[boxw,boxw], np.median(tofit) ]

    fit, res = mpyfit.fit( least, p0, ( ( x, y ), tofit, Gauss2D ) )
    
    centers[i] = np.array( [ fit[2] + starpos[0] - boxw, fit[0] + starpos[1] - boxw ] )

## Part 3 -- Shift images so star is at some x/y position

def register( frame, center, setpos ):

    delx = 1 - ( center[0] - setpos[0] )
    dely = 1 - ( center[1] - setpos[1] )

    xnew = np.linspace( -delx, frame.shape[0] - 1 - delx, frame.shape[0] )
    ynew = np.linspace( -dely, frame.shape[1] - 1 - dely, frame.shape[1] )

    frmint   = scipy.interpolate.interp2d( np.arange( frame.shape[0] ), np.arange( frame.shape[1] ), frame, kind = 'cubic' )
    frmshift = frmint( ynew, xnew )

    return frmshift

def rotate( frame, ang, setpos ):
    ang    = np.radians( ang )
    frmrot = frame.copy()
    frmint = scipy.interpolate.interp2d( np.arange( frame.shape[0] ), np.arange( frame.shape[1] ), frame, kind = 'cubic' )

    xarr  = np.arange( frame.shape[0] )
    yarr  = np.arange( frame.shape[1] )

    for x in xarr:
        for y in yarr:
            xs  = x - setpos[0]
            ys  = y - setpos[1]
            xr  = ( xs * np.cos( ang ) - ys * np.sin( ang ) ) + setpos[0]
            yr  = ( xs * np.sin( ang ) + ys * np.cos( ang ) ) + setpos[1]
            pix = frmint( yr, xr )
            frmrot[x,y] = pix

    return frmrot

regfrms = frms.copy()
rotfrms = frms.copy()

for i in range( frms.shape[0] ):

    PA = heads[i]['PARANG'] + heads[i]['ROTPPOSN'] - heads[i]['EL'] - heads[i]['INSTANGL']
    
    regfrms[i] = register( frms[i], centers[i], [ 511, 511 ] )
    print i, PA
    #rotfrms[i] = rotate( regfrms[i], PA, starpos )
    rotfrms[i] = scipy.ndimage.interpolation.rotate( regfrms[i], PA, reshape = False )

regmed = np.median( regfrms, axis = 0 )
regsum = np.sum( regfrms, axis = 0 )
rotmed = np.median( rotfrms, axis = 0 )
rotsum = np.sum( rotfrms, axis = 0 )

pdb.set_trace()

plt.figure()
plt.imshow( regsum )

plt.figure()
plt.imshow( rotsum )
plt.show()

plt.figure()
plt.imshow( regmed )

plt.figure()
plt.imshow( rotmed )
plt.show()
