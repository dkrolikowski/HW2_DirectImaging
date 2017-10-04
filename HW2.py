import glob, pdb, mpyfit, scipy.interpolate, scipy.ndimage, pickle

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

regfrms = frms.copy()
rotfrms = frms.copy()

for i in range( frms.shape[0] ):

    print i
    PA = heads[i]['PARANG'] + heads[i]['ROTPPOSN'] - heads[i]['EL'] - heads[i]['INSTANGL']
    
    regfrms[i] = register( frms[i], centers[i], [ 511, 511 ] )
    rotfrms[i] = scipy.ndimage.interpolation.rotate( regfrms[i], -(180 + PA), reshape = False )

pickle.dump( regfrms, open( 'ims_reg.pkl', 'wb' ) )
pickle.dump( rotfrms, open( 'ims_rot.pkl', 'wb' ) )

## Tests for radial profile

oneside = np.arange( 1024 ) - 511.5
distmat = np.sqrt( oneside ** 2.0 + np.array( [ oneside ] ).T ** 2.0 )

i = 8
dist1d = distmat.ravel()
frm1d  = rotfrms[8].ravel()

dist1ds = dist1d[np.argsort(dist1d)]
frm1ds  = frm1d[np.argsort(dist1d)]

distsplit = np.unique( dist1ds )
frmsplit  = np.split( frm1ds, np.where( np.diff( dist1ds ) > 0 )[0] + 1 )

meds = np.zeros( distsplit.size )
for i in range( meds.size ):
    meds[i] = np.median( frmsplit[i] )
