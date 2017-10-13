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

## Part 3/4 -- Shift images so star is at some x/y position

def Register( frame, center, setpos ):

    delx     = 1 - ( center[0] - setpos[0] )
    dely     = 1 - ( center[1] - setpos[1] )

    xnew     = np.linspace( -delx, frame.shape[0] - 1 - delx, frame.shape[0] )
    ynew     = np.linspace( -dely, frame.shape[1] - 1 - dely, frame.shape[1] )

    frmint   = scipy.interpolate.interp2d( np.arange( frame.shape[0] ), np.arange( frame.shape[1] ), frame, kind = 'cubic' )
    frmshift = frmint( ynew, xnew )

    return frmshift

# regfrms = frms.copy()
# rotfrms = frms.copy()

# for i in range( frms.shape[0] ):

#     print i
#     PA = heads[i]['PARANG'] + heads[i]['ROTPPOSN'] - heads[i]['EL'] - heads[i]['INSTANGL']

#     regfrms[i] = Register( frms[i], centers[i], [ 511, 511 ], -(180+PA) )
#     rotfrms[i] = scipy.ndimage.interpolation.rotate( regfrms[i], -(180+PA), reshape = False )

# pickle.dump( regfrms, open( 'ims_reg.pkl', 'wb' ) )
# pickle.dump( rotfrms, open( 'ims_rot.pkl', 'wb' ) )

regfrms = pickle.load( open( 'ims_reg.pkl', 'rb' ) )
rotfrms = pickle.load( open( 'ims_rot.pkl', 'rb' ) )

regmed  = np.median( regfrms, axis = 0 )
regsum  = np.sum( regfrms, axis = 0 )
rotmed  = np.median( rotfrms, axis = 0 )
rotsum  = np.sum( rotfrms, axis = 0 )

## Part 5 -- Subtract off the brightness profile

def BrightProfile( frm ):

    oneside = np.arange( 1024 ) - 511.5
    distmat = np.sqrt( oneside ** 2.0 + np.array( [ oneside ] ).T ** 2.0 )

    dist1d  = distmat.ravel()
    frm1d   = frm.ravel()
    dist1ds  = dist1d[np.argsort(dist1d)]
    frm1ds   = frm1d[np.argsort(dist1d)]

    distfit = np.linspace( 1, 723, 724 / 2 )
    medians = []

    for d in distfit:
        region = np.where( ( dist1ds >= d - 1.0 ) & ( dist1ds <= d + 1.0 ) )[0]
        medians.append( np.median( frm1ds[region] ) )

    splfit  = scipy.interpolate.UnivariateSpline( distfit, medians )

    return frm - splfit( distmat )

# subfrms = rotfrms.copy()

# for i in range( subfrms.shape[0] ):
    
#     subfrms[i] = BrightProfile( subfrms[i] )

# pickle.dump( subfrms, open( 'ims_sub.pkl', 'wb' ) )

subfrms = pickle.load( open( 'ims_sub.pkl', 'rb' ) )
    
submed = np.median( subfrms, axis = 0 )
subsum = np.sum( subfrms, axis = 0 )

## Part 6 -- Median PSF subtraction

psfsubfrms = regfrms.copy()

for i in range( regfrms.shape[0] ):

    PA = heads[i]['PARANG'] + heads[i]['ROTPPOSN'] - heads[i]['EL'] - heads[i]['INSTANGL']

    psfsubfrms[i] -= regmed
    psfsubfrms[i] = scipy.ndimage.interpolation.rotate( psfsubfrms[i], -(180+PA), reshape = False )

psfsubmed = np.median( psfsubfrms, axis = 0 )
psfsubsum = np.sum( psfsubfrms, axis = 0 )

plt.imshow( psfsubmed )
plt.show()
