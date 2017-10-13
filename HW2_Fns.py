import mpyfit

import numpy as np

from astropy.io import fits
from scipy.interpolate import interp2d, UnivariateSpline
from scipy import ndimage

def ReadFiles( fnames ):

    testim = fits.open( fnames[0] )[0].data
    frms   = np.zeros( ( len( fnames ), testim.shape[0], testim.shape[1] ) )
    heads  = []

    for i in range( frms.shape[0] ):
        f = fits.open( fnames[i] )[0]

        frms[i] = f.data
        heads.append( f.header )

    return frms, heads

def FindStar( frms, starpos ):

    centers = np.zeros( ( frms.shape[0], 2 ) )
    boxw    = 10
    starsig = 2

    for i in range( frms.shape[0] ):
        tofit = frms[i,starpos[0]-boxw:starpos[0]+boxw,starpos[1]-boxw:starpos[1]+boxw]
        x     = np.arange( tofit.shape[1] ) # x/y are flipped for this, but square so doesn't matter
        y     = np.arange( tofit.shape[0] )

        p0    = [ boxw, starsig, boxw, starsig, tofit[boxw,boxw], np.median(tofit) ]

        fit, res   = mpyfit.fit( least, p0, ( ( x, y ), tofit, Gauss2D ) )
        centers[i] = np.array( [ fit[2] + starpos[0] - boxw, fit[0] + starpos[1] - boxw ] )

    return centers

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

def Register( frms, centers ):

    regfrms = frms.copy()
    setpos  = [ frms.shape[1] / 2. - 0.5, frms.shape[2] / 2. - 0.5 ]
    
    for i in range( frms.shape[0] ):
        
        delx = setpos[0] - centers[i,0]
        dely = setpos[1] - centers[i,1]

        xnew = np.linspace( -delx, frms.shape[1] - 1 - delx, frms.shape[1] )
        ynew = np.linspace( -dely, frms.shape[2] - 1 - dely, frms.shape[2] )

        frminterp  = interp2d( np.arange( frms.shape[1] ), np.arange( frms.shape[2] ), frms[i], kind = 'cubic' )
        regfrms[i] = frminterp( ynew, xnew )

    return regfrms

def Rotate( frms, heads ):

    rotfrms = frms.copy()

    for i in range( frms.shape[0] ):

        posang     = heads[i]['PARANG'] + heads[i]['ROTPPOSN'] - heads[i]['EL'] - heads[i]['INSTANGL']
        rotfrms[i] = ndimage.interpolation.rotate( frms[i], -(180+posang), reshape = False )

    return rotfrms

def StackFrames( frms ):

    medfrm = np.median( frms, axis = 0 )
    sumfrm = np.sum( frms, axis = 0 )

    return medfrm, sumfrm

def SubBrightProf( frms ):

    subfrms = frms.copy()
    
    oneside = np.arange( 1024 ) - 511.5
    distmat = np.sqrt( oneside ** 2.0 + np.array( [ oneside ] ).T ** 2.0 )
    dist1d  = distmat.ravel()
    dist1ds = dist1d[np.argsort(dist1d)]
    distfit = np.linspace( 1, 723, 724/2 )

    for i in range( frms.shape[0] ):
        frm1d   = frms[i].ravel()
        frm1ds  = frm1d[np.argsort(dist1d)]

        medians = []
        for d in distfit:
            region = np.where( ( dist1ds >= d - 1.0 ) & ( dist1ds <= d + 1.0 ) )[0]
            medians.append( np.median( frm1ds[region] ) )

        splfit      = UnivariateSpline( distfit, medians )
        subfrms[i] -= splfit( distmat )

    return subfrms
        
def SubMedPSF( frms, heads, PSF ):

    subfrms = frms.copy()

    for i in range( frms.shape[0] ):

        subfrms[i] -= PSF

    subfrms = Rotate( subfrms, heads )

    return subfrms

        

