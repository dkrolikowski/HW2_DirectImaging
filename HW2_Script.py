import HW2_Fns as fns

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

import glob, pickle
from astropy.io import fits
import pandas as pd

nirc2pixscale = 0.009942

###### ROXs12 and ROXs42 #####

# Part 1 -- Read in data

ROXs12 = glob.glob( 'data/ROXs12/*.fits' )
frms_ROXs12, heads_ROXs12 = fns.ReadFiles( ROXs12 )

ROXs42 = glob.glob( 'data/ROXs42B/*.fits' )
frms_ROXs42, heads_ROXs42 = fns.ReadFiles( ROXs42 )

###

# Part 2/3 -- Find position of star under coronagraph and register the image so that star is in center
# Also write out star positions to a text file

centers_ROXs12 = fns.FindStar( frms_ROXs12, [ 470, 611 ] )
regfrms_ROXs12 = fns.Register( frms_ROXs12, centers_ROXs12 )
regmed_ROXs12, regsum_ROXs12 = fns.StackFrames( regfrms_ROXs12 )

pd.DataFrame( { 'File': ROXs12, 'X': centers_ROXs12[:,0], 'Y': centers_ROXs12[:,1] }).to_csv( 'ROXs12starpos.txt', index = False )

centers_roxs42 = fns.FindStar( frms_ROXs42, [ 470, 611 ] )
regfrms_ROXs42 = fns.Register( frms_ROXs42, centers_ROXs42 )
regmed_ROXs42, regsum_ROXs42 = fns.StackFrames( regfrms_ROXs42 )

pd.DataFrame( { 'File': ROXs42, 'X': centers_ROXs42[:,0], 'Y': centers_ROXs42[:,1] }).to_csv( 'ROXs42starpos.txt', index = False )

###

# Part 4 -- Rotate each image and stack

rotfrms_ROXs12 = fns.Rotate( regfrms_ROXs12, heads_ROXs12 )
rotmed_ROXs12, rotsum_ROXs12 = fns.StackFrames( rotfrms_ROXs12 )

rotfrms_ROXs42 = fns.Rotate( regfrms_ROXs42, heads_ROXs42 )
rotmed_ROXs42, rotsum_ROXs42 = fns.StackFrames( rotfrms_ROXs42 )

###

# Part 5 -- Calculate and subtract off the radial brightness profile

bfsubfrms_ROXs12 = fns.SubBrightProf( rotfrms_ROXs12 )
bfsubmed_ROXs12, bfsubsum_ROXs12 = fns.StackFrames( bfsubfrms_ROXs12 )

bfsubfrms_ROXs42 = fns.SubBrightProf( rotfrms_ROXs42 )
bfsubmed_ROXs42, bfsubsum_ROXs42 = fns.StackFrames( bfsubfrms_ROXs42 )

###

# Part 6 -- Subtract off the median PSF from each image

medpsfsubfrms_ROXs12 = fns.SubMedPSF( regfrms_ROXs12, heads_ROXs12, regmed_ROXs12 )
medpsfsubmed_ROXs12, medpsfsubsum_ROXs12 = fns.StackFrames( medpsfsubfrms_ROXs12 )

medpsfsubfrms_ROXs42 = fns.SubMedPSF( regfrms_ROXs42, heads_ROXs42, regmed_ROXs42 )
medpsfsubmed_ROXs42, medpsfsubsum_ROXs42 = fns.StackFrames( medpsfsubfrms_ROXs42 )

###

# Part 7 -- Subtract off best possible PSF from other star

comppsfsubfrms_ROXs12, bestcomp_ROXs12 = fns.BestCompPSF( regfrms_ROXs12, heads_ROXs12, regfrms_ROXs42 )
comppsfsubmed_ROXs12, comppsfsubsum_ROXs12 = fns.StackFrames( comppsfsubfrms_ROXs12 )

pd.DataFrame( { 'ROXs12File': ROXs12, 'ROXs42 Best Comp': ROXs42[bestcomp_ROXs12] } ).to_csv( 'ROXs12BestPSF.txt' )

comppsfsubfrms_ROXs42, bestcomp_ROXs42 = fns.BestCompPSF( regfrms_ROXs42, heads_ROXs42, regfrms_ROXs12 )
comppsfsubmed_ROXs42, comppsfsubsum_ROXs42 = fns.StackFrames( comppsfsubfrms_ROXs42 )

pd.DataFrame( { 'ROXs42File': ROXs42, 'ROXs12 Best Comp': ROXs12[bestcomp_ROXs42] } ).to_csv( 'ROXs42BestPSF.txt' )

###

# Part 8 -- Measure companion sep/PA

# ROXs12

ROXs12comp = [334, 538]

rotsep, rotPA         = fns.PlanetSepPA( rotmed_ROXs12, ROXs12comp, nirc2pixscale )
bpsep, bpPA           = fns.PlanetSepPA( bfsubmed_ROXs12, ROXs12comp, nirc2pixscale )
medpsfsep, medpsfPA   = fns.PlanetSepPA( medpsfsubmed_ROXs12, ROXs12comp, nirc2pixscale )
comppsfsep, comppsfPA = fns.PlanetSepPA( comppsfsubmed_ROXs12, ROXs12comp, nirc2pixscale )

# ROXs42

ROXs42comp = [512, 393]

rotsep, rotPA         = fns.PlanetSepPA( rotmed_ROXs42, ROXs42comp, nirc2pixscale )
bpsep, bpPA           = fns.PlanetSepPA( bfsubmed_ROXs42, ROXs42comp, nirc2pixscale )
medpsfsep, medpsfPA   = fns.PlanetSepPA( medpsfsubmed_ROXs42, ROXs42comp, nirc2pixscale )
comppsfsep, comppsfPA = fns.PlanetSepPA( comppsfsubmed_ROXs42, ROXs42comp, nirc2pixscale )
