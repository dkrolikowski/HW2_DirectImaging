import HW2_Fns as fns

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

import glob
from astropy.io import fits

###### ROXs12 #####

# Part 1 -- Read in data for ROXs12

ROXs12 = glob.glob( 'data/ROXs12/*.fits' )

frms_ROXs12, heads_ROXs12 = fns.ReadFiles( ROXs12 )

###

# Part 2/3 -- Find position of star under coronagraph and register the image so that star is in center

centers_ROXs12 = fns.FindStar( frms_ROXs12, [ 470, 611 ] )
regfrms_ROXs12 = fns.Register( frms_ROXs12, centers_ROXs12 )

regmed_ROXs12, regsum_ROXs12 = fns.StackFrames( regfrms_ROXs12 )

###

# Part 4 -- Rotate each image and stack

rotfrms_ROXs12 = fns.Rotate( regfrms_ROXs12, heads_ROXs12 )

rotmed_ROXs12, rotsum_ROXs12 = fns.StackFrames( rotfrms_ROXs12 )

###

# Part 5 -- Calculate and subtract off the radial brightness profile

bfsubfrms_ROXs12 = fns.SubBrightProf( rotfrms_ROXs12 )

bfsubmed_ROXs12, bfsubsum_ROXs12 = fns.StackFrames( bfsubfrms_ROXs12 )

###

# Part 6 -- Subtract off the median PSF from each image

medpsfsubfrms_ROXs12 = fns.SubMedPSF( regfrms_ROXs12, heads_ROXs12, regmed_ROXs12 )

medpsfsubmed_ROXs12, medpsfsubsum_ROXs12 = fns.StackFrames( medpsfsubfrms_ROXs12 )

###### ROXs42 #####

# Part 1 -- Read in data for ROXs42

ROXs42 = glob.glob( 'data/ROXs42B/*.fits' )

frms_ROXs42, heads_ROXs42 = fns.ReadFiles( ROXs42 )

###

# Part 2/3 -- Find position of star under coronagraph and register the image so that star is in center

centers_ROXs42 = fns.FindStar( frms_ROXs42, [ 470, 611 ] )
regfrms_ROXs42 = fns.Register( frms_ROXs42, centers_ROXs42 )

regmed_ROXs42, regsum_ROXs42 = fns.StackFrames( regfrms_ROXs42 )

###

# Part 4 -- Rotate each image and stack

rotfrms_ROXs42 = fns.Rotate( regfrms_ROXs42, heads_ROXs42 )

rotmed_ROXs42, rotsum_ROXs42 = fns.StackFrames( rotfrms_ROXs42 )

###

# Part 5 -- Calculate and subtract off the radial brightness profile

bfsubfrms_ROXs42 = fns.SubBrightProf( rotfrms_ROXs42 )

bfsubmed_ROXs42, bfsubsum_ROXs42 = fns.StackFrames( bfsubfrms_ROXs42 )

###

# Part 6 -- Subtract off the median PSF from each image

medpsfsubfrms_ROXs42 = fns.SubMedPSF( regfrms_ROXs42, heads_ROXs42, regmed_ROXs42 )

medpsfsubmed_ROXs42, medpsfsubsum_ROXs42 = fns.StackFrames( medpsfsubfrms_ROXs42 )

plt.imshow( medpsfsubmed_ROXs42 )
plt.show()


#pd.DataFrame( { 'File Name': np.array(ROXs12), 'X': centers_ROXs12[:,0], 'Y': centers_ROXs12[:,1] } ).to_csv( 'ROXs12_starpos.txt', sep = ' ', index = False )
