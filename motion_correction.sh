#!/bin/bash

MOTION_CORRECTION () {
      
    3dvolreg -prefix $3 -base $1 -verbose -1Dfile motion.1D -1Dmatrix_save mats -linear -twopass -float -maxdisp1D rmsabs.1D $2
    3dAFNItoNIFTI $3'+orig.BRIK'
    if [ -f mc_func.nii ]; then
        gzip -1 mc_func.nii
    else
        echo "gzip exists"
    fi
    # 1Dplot -volreg -sep motion.1D
    1dplot -xlabel Time -ylabel "Translation (mm)" -title "3dvolreg translations" -dx 1 -jpgs 640x144 rest_translation 'motion.1D[3..5]'
    1dplot -xlabel Time -ylabel "Rotations (degree)" -title "3dvolreg rotations" -dx 1 -jpgs 640x144 rest_rotation 'motion.1D[0..2]'
}