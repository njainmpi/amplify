#!/bin/bash

BIAS_CORRECTED_IMAGE () {

        input_file=$1
        local b_val=$2

        # Step 1: Creating a mask from 3D mean image
        fslmaths ${input_file} -thrp 30 -bin initial_${input_file}

        # Here we are cleaning the automatic mask to suit it more to the mean EPI image
        echo "Modify the mask"
        echo ""
        echo ""

        echo "Please save the mask by the name 'mask_${input_file}'. "
        fsleyes initial_anatomy.nii.gz ${input_file}


        # Step 2: Coil (B1) inhomogeneity correction of EPI using N4 method 
        # Here we will be using the mask that we created on mean functional image.

        N4BiasFieldCorrection -d 3 -i mc_func.nii.gz -o N4_${input_file} -c [100x100x100,0.0001] -b [$b_val,3] -s 2 -x mask_${input_file}
        # -c [100x100x100,0.0001] means optimizing at 3 scales, each scale has 100 iterations
        # -b [54,3] means start with 32 points scale (equiv 20mm coil divided by 0.375mm resolution) with 3rd order b-spline
        # -s 2 means scale jump by factor of 2 in each iteration


        # Step 3: Applying final mask on 4D motion corrected data to clean the raw time series
        fslmaths N4_${input_file} -mas mask_${input_file} cleaned_${input_file}


}