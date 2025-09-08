#!/bin/bash

BIAS_CORRECTED_IMAGE () {

        input_file=$1
        local b_val=$2

        # Step 1: Creating a mask from 3D mean image
        # fslmaths ${input_file} -thrp 30 -bin initial_${input_file}

        # Here we are cleaning the automatic mask to suit it more to the mean EPI image
        echo "Modify the mask"
        echo ""
        echo ""


        3dTstat -mean -prefix mean_vol500_mc_func_mask.nii.gz ${input_file}'[0..499]'

        echo "Please save the mask by the name 'mask_${input_file}'. "
        fsleyes mean_vol500_mc_func_mask.nii.gz
        
        # if [ -f mask_${input_file} ]; then 
        #         echo -e "\033[32mCleaned Mask already exists.\033[0m "
        #         echo -e "\033[32mConfirm the correctness of mask before proceeding to data analysis.\033[0m "
        #         fsleyes ${input_file} mask_${input_file}
        # else 
        #         echo -e "\033[31mCleaned Mask doesn't exist. Please make a cleaned mask first.\033[0m"
        #         cp initial_mean_mc_func.nii.gz mask_${input_file}
        #         fsleyes ${input_file} mask_${input_file}
        # fi 
        # Step 2: Coil (B1) inhomogeneity correction of EPI using N4 method 
        # Here we will be using the mask that we created on mean functional image.

        N4BiasFieldCorrection -d 3 -i ${input_file} -o N4_${input_file} -c [100x100x100,0.0001] -b [$b_val,3] -s 2 -x mask_${input_file}
        # -c [100x100x100,0.0001] means optimizing at 3 scales, each scale has 100 iterations
        # -b [54,3] means start with 32 points scale (equiv 20mm coil divided by 0.375mm resolution) with 3rd order b-spline
        # -s 2 means scale jump by factor of 2 in each iteration


        # Step 3: Applying final mask on 4D motion corrected data to clean the raw time series
        fslmaths N4_${input_file} -mas mask_${input_file} cleaned_N4_${input_file}
        fslmaths $3 -mas mask_${input_file} cleaned_${3}

}