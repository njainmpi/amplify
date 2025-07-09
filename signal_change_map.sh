#!/bin/sh


SIGNAL_CHANGE_MAPS () {

        input_4d_data=$1
        local file_for_parameter_calculation=$2
        local baseline_duration_in_min=$3
        local duration=$4


        total_reps=$(awk -F'=' '/##\$PVM_NRepetitions=/{print $2}' $file_for_parameter_calculation/method)
        TotalScanTime=$(awk -F'=' '/##\$PVM_ScanTime=/{print $2}' $file_for_parameter_calculation/method)
        slice_count=$(awk -F'=' '/##\$NSLICES=/{print $2}' $file_for_parameter_calculation/acqp)
 
        #here the awk will look at the number of slices acquired using the information located in the methods file   
            
        # 07.08.2024 Estimating Volume TR
        VolTR_msec=$(echo "scale=0; $TotalScanTime/$total_reps" | bc)
        VolTR=$(echo "scale=0; $VolTR_msec/1000" | bc)


        No_of_Vols_in_pre_PACAP_injection=$((baseline_duration_in_min * 60 / VolTR))
        echo $No_of_Vols_in_pre_PACAP_injection

        base_start=$(( 0 + 80 ))
        base_end=$(( No_of_Vols_in_pre_PACAP_injection - 80 ))

        echo $base_start $base_end
        echo $VolTR


        # Step 1: Compute the baseline image (mean of the first 600 repetitions)
        3dTstat -mean -prefix baseline_image.nii.gz ${input_4d_data}"[${base_start}..${base_end}]"

        # Initialize an empty list to store the intermediate processed images
        processed_images=()
        volumes_per_block=$((VolTR * duration))
        total_reps_idx=$((total_reps - 1))


        # Step 2: Process blocks of 60 repetitions, starting from time point 601
        for start_idx in $(seq 0 $volumes_per_block $total_reps_idx); do
                end_idx=$((start_idx + volumes_per_block - 1))
                
                # Create a meaningful label for this block of repetitions
                label="${start_idx}_to_${end_idx}"
                
                # Step 2a: Compute the mean for the current block of 60 repetitions
                3dTstat -mean -prefix block_mean_${label}.nii.gz ${input_4d_data}"[${start_idx}..${end_idx}]"
                
                # Step 2b: Subtract the baseline and divide by the baseline
                3dcalc -a block_mean_${label}.nii.gz -b baseline_image.nii.gz -expr '(a-b)/b' -prefix ratio_processed_${label}.nii.gz

                #Step 2c: Converting into percent by multiplying it by 100
                3dcalc -a ratio_processed_${label}.nii.gz -expr 'a*100' -prefix processed_${label}.nii.gz
                    
                # Add the processed image to the list
                processed_images+=("processed_${label}.nii.gz")
                
                echo "Processed block: ${start_idx} to ${end_idx}"
        done




        # Step 3: Combine all the processed blocks into a single 4D image
        3dTcat -prefix Signal_Change_Map_premask.nii.gz "${processed_images[@]}"

        fslmaths Signal_Change_Map_premask.nii.gz -mas mask_mean_mc_func.nii.gz Signal_Change_Map.nii.gz

        echo "All blocks processed and combined into final 4D image: Signal_Change_Map.nii.gz"


        rm -f *block* processed* ratio_processed*

        # Step 4: Create the output directory for screenshots
        mkdir -p "$output_dir"


        # Step 5: Extract screenshots for each volume in the final processed 4D image
        num_volumes=$(3dinfo -nv Signal_Change_Map.nii.gz)  # Get the number of volumes

        for ((slice_idx=0; slice_idx<slice_count; slice_idx++)); do     

                echo "Making Video for Slice Number $slice_idx"
                python ~/Desktop/Github/amplify/Making_videos_for_SCM.py $7 Signal_Change_Map.nii.gz $slice_idx
        done



}

