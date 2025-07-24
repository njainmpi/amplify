#!/bin/sh


Signal_Change_Map (){

        # input_4d_data=sm_func.nii.gz
        # file_for_parameter_calculation="/Users/njain/Desktop"
        # baseline_duration_in_min=10
        # duration=5
        # injection_duration=10

        input_4d_data=$1
        file_for_parameter_calculation=$2
        baseline_duration_in_min=$3
        duration=$4
        injection_duration=$5


        ##----------------------------------------------------------------------
        ##Defining Colours for Display Output
        # Terminal colors
        RED='\033[0;31m'
        GREEN='\033[0;32m'
        BLUE='\033[0;34m'
        YELLOW='\033[1;33m'
        GRAY='\033[1;37m'
        BOLD='\033[1m'
        NC='\033[0m'
        ##----------------------------------------------------------------------
        # if [ $# -lt 5 ]; then
              
        #         echo -e "  bash $0 <input_4d_data> <file_for_parameter_calculation> <baseline_duration_in_min> <duration> <injection_duration>\n"
               
        #         echo -e "  input_4d_data                → Path to 4D fMRI NIfTI file"
        #         echo -e "  file_for_parameter_calculation → Directory or file used to extract parameters"
        #         echo -e "  baseline_duration_in_min    → Baseline duration in minutes"
        #         echo -e "  duration                    → Total duration in minutes"
        #         echo -e "  injection_duration          → Injection period duration in minutes"
        #         echo ""
                
        # fi


        source ./scm_visual.sh

        total_reps=$(awk -F'=' '/##\$PVM_NRepetitions=/{print $2}' $file_for_parameter_calculation/method)
        TotalScanTime=$(awk -F'=' '/##\$PVM_ScanTime=/{print $2}' $file_for_parameter_calculation/method)
        slice_count=$(awk -F'=' '/##\$NSLICES=/{print $2}' $file_for_parameter_calculation/acqp)
        Sequence="$(grep -A1 'ACQ_protocol_name=( 64 )' "$file_for_parameter_calculation/acqp" \
                | grep -v -e 'ACQ_protocol_name=( 64 )' -e '--' \
                | sed -e 's/[<>]//g' -e 's/_.*$//')"
        TotalScanTime_in_min=$(echo "scale=2; $TotalScanTime/60000" | bc)
       # Convert to integer (drop decimal part)
        TotalScanTime_in_min_integer=$(echo "$TotalScanTime_in_min / 1" | bc)

        #here the awk will look at the number of slices acquired using the information located in the methods file   
           
        # 07.08.2024 Estimating Volume TR
        VolTR_msec=$(echo "scale=0; $TotalScanTime/$total_reps" | bc)
        VolTR=$(echo "scale=0; $VolTR_msec/1000" | bc)

        echo ""
        echo ""
        No_of_Vols_in_pre_PACAP_injection=$((baseline_duration_in_min * 60 / VolTR))
        echo -e "${RED}No of Volumes${NC} in data are ${GREEN}$total_reps${NC} with ${RED}Volume TR${NC} of ${GREEN}$VolTR sec${NC} with ${RED}No of Volumes in Pre_injection${NC} to be ${GREEN}$No_of_Vols_in_pre_PACAP_injection${NC}."
        Seconds_to_be_discarded=80
        No_of_Vols_in_pre_PACAP_injection_to_be_discarded=$((Seconds_to_be_discarded / VolTR))
        echo -e "${RED}No of Volumes to be discard is $No_of_Vols_in_pre_PACAP_injection_to_be_discarded for $Seconds_to_be_discarded Seconds.${NC}"


        base_start=$(( 0 + No_of_Vols_in_pre_PACAP_injection_to_be_discarded ))
        base_end=$(( No_of_Vols_in_pre_PACAP_injection - No_of_Vols_in_pre_PACAP_injection_to_be_discarded ))

        echo ""
        echo -e "${BLUE}For current data, baseline calculation starts at Volume No${NC} ${GREEN}$base_start${NC} ${BLUE}and finishes at${NC} ${GREEN}$base_end${NC}."


        draw_case "Sequence used $Sequence" $TotalScanTime_in_min_integer 1 $baseline_duration_in_min $injection_duration $VolTR $seq

        # Step 1: Compute the baseline image (mean of the first 600 repetitions)
        3dTstat -mean -prefix baseline_image.nii.gz ${input_4d_data}"[${base_start}..${base_end}]"

        # Initialize an empty list to store the intermediate processed images
        processed_images=()
        duration_seconds=$((duration * 60))
        volumes_per_block=$(echo "$duration_seconds / $VolTR" | bc)
        total_reps_idx=$((total_reps))

        # === VALIDATE ===
        if [[ -z "$VolTR" || -z "$total_reps_idx" || -z "$duration" ]]; then
                echo "Usage: $0 <TR_in_seconds> <total_reps_idx> <Block_Duration_in_minutes>"
                exit 1
        fi

        bar_unit="█"  # Can be changed to ▓, ▒, etc.



        # === HEADER ===
        
        echo -e "${GREEN}\nTR = $VolTR sec${NC} | ${BLUE}Block size = $duration min${NC} | ${YELLOW}Volumes/block = $volumes_per_block${NC}"
        echo -e "${RED}Total volumes = $total_reps_idx\n${NC}"
        echo "Bar units represent actual minutes of scan time"
        echo -e "\nBlock Visualization:\n"

        # === LOOP ===
        block_num=0
        for ((start_idx=0; start_idx<total_reps_idx; start_idx+=volumes_per_block)); do
                end_idx=$((start_idx + volumes_per_block - 1))
                        (( end_idx >= total_reps_idx )) && end_idx=$((total_reps_idx - 1))


                # Create a meaningful label for this block of repetitions
                label="${start_idx}_to_${end_idx}"
                block_size=$((end_idx - start_idx + 1))
                block_num=$((block_num + 1))

                # Time duration in minutes
                duration_min=$(echo "scale=2; $block_size * $VolTR / 60" | bc)
                rounded_min=$(printf "%.0f" "$duration_min")

                # Draw one BAR_UNIT per minute
                bar=$(printf "%-${rounded_min}s" "" | tr " " "$bar_unit")

                # Choose color
                if (( rounded_min < duration )); then
                        color="$YELLOW"
                else
                        color="$GREEN"
                fi

                # Output
                printf "Processing ${RED}Block Number %2d${NC}: ${color}%s${NC} [Vols %3d - %3d | ~%4.2f min]\n" \
                "$block_num" "$bar" "$start_idx" "$end_idx" "$duration_min"
                echo ""
                echo ""
                # Step 2a: Compute the mean for the current block of 60 repetitions
                3dTstat -mean -prefix block_mean_${label}.nii.gz ${input_4d_data}"[${start_idx}..${end_idx}]"
                
                #Step 2b: Subtract the baseline and divide by the baseline
                3dcalc -a block_mean_${label}.nii.gz -b baseline_image.nii.gz -expr '(a-b)/b' -prefix ratio_processed_${label}.nii.gz

                #Step 2c: Converting into percent by multiplying it by 100
                3dcalc -a ratio_processed_${label}.nii.gz -expr 'a*100' -prefix processed_${label}.nii.gz
                    
                #Add the processed image to the list
                processed_images+=("processed_${label}.nii.gz")
                
                echo ""
                echo ""
        done

        # === LEGEND ===
        echo -e "\n${GREEN}${bar_unit}${NC} = 1 minute of scan (full block)"
        echo -e "${YELLOW}${bar_unit}${NC} = partial block duration"

        timestamp=$(date +"%Y-%m-%d_%H-%M")

        #Step 3: Combine all the processed blocks into a single 4D image
        3dTcat -prefix Signal_Change_Map_premask_${duration}min_${timestamp}.nii.gz "${processed_images[@]}"

        fslmaths Signal_Change_Map_premask_${duration}min_${timestamp}.nii.gz -mas mask_mean_mc_func.nii.gz Signal_Change_Map_${duration}min_${timestamp}.nii.gz

        echo "All blocks processed and combined into final 4D image: Signal_Change_Map_${duration}min_${timestamp}.nii.gz"


        rm -f *block* processed* ratio_processed*

}