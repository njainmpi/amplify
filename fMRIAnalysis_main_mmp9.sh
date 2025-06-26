#!/bin/bash

#Following script has been made by Naman Jain with following features included in the
#different version upgrades

##Calling all the functions that will be used in the upcoming script

#01.04.2025: Intial Script Planned, all functions called through external script

bash ./toolbox_name.sh
source ./log_execution.sh
source ./missing_run.sh
source ./folder_existence_function.sh
source ./func_parameters_extraction.sh
source ./bias_field_correction.sh
source ./check_spikes.sh
source ./coregistration.sh
source ./data_conversion.sh
source ./motion_correction.sh
source ./quality_check.sh
source ./signal_change_map.sh
source ./smoothing_using_fsl.sh
source ./temporal_snr_using_afni.sh
source ./temporal_snr_using_fsl.sh


##In order to use awk, you need to convert xlsx file to csv file

root_location="/Users/njain/Desktop/fMRI"

cd $root_location/RawData


# Read the CSV file line by line, skipping the header
awk -F ',' 'NR>1 {print $0}' "Animal_Experiments_Sequences_v4.csv" | while IFS=',' read -r col1 dataset_name project_name sub_project_name structural_name functional_name _
do
    # Trim any extra whitespace
    project_name=$(echo "$project_name" | xargs)
    
    if [[ "$project_name" == "Project_MMP9_NJ_MP" ]]; then
        export Project_Name="$project_name"
        export Sub_project_Name="$sub_project_name"
        export Dataset_Name="$dataset_name"
        export structural_run="$structural_name"
        export run_number="$functional_name"
        
        # echo $Structural_Data

        Path_Raw_Data="$root_location/RawData/$project_name/$sub_project_name"
        Path_Analysed_Data="$root_location/AnalysedData/$project_name/$sub_project_name/$Dataset_Name"
    
        # Add your further processing steps here

        datapath=$(find "$Path_Raw_Data" -type d -name "*${Dataset_Name}*" 2>/dev/null)
        # echo "$datapath"     

        echo ""
        echo ""
        echo "Dataset Currently Being Analysed is $Dataset_Name" 
        echo "from $Project_Name with Subproject $Sub_project_Name"
        echo "for Structural run number $structural_run and Functional run number $run_number" 
        echo ""
        echo ""

        if [ -d "$Path_Analysed_Data" ]; then
            echo "$Path_Analysed_Data does exist."
        else
            mkdir $Path_Analysed_Data
        fi

        cd $Path_Analysed_Data

        echo ""
        echo ""

        LOG_DIR="$datapath/Data_Analysis_log" # Define the log directory where you want to store the script.
        user=$(whoami)
        log_execution "$LOG_DIR" || exit 1


        #conversion for structural data
        FUNC_PARAM_EXTARCT $datapath/$structural_run
               
        CHECK_FILE_EXISTENCE "$Path_Analysed_Data/$structural_run""$SequenceName"
        cd $Path_Analysed_Data/$structural_run''$SequenceName

        run_if_missing "G1_cp.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$structural_run" "$datapath/$structural_run/method"
        echo "This data is acquired using $SequenceName"


        #conversion for functional data
        FUNC_PARAM_EXTARCT $datapath/$run_number

        echo $SequenceName

        CHECK_FILE_EXISTENCE "$Path_Analysed_Data/$run_number$SequenceName"
        cd $Path_Analysed_Data/$run_number''$SequenceName
    
        run_if_missing "G1_cp.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$run_number" "$datapath/$run_number/method"
        echo "This data is acquired using $SequenceName"

  
        ## Function to perform motion correction

        echo ""
        echo ""
        echo "Performing Step 1: Motion Correction"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "Motion Correction using AFNI executed on Run Number $run_number acquired using $SequenceName" || exit 1
        run_if_missing "mc_func.nii.gz" "mc_func+orig.HEAD" "mc_func+orig.BRIK" -- MOTION_CORRECTION "$MiddleVolume" G1_cp.nii.gz mc_func

        ## Function to obtain tSNR Maps
        
        echo ""
        echo ""
        echo "Performing Step 2: Obtaining Mean func, Std func and tSNR Maps"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "Temporal SNR estimated on Run Number $run_number acquired using $SequenceName" || exit 1
        run_if_missing  "tSNR_mc_func.nii.gz" "tSNR_mc_func+orig.HEAD" "tSNR_mc_func+orig.BRIK" -- TEMPORAL_SNR_using_AFNI mc_func+orig

        # Function to perform Bias Field Corrections
        
        echo ""
        echo ""
        echo "Performing Step 3: Performing N4 Bias Field Correction of mean_mc_func"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "N4Bias Field Correction on Run Number $run_number acquired using $SequenceName" || exit 1
        run_if_missing  "cleaned_mc_func.nii.gz" -- BIAS_CORRECTED_IMAGE mean_mc_func.nii.gz 100 mc_func.nii.gz
        # -b (Inpout #2 in above command) [54,3] means start with 32 points scale (equiv 20mm coil divided by 0.375mm resolution) with 3rd order b-spline
       
        ## Function to Check for Spikes
        
        echo ""
        echo ""
        echo "Performing Step 4: Checking presence of spikes in the data"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "Checked for presence of spikes in the data on Run Number $run_number acquired using $SequenceName" || exit 1
        run_if_missing "before_despiking_spikecountTC.png" -- CHECK_SPIKES mc_func.nii.gz

        ## Function to remove spikes from the data
        
        echo ""
        echo ""
        echo "Performing Step 5: Removing spikes from the data"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "Checking for Spikes and Despiking Run Number $run_number acquired using $SequenceName" || exit 1
        run_if_missing  "despike_cleaned_mc_func.nii.gz" -- DESPIKE despike_cleaned_mc_func.nii.gz cleaned_mc_func.nii.gz


        ## Function to perform Smoothing and clean the data to get it ready for estimating Signal change maps        
        
        echo ""
        echo ""
        echo "Performing Step 6: Smoothing of the data - 1 or 2 voxel smoothing"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "Smoothing using FSL executed on Run Number $run_number acquired using $SequenceName" || exit 1
        run_if_missing  "sm_despike_cleaned_mc_func.nii.gz" -- SMOOTHING_using_FSL despike_cleaned_mc_func.nii.gz mask_mean_mc_func.nii.gz


        #Function for estimating Signal Change Maps
        
        echo ""
        echo ""
        echo "Performing Step 7: Estimating Signal Change Maps"
        echo ""
        echo ""
        log_function_execution "$LOG_DIR" "Signal Change Map created for Run Number $run_number acquired using $SequenceName" || exit 1

        if [[ "$SequenceName" == *"functionalEPI"* ]]; then
            run_if_missing "Signal_Change_Map.nii.gz" -- \
            SIGNAL_CHANGE_MAPS cleaned_sm_despike_cleaned_mc_func.nii.gz "$datapath/$run_number" 10 10
        elif [[ "$SequenceName" == *"FLASH"* ]]; then
            run_if_missing "Signal_Change_Map.nii.gz" -- \
            SIGNAL_CHANGE_MAPS mc_func.nii.gz 5 12 "$datapath/$run_number" 5 5 mean_mc_func.nii.gz
        else
            echo "Unknown sequence type: $SequenceName — skipping SIGNAL_CHANGE_MAPS."
        fi


        # ##Function to perform coregistration and displaying Signal Change Maps with anatomical underlay

        # echo ""
        # echo ""
        # echo "Performing Step 8: Estimating Transformation Matrix and Coregistration of the data "
        # echo ""
        # echo ""
        
        #     # Loop until user is satisfied
        # while true; do
        #     COREGISTRATION cleaned_N4_mean_mc_func.nii.gz anatomy.nii.gz

        #     read -p "Are you happy with the registration? (yes/no): " response

        # if [[ "$response" == "yes" || "$response" == "y" ]]; then
        #     echo "Proceeding to the next step..."
        #     break
        # else
        #     echo "Re-running coregistration..."
        # fi
        # done

    fi
done



