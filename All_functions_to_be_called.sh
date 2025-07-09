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