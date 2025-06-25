#!/bin/bash 


source ./bias_corrected_image.sh
source ./bias_field_correction.sh
source ./check_spikes.sh
source ./coregistration.sh
source ./data_conversion.sh
source ./fMRIAnalysis_main.sh
source ./log_execution.sh
source ./motion_correction.sh
source ./paramters_extraction.sh
source ./quality_check.sh
source ./signal_change_map.sh
source ./smoothing_using_fsl.sh
source ./temporal_snr_using_afni.sh
source ./temporal_snr_using_fsl.sh
bash ./toolbox_name.sh