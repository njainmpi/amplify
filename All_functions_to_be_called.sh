#!/bin/sh

#07.11.2024: All files with functions located here

source ../fMRI_analysis_pipeline/data_conversion_function.sh #converting data from either Bruker or Dicom format to NIFTI format
source ../fMRI_analysis_pipeline/folder_existence_function.sh #check if folder is present or not
source ../fMRI_analysis_pipeline/motion_correction_function.sh #perform motion correction using AFNI
source ../fMRI_analysis_pipeline/temporal_SNR_spikes_smoothing_function.sh #check presence of spikes, peforms smoothing using either AFNI or NIFTI, caclulates temporal SNR
source ../fMRI_analysis_pipeline/time_series_function.sh
source ../fMRI_analysis_pipeline/activation_maps.sh # to map areas of activation using AFNI, also generates signal change maps
source ../fMRI_analysis_pipeline/outlier_count.sh #14.08.2024 new function to perfom slice timing correction and outlier estimate before and after slice timing correction
source ../fMRI_analysis_pipeline/video_making.sh #19.08.2024 new function to make videos of the signal change maps
source ../fMRI_analysis_pipeline/func_parameters_extraction.sh #07.11.2024 new function to extract parameters for fMRI analysis
source ../fMRI_analysis_pipeline/bash_log_create.sh
source ../fMRI_analysis_pipeline/Signal_Change_Map.sh
source ../fMRI_analysis_pipeline/missing_run.sh
bash ../fMRI_analysis_pipeline/toolbox_name.sh