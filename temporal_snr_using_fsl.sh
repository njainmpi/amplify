#!/bin/sh

#Function 1
TEMPORAL_SNR_using_FSL () {
    echo "******* Computing Temporal SNR *******"
    fslmaths ${1} -Tmean mean_${1}
    fslmaths ${1} -Tstd std_${1}
    fslmaths mean_${1} -div std_${1} tSNR_mc_func.nii.gz
}