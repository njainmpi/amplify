#!/bin/sh

#Function 1
TEMPORAL_SNR_using_FSL () {
    echo "******* Computing Temporal SNR *******"
    fslmaths $1 -Tmean rG1_fsl_mean
    fslmaths $1 -Tstd rG1_fsl_std
    fslmaths rG1_fsl_mean.nii.gz -div rG1_fsl_std.nii.gz rG1_fsl_tSNR
}