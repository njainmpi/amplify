#!/bin/bash


TEMPORAL_SNR_using_AFNI () {
    3dTstat -mean -prefix 'mean_'$1 $1
    3dTstat -stdev -prefix 'std_'$1 $1
    3dcalc -a 'mean_'$1 -b 'std_'$1 -expr 'a/b' -prefix 'tSNR_'$1
    3dAFNItoNIFTI "mean_$1" && gzip "mean_${1/+orig/.nii}"
    3dAFNItoNIFTI "tSNR_$1" && gzip "tSNR_${1/+orig/.nii}"

}
