#!/bin/bash

SMOOTHING_using_FSL () {
    fslmaths ${1} -s 0.24 sm_${1}
    fslmaths sm_${1} -mas $2 cleaned_sm_${1}
}
