#!/bin/bash

SMOOTHING_using_FSL () {
    fslmaths $1 -s 1.1774 'sm_'$1
}
