#!/bin/bash


CHECK_SPIKES () {
    3dToutcount -automask -fraction -polort 3 -legendre $1 > spikecountTC.1D
    fsl_tsplot -i spikecountTC.1D -o spikecountTC -t 'spike count' -x $Runid'_'$FunctionalRunName -y 'fraction of voxels' -h 650 -w 1800
}