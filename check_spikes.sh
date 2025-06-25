#!/bin/bash


CHECK_SPIKES () {
     
    3dToutcount -automask -fraction -polort 3 -legendre $1 > before_despiking_spikecountTC.1D
    fsl_tsplot -i before_despiking_spikecountTC.1D -o before_despiking_spikecountTC -t 'spike count' -x $Runid'_'$FunctionalRunName -y 'fraction of voxels' -h 650 -w 1800
}


DESPIKE () {

    3dDespike -prefix $1 -ssave rest_spikemap.nii.gz -NEW $2 -localedit
    3dToutcount -automask -fraction -polort 3 -legendre $1 > after_despiking_spikecountTC.1D
    fsl_tsplot -i after_despiking_spikecountTC.1D -o after_despiking_spikecountTC -t 'spike count' -x $Runid'_'$FunctionalRunName -y 'fraction of voxels' -h 650 -w 1800
}
