#!/bin/bash


Static_Map () {


    input_4d_data=$1
    local base_start=$2
    local base_end=$3
    local sig_start=$4
    local sig_end=$5  


    base_label="${base_start}_to_${base_end}"
    sig_label="${sig_start}_to_${sig_end}"    

    echo $base_start $base_end $sig_start $sig_end

    3dTstat -mean -prefix baseline_image_${base_label}.nii.gz ${input_4d_data}"[${base_start}..${base_end}]"
    3dTstat -mean -prefix signal_image_${sig_label}.nii.gz ${input_4d_data}"[${sig_start}..${sig_end}]"

    #Step 2b: Subtract the baseline and divide by the baseline
    3dcalc -a signal_image_static.nii.gz -b baseline_image_static.nii.gz -expr '(a-b)/b' -prefix signal_processed.nii.gz

    #Step 2c: Converting into percent by multiplying it by 100
    3dcalc -a signal_processed.nii.gz -expr 'a*100' -prefix Static_SCM_${base_label}_${sig_label}.nii.gz

}

