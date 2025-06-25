#!/bin/bash

COREGISTRATION () {


    echo "Please save your transformation matrix by the name 'anatomy_to_epi_mean.txt' "
    
    # Step 1: Launch ITK-SNAP to align manually
    itksnap -g ${1} -o ${2}

    # Step 2: Apply the transformation
    antsApplyTransforms -d 3 -i ${1} -r ${2} -o epi_mean_to_anatomy.nii.gz -t [anatomy_to_epi_mean.txt,1] -n BSpline

    # Step 3: Visualize result in FSLeyes
    fsleyes epi_mean_to_anatomy.nii.gz ${1}



}
