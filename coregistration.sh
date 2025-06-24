#!/bin/bash

COREGISTRATION () {


    echo "Please save your transformation matrix by the name 'anatomy_to_epi_mean.' "
    
    # Step 1: Launch ITK-SNAP to align manually
    itksnap -g mc_func.nii.gz -o anatomy.nii.gz

    # Step 2: Apply the transformation
    antsApplyTransforms -d 3 -i ${mc_func.nii.gz} -r ${anatomy.nii.gz} -o epi_mean_to_anatomy.nii.gz -t [anatomy_to_epi_mean.txt,1] -n BSpline

    # Step 3: Visualize result in FSLeyes
    fsleyes epi_mean_to_anatomy.nii.gz mc_func.nii.gz



}




# Loop until user is satisfied
while true; do
    COREGISTRATION

    read -p "Are you happy with the registration? (yes/no): " response

    if [[ "$response" == "yes" || "$response" == "y" ]]; then
        echo "Proceeding to the next step..."
        break
    else
        echo "Re-running coregistration..."
    fi
done

# Add the next steps here
echo "Continuing with the rest of the pipeline..."
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      