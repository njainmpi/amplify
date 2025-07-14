#!/bin/sh


#06th June 2024: $$Naman Jain$$ This function is created to converting Bruker and Dicom file to NIFTI


#Function 1
DICOM_to_NIFTI () {
    
    NoOfSlices="`grep -A1 'PVM_SPackArrNSlices=( 1 )' $1 | grep -v -e 'PVM_SPackArrNSlices=( 1 )' -e '--'`"
    #here the grep command will read the no of slices from the Bruker raw file acquired during acquisiton
    NoOfRepetitions=$(awk '/PVM_NRepetitions=/ {print substr($0,21,3)}' $1)
    #here the awk will look at the number of slices acquired using the information located in the methods file
    NoOfEchoImages=$(awk '/PVM_NEchoImages=/ {print substr($0,20,2)}' $1)
    TotalNumberOfRepetitions=$(($NoOfRepetitions * $NoOfEchoImages))
    VolumeTR=$(awk '/PVM_RepetitionTime=/ {print substr($0,23,5)}' $1)
    #here the awk will look at the number of slices acquired using the information located in the methods file
    
    echo "DICOM file used to get NIFTI file using AFNI"
    
    if [ $NoOfEchoImages == 1 ]; then 
        to3d -prefix $2 -time:zt $NoOfSlices $TotalNumberOfRepetitions $VolumeTR alt+z $3
    else
        to3d -prefix $2 -time:tz $TotalNumberOfRepetitions $NoOfSlices $VolumeTR alt+z $3
    fi

    3dAFNItoNIFTI $2*.BRIK
    cp $2.nii G1_cp.nii.gz
    fslhd G1_cp.nii.gz > NIFTI_file_header_info.txt
}

#Function 2
BRUKER_to_NIFTI () {
echo "Using brkraw to convert Bruker format to NIFTI"
    format="Bruker"
    brkraw tonii $1/ -s $2

    if grep -q "PVM_NEchoImages" "$3"; then
        NoOfEchoImages=$(awk '/PVM_NEchoImages=/ {print substr($0,20,2)}' $3)

        echo "No of echo images $NoOfEchoImages"
        if [ $NoOfEchoImages == 1 ]; then 
            cp *$2* G1_cp.nii.gz
        else
            fslmerge -t $2'_combined_images' *$2*
            cp $2'_combined_images'* G1_cp.nii.gz
        fi
    else
        echo "No of Echoes not present"
        cp *$2* G1_cp.nii.gz
    fi

    fslhd G1_cp.nii.gz > NIFTI_file_header_info.txt
}

