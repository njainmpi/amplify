#!/bin/bash



COREGISTRATION_ROI () {
    
   
    # Help Message
    if [[ $# -lt 3 || $1 == "--help" || $1 == "-h" ]]; then
        echo "Usage: COREGISTRATION_ROI <input_roi> <ref_image> <transform>"
        echo ""
        echo "Description:"
        echo "  Applies a spatial transformation to an ROI image, thresholds the transformed ROI,"
        echo "  and extracts the mean time series from a functional image using the cleaned ROI mask."
        echo ""
        echo "Positional Arguments:"
        echo "  input_roi    = Path to ROI image marked on the structural image"
        echo "  ref_image    = Mean functional image where you want to move your ROI"
        echo "  transform    = Transformation matrix (e.g., from ANTs registration)"
        echo ""
        echo "Example:"
        echo "  COREGISTRATION_ROI roi.nii.gz mean_func.nii.gz transformation_matrix.txt"
        return
    fi

    # Assign input parameters
    input_roi="$1"
    ref_image="$2"
    transform="$3"

    input_name=$(basename "${input_roi%.nii.gz}")
    ref_name="${ref_image%.nii.gz}"

    # Apply transform
    antsApplyTransforms \
        -i "$input_roi" \
        -r "$ref_image" \
        -o "${input_name}_to_${ref_name}.nii.gz" \
        -t "$transform" \
        -n BSpline


    echo -e "Your \033[31mCoregistered ROI\033[0m is saved by the name \033[32m${input_name}_to_${ref_name}.nii.gz\033[0m"

    # Threshold transformed ROI
    fslmaths "${input_name}_to_${ref_name}.nii.gz" -thr 0.4 -uthr 1.1 "${input_name}_to_${ref_name}_cleaned.nii.gz"

    # Extract time series
    fslmeants -i mc_func.nii.gz -m "${input_name}_to_${ref_name}_cleaned.nii.gz" -o "ts_${input_name}.txt"
    echo -e "Your \033[31mTime Series\033[0m is saved by the name \033[32mts_${input_name}.txt. \033[0m"

    python $ts_roi_python_script ts_${input_name}.svg ts_${input_name}.txt
    echo -e "Your \033[31mVectorised TS Graph\033[0m is saved by the name \033[32mts_${input_name}.svg. \033[0m"
}


COREGISTRATION_UPSAMPLING () {


        echo ""
        echo "Description:"
        echo "  This function performs manual alignment (coregistration) between a mean functional image"
        echo "  and a high-resolution anatomical image using ITK-SNAP for visualization and ANTs for transformation."
        echo ""
        echo "Steps:"
        echo "  1. Opens ITK-SNAP with 'mean_func' as the background and 'anatomy.nii.gz' as the overlay image."
        echo "     You are expected to perform a manual rigid alignment and save the transformation matrix"
        echo "     as 'anatomy_to_epi_mean.txt'."
        echo ""
        echo "  2. Applies the saved transformation using ANTs to bring 'mean_func' into anatomical space."
        echo "     Output file is saved as 'epi_mean_to_anatomy.nii.gz'."
        echo ""
        echo "  3. Optionally visualizes the result in FSLeyes, overlaid with another image if provided"
        echo "     as a positional argument to this function."
        echo ""
        echo "Notes:"
        echo "  - Ensure 'mean_func' and 'anatomy.nii.gz' are in your working directory."
        echo "  - The transformation matrix must be saved manually in ITK-SNAP as 'anatomy_to_epi_mean.txt'."
        echo "  - This function uses ANTs (antsApplyTransforms) and FSLeyes. Ensure they are installed and available."
        echo ""


    # Help Message
    if [[ "$1" == "--help" || "$1" == "-h" ]]; then
        echo "Usage: COREGISTRATION_UPSAMPLING"
        echo ""
        echo "Description:"
        echo "  This function performs volume-wise upsampling of a 3D/ 4D functional image by:"
        echo "  applying an inverse spatial transform to bring each volume into anatomical space."
        echo ""
        echo "Positional Arguements (can be modified inside the function):"
        echo "  input_4D     = Signal_Change_Map.nii.gz (4D functional ROI image)"
        echo "  ref_image    = ../11anatomy/G1_cp.nii.gz (anatomical reference image)"
        echo "  transform    = anatomy_to_func.txt (transform matrix, used in inverse mode)"
        echo ""
        echo "Output:"
        echo "  - Coregistered 3D/ 4D image: Coregistered_SCM.nii.gz"
        echo ""
        echo "Note:"
        echo "  Ensure FSL and ANTs are sourced and available in your environment."
        return
    fi

    # Function body
    input_4D="$1"
    ref_image="$2"
    transform="$3"

    n_vols=$(fslval "$input_4D" dim4)

    for ((i=0; i<n_vols; i++)); do
        volname=$(printf "vol%04d" "$i")
       
        fslroi $input_4D ${volname}.nii.gz $i 1

        antsApplyTransforms \
            -i "${volname}.nii.gz" \
            -r "$ref_image" \
            -o "coreg_${volname}_anatomy.nii.gz" \
            -t "[$transform,1]" \
            -n BSpline
    done

    fslmerge -t Coregistered_SCM.nii.gz coreg_vol*

    echo -e "Your \033[31mCoregistered Signal Change Map\033[0m is saved by the name \033[32mCoregistered_SCM.nii.gz. \033[0m"
    rm -rf coreg_vol0*
    rm -rf vol0*
}
