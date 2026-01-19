#!/bin/bash

# represents main section of the data analysis pipeline
## respresents subsections of the main section
### represents individual steps within the subsections
#### represents comments for better understanding of the code

# ============================================================
# Terminal safety
# ============================================================
reset_terminal() {
  tput sgr0
  stty sane
}
trap reset_terminal EXIT

# ==============================================================================================================================================================================================================
# Importing in-house made functions to be used throughout the analysis pipeline
# ==============================================================================================================================================================================================================

## Prime functions

bruker_to_nifti(){
    
    local in_path="$1"
    local scan_number="$2"
    local out_name="${3}"

    scan_dir="$in_path/$scan_number"
    method_file="$scan_dir/method"

    brkraw tonii $in_path/ -s $scan_number

    if grep -q "PVM_NEchoImages" "$method_file"; then
        NoOfEchoImages=$(awk -F= '/PVM_NEchoImages/ {gsub(/[^0-9]/,"",$2); print $2; exit}' "$method_file")

        echo "No of echo images $NoOfEchoImages"
        if [ $NoOfEchoImages == 1 ]; then 
            cp *${scan_number}* G1_cp.nii.gz
        else
            fslmerge -t ${scan_number}'_combined_images' *${scan_number}*
            cp ${scan_number}'_combined_images'* G1_cp.nii.gz
        fi
    else
        echo "No of Echoes not present"
        cp *${scan_number}* G1_cp.nii.gz
    fi

    echo ""
    echo "Fixing orientation to LPI"
    echo ""
    3dresample -orient LPI -inset G1_cp.nii.gz -prefix G1_cp_resampled.nii.gz -overwrite

    cp G1_cp_resampled.nii.gz ${out_name}
    fslhd ${out_name} > NIFTI_file_header_info.txt
}

masking_file(){

    local input_file="${1}"
    local out_file="${2}"
    local mask_name="mask_${input_file}"

    if [[ ! -f "$mask_name" ]]; then
        echo "❌ Mask file $mask_name not found. Please create a mask"
        fsleyes "$input_file"
        return 1
    else
        echo "✔ Mask file found: mask_${input_file}"
    fi

    fslmaths ${input_file} -mas ${mask_name} ${out_file}    
}

extract_middle_volume(){

    local input_file="${1}"
    local out_file="${2}"
    local middle_vol="${3}"
    fslroi ${input_file} ${out_file} ${middle_vol} 1
}

motion_correction () {

    local base_image="${1}"
    local input_image="${2}"
    local output_prefix="${3}"

    3dvolreg -prefix ${output_prefix} -base ${base_image} -verbose -1Dfile motion.1D -1Dmatrix_save mats -linear -float -maxdisp1D rmsabs.1D ${input_image}
    3dAFNItoNIFTI ${output_prefix}'+orig.BRIK'
    if [ -f ${output_prefix}.nii ]; then
        gzip -1 ${output_prefix}.nii
    else
        echo "gzip exists"
    fi
    # 1Dplot -volreg -sep motion.1D
    1dplot -xlabel Time -ylabel "Translation (mm)" -title "3dvolreg translations" -dx 1 -jpgs 640x144 rest_translation 'motion.1D[3..5]'
    1dplot -xlabel Time -ylabel "Rotations (degree)" -title "3dvolreg rotations" -dx 1 -jpgs 640x144 rest_rotation 'motion.1D[0..2]'

    rm -f ${3}'+orig.BRIK' ${3}'+orig.HEAD'
}

tSNR () {
    
    local input_file="${1}"

    echo "******* Computing Temporal SNR *******"

    fslmaths ${input_file} -Tmean mean_${input_file}
    fslmaths ${input_file} -Tstd std_${input_file}
    fslmaths mean_${input_file} -div std_${input_file} tSNR_${input_file}
    rm -f mean_${input_file} std_${input_file}
}

signal_change_map(){

    local input="${1}"
    local base_start="${2}"
    local base_end="${3}"
    local sig_start="${4}"
    local sig_end="${5}"
    local out_psc="${6}"
    
    dim_inplane=$(fslhd ${input} | awk '/^pixdim4/ {print $2}')
    
    3dTstat -mean -prefix "signal_image_${sig_start}_to_${sig_end}.nii.gz" "${input}[${sig_start}..${sig_end}]"
    3dTstat -mean -prefix "baseline_image_${base_start}_to_${base_end}.nii.gz" "${input}[${base_start}..${base_end}]"

    echo "Estimating Signal Change Map"
    fslmaths "signal_image_${sig_start}_to_${sig_end}.nii.gz" -sub "baseline_image_${base_start}_to_${base_end}.nii.gz" -div "baseline_image_${base_start}_to_${base_end}.nii.gz" -mul 100 "${out_psc}.nii.gz" #### Estimates Static SCM ####
    
    echo "Estimating Signal Change Time Series"
    fslmaths ${input} -sub "baseline_image_${base_start}_to_${base_end}.nii.gz" -div "baseline_image_${base_start}_to_${base_end}.nii.gz" -mul 100 "time_series_${out_psc}.nii.gz" #### Estimates SC Time Series ####


    rm -f "signal_image_${sig_start}_to_${sig_end}.nii.gz" "baseline_image_${base_start}_to_${base_end}.nii.gz"
}

coregistration_afni() {

  # -------- Arguments (named via env or positional) --------
  local input_file1="$1"
  local input_file2="$2"
  local reference_file="$3"
  local output_file1="$4"
  local output_file2="$5"

  local estimate_affine="${6:-true}"
  local apply_affine="${7:-true}"
  local affine_mat="${8:-mean_func_struct_aligned.aff12.1D}"

  # -------- Safety checks --------
  if [[ ! -f "$reference_file" ]]; then
    echo "[ERROR] reference_file must be provided" >&2
    return 1
  fi

  # -------- STEP 1: Estimate affine --------
  if [[ "$estimate_affine" == "true" ]]; then

    [[ -z "$input_file1" ]] && { echo "[ERROR] input_file1 required for affine estimation"; return 1; }
    [[ -z "$output_file1" ]] && { echo "[ERROR] output_file1 required for affine estimation"; return 1; }

    3dAllineate -input "$input_file1" -base "$reference_file" -prefix "$output_file1" -1Dmatrix_save "$affine_mat" -cost lpa -twopass -1Dparam_save params.1D -verb

    if [[ $? -ne 0 ]]; then
      echo "[ERROR] Affine estimation failed" >&2
      return 1
    fi

    echo "[OK] Affine estimated and saved → $affine_mat"
    echo "[OK] Coregistered image (step 1) → $output_file1"
  fi

  # -------- STEP 2: Apply affine --------
  if [[ "$apply_affine" == "true" ]]; then

    [[ -z "$input_file2" ]] && { echo "[ERROR] input_file2 required for affine application"; return 1; }
    [[ -z "$output_file2" ]] && { echo "[ERROR] output_file2 required for affine application"; return 1; }
    [[ ! -f "$affine_mat" ]] && { echo "[ERROR] Affine matrix not found: $affine_mat"; return 1; }

    3dAllineate -input "$input_file2" -base "$reference_file" -master "$reference_file" -prefix "$output_file2" -1Dmatrix_apply "$affine_mat" -final linear -verb

    if [[ $? -ne 0 ]]; then
      echo "[ERROR] Applying affine failed" >&2
      return 1
    fi

    echo "[OK] Affine applied → $output_file2"
  fi
}

smooth_movavg() {
  local in_file="$1"
  local out_file="$2"
  local win_sec_duration="$3"
  local tr="$4"

  python3 <<EOF
import numpy as np
import nibabel as nib

win_sec_duration = float("$win_sec_duration")
tr = float("$tr")
in_file = "$in_file"
out_file = "$out_file"

win = max(1, int(round(win_sec_duration / tr)))

def moving_average_1d(x, win):
    k = np.ones(win, dtype=float) / win
    xpad = np.pad(x, (win//2, win-1-win//2), mode='edge')
    return np.convolve(xpad, k, mode='valid')

img = nib.load(in_file)
data = img.get_fdata()
T = data.shape[-1]
flat = data.reshape(-1, T)
sm = np.vstack([moving_average_1d(ts, win) for ts in flat]).reshape(data.shape)

nib.Nifti1Image(sm, img.affine, img.header).to_filename(out_file)
print(f"Wrote: {out_file}  (tr={tr}s, window={win_sec_duration}s => {win} vols)")
EOF
}


roi_analysis() {

  local func_in_file="$1"
  local roi_file="$2"
  local base_start="$3"
  local base_end="$4"
  local sig_start="$5"
  local sig_end="$6"
  local TR="$7"

  local tc_file="time_course_${roi_file%.nii.gz}.txt"
  local psc_file="PSC_time_series_${roi_file%.nii.gz}.txt"
  local graph_file="PSC_Time_Series_${roi_file%.nii.gz}.svg"

  echo "[INFO] Processing ROI: $roi_file" >&2

  # ---- Time course extraction ----
  fslmeants -i "$func_in_file" -m "$roi_file" -o "$tc_file" || return 1
  echo "[OK] Time course saved → $tc_file" >&2

  # ---- Baseline mean (FIXED: 0-based → 1-based indexing) ----
  baseline=$(awk -v s="$base_start" -v e="$base_end" '
    NR >= (s+1) && NR <= e {sum+=$1; n++}
    END { if (n>0) print sum/n; else print "nan" }' "$tc_file")

  if [[ "$baseline" == "nan" || "$baseline" == "0" || -z "$baseline" ]]; then
    echo "[ERROR] Invalid baseline ($baseline) for $roi_file" >&2
    return 1
  fi

  echo "[DEBUG] Baseline = $baseline" >&2

  # ---- Percent Signal Change ----
  awk -v b="$baseline" '{print (($1-b)/b)*100}' "$tc_file" > "$psc_file"
  echo "[OK] PSC calculated → $psc_file" >&2

  # ---- Mean PSC in signal window (FIXED indexing) ----
  mean_signal=$(awk -v s="$sig_start" -v e="$sig_end" '
    NR >= (s+1) && NR <= e {sum+=$1; n++}
    END { if (n>0) print sum/n; else print "nan" }' "$psc_file")

  # ---- LIVE TERMINAL PLOT ----
  echo "[INFO] Live PSC plot (terminal)" >&2
  gnuplot << EOF
set terminal dumb size 120,30
set title "LIVE PSC – $roi_file"
set xlabel "Time (minutes)"
set ylabel "PSC (%)"
set grid
plot "$psc_file" using (\$0*$TR/60.0):1 with lines title "PSC"
EOF

  # ---- SVG PLOT ----
  gnuplot << EOF
set terminal svg size 900,400
set output "$graph_file"
set title "Percent Signal Change – $roi_file"
set xlabel "Time (minutes)"
set ylabel "Signal Change (%)"
set grid

set style rect fc rgb "green" fs solid 0.2 noborder
set object rect from ($base_start*$TR/60.0), graph 0 \
                 to ($base_end*$TR/60.0), graph 1

set style rect fc rgb "blue" fs solid 0.2 noborder
set object rect from ($sig_start*$TR/60.0), graph 0 \
                 to ($sig_end*$TR/60.0), graph 1

plot "$psc_file" using (\$0*$TR/60.0):1 with lines lw 2 title "PSC"
EOF

  echo "[OK] SVG saved → $graph_file" >&2

  # ---- RETURN VALUE ----
  echo "$mean_signal"
}

## Helper functions

func_param_extract() {

  local scan_dir="$1"
  local export_env="${2:-true}"

  local acqp_file="$scan_dir/acqp"
  local method_file="$scan_dir/method"

  # -----------------------------
  # File checks
  # -----------------------------
  if [[ ! -f "$acqp_file" || ! -f "$method_file" ]]; then
    echo "❌ acqp or method file not found in $scan_dir" >&2
    return 1
  fi

  # -----------------------------
  # Extract Sequence Name
  # -----------------------------
  SequenceName=$(sed -nE \
    's/.*ACQ_protocol_name=\([[:space:]]*64[[:space:]]*\)[[:space:]]*<([^>]+)>.*/\1/p' \
    "$acqp_file")

  # -----------------------------
  # Helper to extract numeric values
  # -----------------------------
  extract_val() {
    local pattern="$1"
    local file="$2"
    grep -oE "$pattern" "$file" | grep -oE '[0-9]+'
  }

  NoOfRepetitions=$(extract_val '##\$PVM_NRepetitions=[[:space:]]*[0-9]+' "$method_file")
  TotalScanTime=$(extract_val '##\$PVM_ScanTime=[[:space:]]*[0-9]+' "$method_file")

  Baseline_TRs=$(extract_val 'PreBaselineNum=[[:space:]]*[0-9]+' "$method_file")
  StimOn_TRs=$(extract_val 'StimNum=[[:space:]]*[0-9]+' "$method_file")
  StimOff_TRs=$(extract_val 'InterStimNum=[[:space:]]*[0-9]+' "$method_file")
  NoOfEpochs=$(extract_val 'NEpochs=[[:space:]]*[0-9]+' "$method_file")

  # -----------------------------
  # Derived values
  # -----------------------------
  VolTR_msec=""
  VolTR=""
  MiddleVolume=""

  if [[ -n "$NoOfRepetitions" && -n "$TotalScanTime" ]]; then
    VolTR_msec=$(echo "scale=4; $TotalScanTime / $NoOfRepetitions" | bc)
    VolTR=$(echo "scale=4; $VolTR_msec / 1000" | bc)
    MiddleVolume=$(echo "scale=2; $NoOfRepetitions / 2" | bc)
  fi

  # -----------------------------
  # Export to environment (optional)
  # -----------------------------
  if [[ "$export_env" == "true" ]]; then
    export SequenceName
    export NoOfRepetitions
    export TotalScanTime
    export VolTR_msec
    export VolTR
    export Baseline_TRs
    export StimOn_TRs
    export StimOff_TRs
    export NoOfEpochs
    export MiddleVolume
  fi

  # -----------------------------
  # Print key=value (like return)
  # -----------------------------
  cat <<EOF
SequenceName=$SequenceName
NoOfRepetitions=$NoOfRepetitions
TotalScanTime=$TotalScanTime
VolTR_msec=$VolTR_msec
VolTR=$VolTR
Baseline_TRs=$Baseline_TRs
StimOn_TRs=$StimOn_TRs
StimOff_TRs=$StimOff_TRs
NoOfEpochs=$NoOfEpochs
MiddleVolume=$MiddleVolume
EOF
}

is_step_enabled() {
  local key="$1"
  echo "$SELECTED_STEPS" | grep -qw "$key"
}


# ======================================================================
# Terminal safety (CRITICAL)
# ======================================================================

# ======================================================================
# Asking what steps need to be executed
# ======================================================================

PIPELINE_STEPS=(
  "bruker_to_nifti|Bruker → NIfTI|on"
  "motion_correction|Motion Correction|on"
  "smooth_movavg|Temporal Smoothing (MovAvg)|on"
  "spatial_smoothing|Spatial Smoothing|on"
  "masking_file|Masking|on"
  "tSNR|Temporal SNR Estimation|on"
  "signal_change_map|Signal Change Map|on"
  "roi_analysis|ROI Analysis|on"
  "coregistration_afni|Coregistration (AFNI)|off"
)

DIALOG_ARGS=()
for step in "${PIPELINE_STEPS[@]}"; do
  IFS="|" read -r key label default <<< "$step"
  DIALOG_ARGS+=("$key" "$label" "$default")
done

printf '\e[?1l\e>'   # reset keypad mode before dialog

SELECTED=$(dialog \
  --clear \
  --backtitle "Neuroimaging Pipeline" \
  --title "Pipeline Step Selection" \
  --checklist "Select preprocessing steps to run:" \
  20 80 12 \
  "${DIALOG_ARGS[@]}" \
  3>&1 1>&2 2>&3)

STATUS=$?
reset_terminal

if [[ $STATUS -ne 0 ]]; then
  echo "Pipeline selection cancelled."
  exit 1
fi

SELECTED_STEPS=$(echo "$SELECTED" | tr -d '"')

echo "✔ Selected pipeline steps:"
for step in $SELECTED_STEPS; do
  echo "• $step"
done
echo

# ======================================================================
# Selecting raw data folder
# ======================================================================

root_location="/Users/njain/Desktop/test"

printf '\e[?1l\e>'   # reset keypad mode before dialog

in_path=$(dialog \
  --clear \
  --no-shadow \
  --title "Select Raw Data Folder" \
  --dselect "$root_location/" 18 90 \
  3>&1 1>&2 2>&3)

STATUS=$?
reset_terminal

if [[ $STATUS -ne 0 ]]; then
  echo "❌ Folder selection cancelled."
  exit 1
fi

# Normalize path
[[ "$in_path" != */ ]] && in_path="${in_path}/"

if [[ ! -d "$in_path" ]]; then
  echo "❌ Invalid directory selected."
  exit 1
fi

echo "✔ Raw data folder selected:"
echo "  $in_path"
echo

## Final summary about the functions selected here

echo "======================================"
echo "✔ Pipeline Configuration Summary"
echo "======================================"
echo "Raw Data Path : $in_path"
echo
echo "Selected Steps:"
echo "--------------------------------------"
for step in $SELECTED_STEPS; do
  echo "• $step"
done
echo "======================================"


# ==============================================================================================================================================================================================================
# Creating and assigning directories to the variable
# ==============================================================================================================================================================================================================

## Using dialog box here, you enter the scan parameters and the values needed for analysis

FORM_VALUES=$(dialog \
  --clear \
  --title "Scan Parameters" \
  --form "Enter scan parameters:" \
  22 110 7 \
  "Structural Scan Number:"                                     1  1 "10"   1  50 10 0 \
  "Functional Scan Number:"                                     2  1 "1"    2  50 10 0 \
  "Window Duration (in vols):"                                  3  1 "100"  3  50 10 0 \
  "Temporal Smooth (in vols):"                                  4  1 "60"   4  50 10 0 \
  "Spatial Smoothing for Functional Data (num of voxels):"      5  1 "300"  5  50 10 0 \
  "Spatial Smoothing for Coregistered Data (num of voxels):"    6  1 "300"  6  50 10 0 \
  3>&1 1>&2 2>&3)



#### below piece of code is added to ensure terminal safety even if dialog is exited without entering values i.e. being cancalled ####
STATUS=$?
reset_terminal

if [[ $STATUS -ne 0 ]]; then
  echo "❌ Parameter entry cancelled."
  exit 1
fi

## Parse dialog output CORRECTLY (multi-line safe) ##

mapfile -t FORM_ARRAY <<< "$FORM_VALUES" # mapfile reads lines into an array and works only with bash v4.0 and above

struct_scan_number="${FORM_ARRAY[0]}"
func_scan_number="${FORM_ARRAY[1]}"
win_vol="${FORM_ARRAY[2]}"
temp_smooth_factor="${FORM_ARRAY[3]}"
func_smooth_factor="${FORM_ARRAY[4]}"
coreg_smooth_factor="${FORM_ARRAY[5]}"


#### Getting subject ID ####
subject_id=$(awk '/^<.*>$/ {gsub(/[<>]/,""); print; exit}' "${in_path}/subject")


echo "Data is analysed for ${subject_id} for Functional Scan Number: ${func_scan_number} and Structural Scan Number: ${struct_scan_number} and a termporal smoothing of ${temp_smooth_factor} vols and spatial smoothing of ${func_smooth_factor} for functional data and ${coreg_smooth_factor} for coregistered data is applied."
## Setting up path and assigning locations to variables

path_raw_struct="${in_path}/${struct_scan_number}"  # Path for raw structural data
path_raw_func="${in_path}/${func_scan_number}"  # Path for raw functional data

sequence_name_struct=$(func_param_extract "$path_raw_struct" false | grep '^SequenceName=' | cut -d= -f2)
sequence_name_func=$(func_param_extract "$path_raw_func" false | grep '^SequenceName=' | cut -d= -f2)

## Setting up path for storing analysed data overall, strucutral data, functional data and processed functional data

#### Below piece of code replaced RawData with AnalysedData in the input path to create analysed data path ####
analysed_path="$(dirname "${in_path/RawData/AnalysedData}")/$subject_id"

#### Structural and Functional data analysed data directories ####
analysed_struct_dir=${analysed_path}/${struct_scan_number}${sequence_name_struct}
analysed_func_dir=${analysed_path}/${func_scan_number}${sequence_name_func}

#### processed functional data ####
timestamp="$(date +"%Y_%m_%d_%H%M%S")"
user="$(whoami)"
folder_created="${timestamp}_${user}"
analysis_func_subdir="${analysed_func_dir}/${folder_created}"
if [[ -e "$analysis_func_subdir" ]]; then
  echo "❌ Analysis directory already exists:"
  echo "   $analysis_func_subdir"
  exit 1
fi
mkdir -p "$analysis_func_subdir"


# src_dir = Path.cwd()

#### PSQL-style table header ####

printf "+----+----------------------------------------------+---------------------------+------------------------------------------------------------+\n"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "ID" "Purpose" "Variable Name" "Value"
printf "+----+----------------------------------------------+---------------------------+------------------------------------------------------------+\n"


printf "| %-2s | %-44s | %-25s | %-58s |\n" "a" "Location of all Raw Data"                      "root_location"        "$root_location"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "b" "Location of Raw Data to be analysed"           "in_path"              "$in_path"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "c" "Location of Raw Structural Data"               "path_raw_struct"      "$path_raw_struct"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "d" "Location of Raw Functional Data"               "path_raw_func"        "$path_raw_func"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "e" "Location of analysed structural data"          "analysed_struct_dir"  "$analysed_struct_dir"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "f" "Location of analysed functional data"          "analysed_func_dir"    "$analysed_func_dir"
printf "| %-2s | %-44s | %-25s | %-58s |\n" "g" "Location of final processed functional data"   "analysis_func_subdir" "$analysis_func_subdir"

printf "+----+----------------------------------------------+---------------------------+------------------------------------------------------------+\n"


echo "For your information: all generated files after processing stored in the directory will be utilising the below nomenclature:"

echo "After applying motion correction, file will be saved as mc_input_file"
echo "After applying temporal smoothing, file will be saved as ts_input_file"
echo "After applying spatial smoothing, file will be saved as sm_input_file"
echo "After applying temporal smoothing, file will be saved as ts_input_file"
echo "After applying temporal SNR, file will be saved as tsnr_input_file"



# ==============================================================================================================================================================================================================
# Converting Structural Data into NIFTI and Processing Structural Data
# ==============================================================================================================================================================================================================

## Converting Structural Data into NIFTI ##
cd ${analysed_struct_dir}

if is_step_enabled "bruker_to_nifti"; then
    echo "Converting Bruker to NIFTI: Structural Data"
    
    nifti_file_name="struct.nii.gz"
    if [[ -f "$nifti_file_name" ]]; then
        echo "✔ NIfTI file already exists: $nifti_file_name"
    else
        bruker_to_nifti ${in_path} ${struct_scan_number} struct.nii.gz

    fi
    ## Cleaning the structural image by masking it with a manually created mask ##

    cleaned_struct_file_name="cleaned_struct.nii.gz"

    if is_step_enabled "masking_file"; then
        echo "Cleaning the structural image by manually creating mask"
        
        if [[ -f ${cleaned_struct_file_name} ]]; then
            echo "✔ Structural Image for Coregistration exists."
        else
            masking_file struct.nii.gz ${cleaned_struct_file_name} || exit 1
        fi
        structural_file_for_coregistration="${analysed_struct_dir}/cleaned_struct.nii.gz"    
    else
        echo "⏭️ Masks not generated, will not get cleaned data in further pipeline"    
        structural_file_for_coregistration="${analysed_struct_dir}/struct.nii.gz"    
    fi
else
    echo "⏭️ Data Conversion Step not executed, may not have data in further pipeline"
fi


# ==============================================================================================================================================================================================================
# Converting Functional Data into NIFTI and Processing Functional Data
# ==============================================================================================================================================================================================================

## Converting Functional Data into NIFTI ##
cd ${analysed_func_dir}

if is_step_enabled "bruker_to_nifti"; then
    echo "Converting Bruker to NIFTI: Functional Data"
    
    nifti_file_name="func.nii.gz"
    if [[ -f "$nifti_file_name" ]]; then
        echo "✔ NIfTI file already exists: $nifti_file_name"
    else
        bruker_to_nifti ${in_path} ${func_scan_number} func.nii.gz
    fi
else
    echo "⏭️ Data Conversion Step not executed, may not have data in further pipeline"
fi

input_name_for_next_step="func"
# ==============================================================================================================================================================================================================
# Generating Masks
# ==============================================================================================================================================================================================================

## Extracting functional parameters needed for to extract middle volume ##
n_vols=$(func_param_extract "$path_raw_func" false | grep '^NoOfRepetitions=' | cut -d= -f2)
tr=$(func_param_extract "$path_raw_func" false | grep '^VolTR=' | cut -d= -f2)
middle_vol=$(func_param_extract "$path_raw_func" false | grep '^MiddleVolume=' | cut -d= -f2)

## Extracted middle volume for generating masks and for motion correction ##
extract_middle_volume ${input_name_for_next_step} "middle_vol.nii.gz" ${middle_vol}

mask_name_wo_cannulas="mask_func.nii.gz"
mask_name_with_cannulas="mask_func_with_cannulas.nii.gz"

## Generating masks for cleaning the data ##
if is_step_enabled "masking_file"; then
    echo "Generating masks for functional data"

    if [[ -f ${mask_name_wo_cannulas} && -f ${mask_name_with_cannulas} ]]; then
        echo "✔ Mask file for functional data exists."
    else
        echo "Create and save masks using fsleyes now and save them as ${mask_name_wo_cannulas} and ${mask_name_with_cannulas}"
        fsleyes "middle_vol.nii.gz"
    fi

    fslmaths middle_vol.nii.gz -mas ${mask_name_with_cannulas} "cleaned_middle_vol.nii.gz"
    mean_image_for_coregistration="cleaned_middle_vol.nii.gz"
else
    echo "⏭️ Masks not generated, will not have cleaned data in further pipeline"
    mean_image_for_coregistration="middle_vol.nii.gz"    
fi


# ==============================================================================================================================================================================================================
# Applying Motion Correction on raw functional data and plotting motion parameters
# ==============================================================================================================================================================================================================

if is_step_enabled "motion_correction"; then
     echo "Applying Motion Correction on raw functional data and plotting motion parameters"

    if [[ -f "mc_${input_name_for_next_step}.nii.gz" ]]; then
        echo "✔ Motion Corrected functional data exists. Skipping motion correction."
    else
        motion_correction "middle_vol.nii.gz" "${input_name_for_next_step}.nii.gz" "mc_${input_name_for_next_step}"
    fi

    input_name_for_next_step="mc_${input_name_for_next_step}"
else
    echo "⏭️ Motion Correction not executed"
fi
# ==============================================================================================================================================================================================================
# Move into further subdirectory with in analysed data folder as well
# ==============================================================================================================================================================================================================

cd ${analysis_func_subdir}
cp ../${input_name_for_next_step}.nii.gz ./

#### Copying mask files into the working directory for later use ####
for f in ../mask_func.nii.gz ../mask_func_with_cannulas.nii.gz; do
  if [[ -f "$f" ]]; then
    cp "$f" ./
  fi
done

# ==============================================================================================================================================================================================================
# Applying temporal smoothing to the motion corrected functional data to see the temporal signatures
# ==============================================================================================================================================================================================================

## Applying temporal smoothing to the motion corrected functional data to see the temporal signatures ##
if is_step_enabled "smooth_movavg"; then
    echo "Applying temporal smoothing to the motion corrected functional data to see the temporal signatures"
    smooth_movavg "${input_name_for_next_step}.nii.gz" "ts_${input_name_for_next_step}.nii.gz" "${temp_smooth_factor}" "${tr}"

    if is_step_enabled "masking_file"; then
        echo "Cleaning the temporally smoothed motion-corrected functional data by applying mask"
        
        cleaned_ts_func_file_name="cleaned_ts_${input_name_for_next_step}.nii.gz"

        if [[ -f ${cleaned_ts_func_file_name} ]]; then
            echo "✔ Temporally Smoothed Motion Corrected Functional Image for further analysis exists."
        else
            masking_file "ts_${input_name_for_next_step}.nii.gz" ${cleaned_ts_func_file_name} || exit 1
        fi
    else
        echo "⏭️ Masks not generated, will not get cleaned data in further pipeline"    
    fi

    input_name_for_next_step="ts_${input_name_for_next_step}"    
else
    echo "⏭️ Temporal Smoothing not executed"
fi


## Based on the temporal smoothing, selecting baseline and signal volumes ##
echo "Choose your baseline and signal volumes from the temporally smoothed motion-corrected functional data."
fsleyes "${input_name_for_next_step}.nii.gz"
idx_vals=$(dialog \
  --clear \
  --title "Baseline / Signal Selection" \
  --form "Enter volume indices:" \
  10 60 2 \
  "Baseline Start Index:" 1 1 "10" 1 30 10 0 \
  "Signal Start Index:"   2 1 "10" 2 30 10 0 \
  3>&1 1>&2 2>&3
)

mapfile -t idx_vals <<< "$idx_vals"

base_start_idx="${idx_vals[0]}"
sig_start_idx="${idx_vals[1]}"

base_end_idx=$((base_start_idx + win_vol))
sig_end_idx=$((sig_start_idx + win_vol))

echo "Baseline Start Index : $base_start_idx to $base_end_idx"
echo "Signal Start Index   : $sig_start_idx to $sig_end_idx"

# ==============================================================================================================================================================================================================
# tSNR estimation
# ==============================================================================================================================================================================================================

if is_step_enabled "tSNR"; then
    echo "Computing Temporal SNR on the temporally smoothed motion-corrected functional data "
    tSNR "${input_name_for_next_step}.nii.gz"
else
    echo "⏭️ Temporal SNR not executed"

fi
# ==============================================================================================================================================================================================================
# Applying spatial smoothing
# ==============================================================================================================================================================================================================

## Applying spatial smoothing to the temporally smoothed motion-corrected functional data ##

if is_step_enabled "spatial_smoothing"; then
    echo "Applying spatial smoothing to the temporally smoothed motion-corrected functional data"
    
    dim_inplane=$(fslhd "${input_name_for_next_step}.nii.gz" | awk '/^pixdim1/ {print $2}')

    sigma_val=$(echo "scale=6; ($func_smooth_factor * $dim_inplane) / 2.3548" | bc)

    fslmaths "${input_name_for_next_step}.nii.gz" -s ${sigma_val} "sm_${input_name_for_next_step}.nii.gz"
    
    if is_step_enabled "masking_file"; then
        echo "Cleaning the spatially smoothed temporally smoothed motion-corrected functional data by applying mask"
        
        func_file_name_for_scm="cleaned_sm_${input_name_for_next_step}.nii.gz"

        if [[ -f ${func_file_name_for_scm} ]]; then
            echo "✔ Spatially Smoothed Temporally Smoothed Motion Corrected Functional Image for further analysis exists."
        else
            fslmaths "sm_${input_name_for_next_step}.nii.gz" -mas ${mask_name_with_cannulas} ${func_file_name_for_scm}
        fi
    else
        echo "⏭️ Masks not generated, will not get cleaned data in further pipeline"   
        func_file_name_for_scm="sm_${input_name_for_next_step}.nii.gz" 
    fi

    input_name_for_next_step="sm_${input_name_for_next_step}"
    
else
    echo "⏭️ Spatial Smoothing not executed, signal change maps may be noisy"
fi


# ----------------------------------------------------------
# 🔐 SAFETY FALLBACK 
# ----------------------------------------------------------
if [[ -z "${func_file_name_for_scm:-}" ]]; then
  func_file_name_for_scm="${input_name_for_next_step}.nii.gz"
fi



# ==============================================================================================================================================================================================================
# Generating Signal Change Maps and Signal Change Time Series
# ==============================================================================================================================================================================================================

## Generating Signal Change Maps and Signal Change Time Series ##

if is_step_enabled "signal_change_map"; then

    echo "Generating Signal Change Map and Signal Change Time Series"
    signal_change_map "${func_file_name_for_scm}" "$base_start_idx" "$base_end_idx" "$sig_start_idx" "$sig_end_idx" scm_mc_func


    ## Cleaning the Signal Change Map and Signal Change Time Series by applying mask ##
    if is_step_enabled "masking_file"; then
        echo "Cleaning the Signal Change Map and Signal Change Time Series by applying mask"
        
        cleaned_scm_file_name="cleaned_scm_mc_func.nii.gz"
        cleaned_scm_ts_file_name="cleaned_time_series_scm_mc_func.nii.gz"

        ### Checking if cleaned files already exist ###
        if [[ -f ${cleaned_scm_file_name} && -f ${cleaned_scm_ts_file_name} ]]; then
            echo "✔ Cleaned Signal Change Map and Signal Change Time Series for further analysis exists."
        else
            masking_file scm_mc_func.nii.gz ${cleaned_scm_file_name} || exit 1
            masking_file time_series_scm_mc_func.nii.gz ${cleaned_scm_ts_file_name} || exit 1
        fi
    else
        echo "⏭️ Masks not generated, will not get cleaned data in further pipeline"    
    fi
else
    echo "⏭️ Signal Change Map not estimated"
fi

# ==============================================================================================================================================================================================================
# Coregistering functional time series and functional signal change map to structural images using AFNI
# ==============================================================================================================================================================================================================



affine_matrix_file="mean_func_struct_aligned.aff12.1D"
reference_file=${structural_file_for_coregistration}
input_file1=${mean_image_for_coregistration}
apply_affine="true" #### will always be set to true as we want to apply the affine estimated from mean functional image to SCM and functional time series ####

#### Deciding input files #2 i.e cleaned or uncleaned SCM based on whether masking was applied or not ####
if is_step_enabled "masking_file"; then
    input_file2="cleaned_scm_mc_func.nii.gz"
else
    input_file2="scm_mc_func.nii.gz"
fi

uncleaned_coregisterd_scm_file="signal_change_map_coregistered_structural_space.nii.gz"
coregistered_time_series_file="fMRI_coregistered_to_struct.nii.gz"

#### Coregistration to be applied or not ####
if is_step_enabled "coregistration_afni"; then
    echo "Coregistering functional time series and functional signal change map to structural images using AFNI"
    
    #### Based on presence of affine matrix, decide if affine matrix needs to be estimated or not ####
    if [[ -f ${affine_matrix_file} ]]; then
        echo "✔ Affine matrix for SCM coregistration already exists. Reusing it."
        estimate_affine="false"
    else
        echo "Estimating affine matrix for SCM coregistration."
        estimate_affine="true"
    fi
    
    #### Coregistering functional SCM to structural images ####
    coregistration_afni ${input_file1} ${input_file2} ${reference_file} mean_func_coreg_to_struct.nii.gz scm_coreg_to_struct.nii.gz ${estimate_affine} ${apply_affine} ${affine_matrix_file}
    
    #### Coregistering functional time series to structural images ####
    coregistration_afni "" ${func_file_name_for_scm} ${reference_file} "" fMRI_coregistered_to_struct.nii.gz false true ${affine_matrix_file}
    
    mean_coregistered_image="mean_func_coreg_to_struct.nii.gz"
    fslmaths fMRI_coregistered_to_struct.nii.gz -Tmean ${mean_coregistered_image}

    #### Generating Signal Change Maps and Signal Change Time Series of the coregistered functional data ####
    if is_step_enabled "signal_change_map"; then

        echo "Generating Signal Change Map and Signal Change Time Series of the coregistered functional data"
        signal_change_map fMRI_coregistered_to_struct.nii.gz ${base_start_idx} ${base_end_idx} ${sig_start_idx} ${sig_end_idx} scm_coregistered_structural_space
    
        if is_step_enabled "masking_file"; then
            echo "Cleaning the Signal Change Map generated from the coregistered functional data by applying mask"
            
            cleaned_coreg_scm_file_name="cleaned_scm_coregistered_structural_space.nii.gz"
            cleaned_coreg_scm_ts_file_name="cleaned_time_series_scm_coregistered_structural_space.nii.gz"
    
            ### Checking if cleaned files already exist ###
            if [[ -f ${cleaned_coreg_scm_file_name} && -f ${cleaned_coreg_scm_ts_file_name} ]]; then
                echo "✔ Cleaned Signal Change Map and Signal Change Time Series of the coregistered functional data for further analysis exists."
            else
                masking_file scm_coregistered_structural_space.nii.gz ${cleaned_coreg_scm_file_name} || exit 1
                masking_file time_series_scm_coregistered_structural_space.nii.gz ${cleaned_coreg_scm_ts_file_name} || exit 1
            fi
        else
            echo "⏭️ Masks not generated, will not get cleaned data in further pipeline"
    fi

else
    echo "⏭️ Coregistration not executed"
fi

# ==============================================================================================================================================================================================================
# Marking ROIs and saving time courses
# ==============================================================================================================================================================================================================


if is_step_enabled "roi_analysis"; then
    echo "Marking ROIs and saving time courses"

    
    ### Detecting ROIs safely ###
    
    shopt -s nullglob
    roi_files=(roi_*.nii.gz)
    shopt -u nullglob

    if [[ ${#roi_files[@]} -eq 0 ]]; then
        echo "✔ ROI files not found."
        echo "Please create ROI files named as:"
        echo "roi_{Neurotransmitter}_{Direct or AAV}_{Concentration}_{Functional or Coregistered}_{Left or Right}.nii.gz"

        #### Display images for ROI creation ###
        if [[ -f "${mean_coregistered_image:-}" ]]; then
            fsleyes "${structural_file_for_coregistration}" "${mean_coregistered_image}" &
        else
            fsleyes "${structural_file_for_coregistration}" &
        fi

        fsleyes ../middle_vol.nii.gz &
        exit 0
    else
        echo "✔ ROI files found. Proceeding to data analysis."
    fi
    ### ROI Loop ###
    roi_left=""
    roi_right=""

    roi_names=()
    roi_hemis=()
    roi_means=()

    for roi in "${roi_files[@]}"; do

        roi_base=$(basename "$roi")
        #### Determine hemisphere ####
        if [[ "$roi_base" == *left.nii.gz ]]; then
            hemi="left"
            roi_left="$roi"
        elif [[ "$roi_base" == *right.nii.gz ]]; then
            hemi="right"
            roi_right="$roi"
        else
            echo "[WARN] Cannot determine hemisphere for $roi_base, skipping"
            continue
        fi

        #### Determine correct functional file ####
        if [[ "$roi_base" == *functional* ]]; then
            func_to_use="${func_file_name_for_scm}"
        elif [[ "$roi_base" == *coregistered* ]]; then
            func_to_use="fMRI_coregistered_to_struct.nii.gz"
        else
            echo "[WARN] Cannot determine space for $roi_base, skipping"
            continue
        fi

        echo
        echo "[INFO] Running ROI analysis"
        echo "       ROI  : $roi_base"
        echo "       Hemi : $hemi"
        echo "       Func : $func_to_use"

        mean_psc=$(roi_analysis "$func_to_use" "$roi" "$base_start_idx" "$base_end_idx" "$sig_start_idx" "$sig_end_idx" "$tr") || {echo "[ERROR] ROI analysis failed for $roi_base"; continue;}

        roi_names+=("$roi_base")
        roi_hemis+=("$hemi")
        roi_means+=("$mean_psc")

        echo "[RESULT] Mean PSC ($roi_base) = $mean_psc"

        ## live terminal plot for individual ROIs ##
        
        gnuplot << EOF
set terminal dumb size 120,30
set title "LIVE PSC – $roi_base"
set xlabel "Time (minutes)"
set ylabel "Signal Change (%)"
set grid
plot "PSC_time_series_${roi_base%.nii.gz}.txt" using (\$0*$tr/60.0):1 with lines title "$roi_base"
EOF

    done

    ### Combining left and right ROI plots ###    

    if [[ -n "$roi_left" && -n "$roi_right" ]]; then

        combined_plot="PSC_Time_Series_LEFT_vs_RIGHT.svg"

        gnuplot << EOF
set terminal svg size 900,400
set output "$combined_plot"
set title "Percent Signal Change – Left vs Right ROI"
set xlabel "Time (minutes)"
set ylabel "Signal Change (%)"

set style rect fc rgb "green" fs solid 0.2 noborder
set object rect from ($base_start_idx*$tr/60.0), graph 0 to ($base_end_idx*$tr/60.0), graph 1

set style rect fc rgb "blue" fs solid 0.2 noborder
set object rect from ($sig_start_idx*$tr/60.0), graph 0 to ($sig_end_idx*$tr/60.0), graph 1

plot \
  "PSC_time_series_${roi_left%.nii.gz}.txt"  using (\$0*$tr/60.0):1 with lines lw 2 lc rgb "red"  title "Left ROI", \
  "PSC_time_series_${roi_right%.nii.gz}.txt" using (\$0*$tr/60.0):1 with lines lw 2 lc rgb "blue" title "Right ROI"
EOF

        echo "[OK] Combined plot saved → $combined_plot"
    else
        echo "[INFO] Skipping combined plot (left or right ROI missing)"
    fi

else
    echo "⏭️ ROI Analysis not executed"
fi



# ==============================================================================================================================================================================================================
# Voxel-wise Correlation Analysis
# ==============================================================================================================================================================================================================

if [[ "$sequence_name_func" == "functionalEPI" ]]; then 

    curl -f -O https://raw.githubusercontent.com/njainmpi/amplify/main/seed_voxelwise_correlation.py || { echo "❌ Download failed"; exit 1; }
    curl -f -O https://raw.githubusercontent.com/njainmpi/amplify/main/compare_conditions.py || { echo "❌ Download failed"; exit 1; }
    chmod +x seed_voxelwise_correlation.py compare_conditions.py

    seed_used="roi_seed.nii.gz"

    ## Establishing correlation analysis before injection ##
    python3 seed_voxelwise_correlation.py --func "${input_name_for_next_step}.nii.gz" --seed-mask "${seed_used}" --target-mask "${mask_name_wo_cannulas}" --start 100 --stop 500 --output-prefix before

    ## Establishing correlation analysis during injection ##
    python3 seed_voxelwise_correlation.py --func "${input_name_for_next_step}.nii.gz" --seed-mask "${seed_used}" --target-mask "${mask_name_wo_cannulas}" --start 700 --stop 1100 --output-prefix during

    ## Establishing correlation analysis after injection ##
    python3 seed_voxelwise_correlation.py --func "${input_name_for_next_step}.nii.gz" --seed-mask "${seed_used}" --target-mask "${mask_name_wo_cannulas}" --start 1300 --stop 2300 --output-prefix after

    #### Comparing of three conditions ####
    python3 compare_conditions.py --before before_voxelwise_correlation.txt --during during_voxelwise_correlation.txt --after  after_voxelwise_correlation.txt --output-prefix injection_effect

    rm -f seed_voxelwise_correlation.py compare_conditions.py
else
    echo "⏭️ Voxel-wise correlation skipped (not functionalEPI)"
fi


# ==============================================================================================================================================================================================================
# Entering all data analysis record in the SQL Database
# ==============================================================================================================================================================================================================

DB="${root_location}/analysis_record.db"

#############################################
# CREATE DATABASE TABLES (SAFE TO RE-RUN)
#############################################

echo "[INFO] Initializing SQLite database"

sqlite3 "$DB" <<EOF
CREATE TABLE IF NOT EXISTS Analysis_Run (
    run_id INTEGER PRIMARY KEY AUTOINCREMENT,
    analysis_done_by TEXT,
    root_location TEXT,
    subject_id TEXT,
    analysis_folder_name TEXT,
    functional_run INTEGER,
    structural_run INTEGER,
    temporal_res REAL,
    window_duration INTEGER,
    spatial_smoothing REAL,
    spatial_smoothing_coreg_image REAL,
    start_idx_correlation INTEGER,
    end_idx_correlation INTEGER,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ROI_Results (
    roi_id INTEGER PRIMARY KEY AUTOINCREMENT,
    analysis_done_by TEXT,
    root_location TEXT,
    subject_id TEXT,
    analysis_folder_name TEXT,
    run_id INTEGER,
    roi_name TEXT,
    hemisphere TEXT,
    mean_signal_change REAL,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY(run_id) REFERENCES Analysis_Run(run_id)
);
EOF

#############################################
# INSERT ANALYSIS RUN METADATA
#############################################

echo "[INFO] Writing analysis metadata to database"

run_id=$(sqlite3 "$DB" <<EOF
INSERT INTO Analysis_Run (
  analysis_done_by,
  root_location,
  subject_id,
  analysis_folder_name,
  functional_run,
  structural_run,
  temporal_res,
  window_duration,
  spatial_smoothing,
  spatial_smoothing_coreg_image
)
VALUES (
  '${user}',
  '${root_location}',
  '${subject_id}',
  '${folder_created}',
  $func_scan_number,
  $struct_scan_number,
  $tr,
  $win_vol,
  $func_smooth_factor,
  $coreg_smooth_factor
);
SELECT last_insert_rowid();
EOF
)

echo "[OK] Analysis_Run inserted (run_id=$run_id)"

#############################################
# INSERT ROI RESULTS (ANY NUMBER OF ROIS)
#############################################

echo "[INFO] Writing ROI results to database"

for i in "${!roi_names[@]}"; do
    roi="${roi_names[$i]}"
    hemi="${roi_hemis[$i]}"
    mean="${roi_means[$i]}"

    if [[ -z "$mean" ]]; then
        mean_sql="NULL"
    else
        mean_sql="$mean"
    fi

    sqlite3 "$DB" <<EOF
INSERT INTO ROI_Results (
  analysis_done_by,
  root_location,
  subject_id,
  analysis_folder_name,
  run_id,
  roi_name,
  hemisphere,
  mean_signal_change
)
VALUES (
  '$user',
  '$root_location',
  '$subject_id',
  '$folder_created',
  $run_id,
  '$roi',
  '$hemi',
  $mean_sql
);
EOF
done

echo "[OK] ROI results written to database"

echo "📁 Database location: $DB"


