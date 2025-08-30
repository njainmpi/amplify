#!/bin/bash

#Following script has been made by Naman Jain with following features included in the
#different version upgrades

##Calling all the functions that will be used in the upcoming script

#01.04.2025: Intial Script Planned, all functions called through external script


bash ../toolbox_name.sh
source ../log_execution.sh
source ../missing_run.sh
source ../folder_existence_function.sh
source ../func_parameters_extraction.sh
source ../bias_field_correction.sh
source ../check_spikes.sh
source ../coregistration.sh
source ../data_conversion.sh
source ../motion_correction.sh
source ../quality_check.sh
source ../signal_change_map.sh
source ../smoothing_using_fsl.sh
source ../temporal_snr_using_afni.sh
source ../temporal_snr_using_fsl.sh
source ../scm_visual.sh
source ../print_function.sh

currentpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd ..
path_for_python_script_time_course="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ts_roi_python_script=$path_for_python_script_time_course/time_course_single_subject.py

##In order to use awk, you need to convert xlsx file to csv file


identity="$(whoami)@$(hostname)"


# Path to your data file
datafile="individual_project_based_scripts/path_definition.txt"  # replace with actual path

# Extract the matching path (3rd column)

matched_path=$(awk -v id="$identity" -F',' '
    $2 == id {
        path = $3
        for (i=4; i<=NF; i++) {
            path = path "," $i
        }
        print path
    }
' "$datafile")

# echo "$matched_path"

echo "Running data anaylysis for $(whoami) on system $(hostname) with $matched_path as root location."
############################################
# CONFIG
############################################
# root_location="$matched_path"                 # server path (uncomment when needed)
root_location="/Volumes/Extreme_Pro/fMRI"       # local troubleshooting path
csv="Animal_Experiments_Sequences.csv"

# Your CSV shows 2 header/notes lines to hide from preview/selection
header_lines=2

# Safe delimiter for passing parsed fields (bash cannot store NULs)
DELIM=$'\x1f'   # ASCII Unit Separator

############################################
# SETUP
############################################
cd "$root_location/RawData"
[[ -f "$csv" ]] || { echo "ERROR: $csv not found in $(pwd)"; exit 1; }

############################################
# PREVIEW (first 8 columns; aligned; real CSV line numbers)
############################################
echo
echo "== Available rows (CSV line | first 8 columns) =="

python3 - "$csv" "$header_lines" <<'PY' | column -t -s'|'
import csv, sys

path = sys.argv[1]
skip_n = int(sys.argv[2])

def try_decode(b: bytes) -> str:
    for enc in ('utf-8-sig', 'cp1252', 'latin-1'):
        try:
            return b.decode(enc)
        except UnicodeDecodeError:
            continue
    return b.decode('utf-8', errors='replace')

def clean(cell: str) -> str:
    return (cell or '').replace('|','¦').strip()

# Header for preview: line number + first 8 CSV columns
print("Line|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8")

with open(path, 'rb') as f:
    for lineno, raw in enumerate(f, start=1):
        if lineno <= skip_n:
            continue
        if not raw.strip():              # skip blank lines in display
            continue
        row = next(csv.reader([try_decode(raw)]))
        row = [clean(c) for c in row]
        first8 = (row + [""]*8)[:8]
        print("|".join([str(lineno)] + first8))
PY

############################################
# SELECTION (by physical CSV line number)
############################################
echo
read -rp "Enter the CSV LINE NUMBER to run (e.g., 5, 7, 10), or q to quit: " sel
case "$sel" in q|Q) echo "Aborted."; exit 0 ;; esac
if ! [[ "$sel" =~ ^[0-9]+$ ]]; then
  echo "Invalid selection."; exit 1
fi
if (( sel <= header_lines )); then
  echo "Line $sel is within the header region (1..$header_lines). Pick a data line."; exit 1
fi

############################################
# PARSE SELECTED LINE DIRECTLY IN PYTHON
############################################
parsed="$(
python3 - "$csv" "$sel" <<'PY'
import csv, sys

US = '\x1f'
path = sys.argv[1]
target = int(sys.argv[2])

def try_decode(b: bytes) -> str:
    for enc in ('utf-8-sig', 'cp1252', 'latin-1'):
        try:
            return b.decode(enc)
        except UnicodeDecodeError:
            continue
    return b.decode('utf-8', errors='replace')

with open(path, 'rb') as f:
    for lineno, raw in enumerate(f, start=1):
        if lineno == target:
            # reject blank/comma-only rows right here
            txt = try_decode(raw)
            if not raw.strip() or txt.strip(', \t\r\n') == '':
                sys.stdout.write('')
                sys.exit(0)
            row = next(csv.reader([txt]))
            row = [(c or '').replace(US, ' ') for c in row]
            sys.stdout.write(US.join(row))
            break
PY
)"

if [[ -z "$parsed" ]]; then
  echo "Line $sel is blank or contains only commas. Pick a non-blank data line."
  exit 1
fi

############################################
# Helper: get Nth field from $parsed (1-based)
############################################
get_n_field() {
  local n="$1"
  # 'cut' handles non-printable delimiters robustly
  echo "$parsed" | cut -d"$DELIM" -f"$n"
}

# Trim helper
trim() { printf '%s' "$1" | xargs; }

############################################
# MAP FIELDS to your CSV schema:
# 3 → Project, 4 → SubProject, 2 → Animal Id (Dataset),
# 5 → Structural, 6 → Functional, 7 → Coreg,
# 13 → Baseline, 14 → Injection
############################################
project_name=$(trim "$(get_n_field 3)")
sub_project_name=$(trim "$(get_n_field 4)")
dataset_name=$(trim "$(get_n_field 2)")
structural_name=$(trim "$(get_n_field 5)")
functional_name=$(trim "$(get_n_field 6)")
struc_coregistration=$(trim "$(get_n_field 7)")
baseline_duration=$(trim "$(get_n_field 13)")
injection_duration=$(trim "$(get_n_field 14)")

############################################
# EXPORT ENV VARS (for downstream pipeline)
############################################
export Project_Name="$project_name"
export Sub_project_Name="$sub_project_name"
export Dataset_Name="$dataset_name"
export structural_run="$structural_name"
export run_number="$functional_name"
export str_for_coreg="$struc_coregistration"
export baseline_duration_in_min="$baseline_duration"
export injection_duration_in_min="$injection_duration"

############################################
# SUMMARY (only — no raw field dump)
############################################
echo
echo "== Selection summary =="
echo "CSV line:                  $sel"
echo "Project_Name (col3):       $Project_Name"
echo "Sub_project_Name (col4):   $Sub_project_Name"
echo "Dataset_Name (col2):       $Dataset_Name"
echo "structural_run (col5):     $structural_run"
echo "run_number (col6):         $run_number"
echo "str_for_coreg (col7):      $str_for_coreg"
echo "baseline (min) (col13):    $baseline_duration_in_min"
echo "injection (min) (col14):   $injection_duration_in_min"
echo


# echo $Structural_Data

Path_Raw_Data="$root_location/RawData/$project_name/$sub_project_name"
Path_Analysed_Data="$root_location/AnalysedData/$project_name/$sub_project_name/$Dataset_Name"

datapath=$(find "$Path_Raw_Data" -type d -name "*${Dataset_Name}*" 2>/dev/null)
echo "$datapath"     

echo ""
echo ""
echo "Dataset Currently Being Analysed is $Dataset_Name" 
echo "from $Project_Name with Subproject $Sub_project_Name"
echo "for Structural run number $structural_run and Functional run number $run_number" 
echo ""
echo ""

if [ -d "$Path_Analysed_Data" ]; then
    echo "$Path_Analysed_Data does exist."
else
    mkdir $Path_Analysed_Data
fi

cd $Path_Analysed_Data

echo ""
echo ""

#conversion for structural data
FUNC_PARAM_EXTARCT $datapath/$structural_run
               
CHECK_FILE_EXISTENCE "$Path_Analysed_Data/$structural_run""$SequenceName"
cd $Path_Analysed_Data/$structural_run''$SequenceName

run_if_missing "anatomy.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$structural_run" "$datapath/$structural_run/method"
cp G1_cp.nii.gz anatomy.nii.gz
    
#conversion for functional data
FUNC_PARAM_EXTARCT $datapath/$run_number
CHECK_FILE_EXISTENCE "$Path_Analysed_Data/$run_number$SequenceName"
cd $Path_Analysed_Data/$run_number''$SequenceName
    
run_if_missing "G1_cp.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$run_number" "$datapath/$run_number/method"
        
## Function to perform motion correction
PRINT_YELLOW "Performing Step 1: Motion Correction"
run_if_missing "mc_func.nii.gz" "mc_func+orig.HEAD" "mc_func+orig.BRIK" -- MOTION_CORRECTION "$MiddleVolume" G1_cp.nii.gz mc_func

       
## Function to obtain tSNR Maps        
PRINT_YELLOW "Performing Step 2: Obtaining Mean func, Std func and tSNR Maps"
run_if_missing  "tSNR_mc_func.nii.gz" "tSNR_mc_func+orig.HEAD" "tSNR_mc_func+orig.BRIK" -- TEMPORAL_SNR_using_AFNI mc_func+orig

        
# Function to perform Bias Field Corrections
PRINT_YELLOW "Performing Step 3: Performing N4 Bias Field Correction of mean_mc_func"
run_if_missing  "cleaned_mc_func.nii.gz" -- BIAS_CORRECTED_IMAGE mean_mc_func.nii.gz 100 mc_func.nii.gz
# -b (Inpout #2 in above command) [54,3] means start with 32 points scale (equiv 20mm coil divided by 0.375mm resolution) with 3rd order b-spline


# Obtaining mean images across three regimes: Pre Injection, During Injection, Post Injection

if [ -f mean_image_min_1_to_10.nii.gz ]; then
   echo "Mean files exists."
else
   echo "Mean file does not exist."
   3dTstat -mean -prefix mean_image_min_1_to_10.nii.gz cleaned_mc_func.nii.gz"[0..599]"
   3dTstat -mean -prefix mean_image_min_11_to_20.nii.gz cleaned_mc_func.nii.gz"[600..1199]"
   3dTstat -mean -prefix mean_image_min_21_to_30.nii.gz cleaned_mc_func.nii.gz"[1199..1799]"
   3dTstat -mean -prefix mean_image_min_31_to_40.nii.gz cleaned_mc_func.nii.gz"[1800..2399]"
fi