#!/usr/bin/env bash
# ============================================================================
# data_processing.sh
# Author: Naman Jain
# Runs AMPLIFY data-processing by sourcing helper scripts directly from GitHub.
# Compatible with macOS, Linux, and WSL.
# ============================================================================

set -Eeuo pipefail

# --------------------------- Minimal dependency checks -----------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH." >&2; exit 127; }; }
need curl
need awk
need python3
need column || true   # nice-to-have for preview formatting

trap 'echo "ERROR at line $LINENO. Exiting."; exit 1' ERR

# ------------------------ Robust GitHub sourcing helper ----------------------
# Tries multiple candidate paths within the repo.
# Verifies HTTP 200 AND non-empty payload before sourcing.
gh_source() {
    local fname="$1"
    local repo_base="https://raw.githubusercontent.com/njainmpi/amplify/main"

    # Add/remove candidates to match your repo layout as needed
    local try_paths=(
        "$fname"
        "individual_project_based_scripts/$fname"
        "toolbox/$fname"
        "scripts/$fname"
    )

    local url status
    for p in "${try_paths[@]}"; do
        url="$repo_base/$p"
        status="$(curl -sS -L -w '%{http_code}' -o /tmp/ghsrc.$$ "$url" || echo 000)"
        if [[ "$status" == "200" ]] && [[ -s /tmp/ghsrc.$$ ]]; then
            # shellcheck source=/dev/null
            source /tmp/ghsrc.$$
            rm -f /tmp/ghsrc.$$
            echo "Sourced: $p"
            return 0
        fi
    done

    echo "ERROR: Could not fetch '$fname' from any known path." >&2
    echo "       Tried: ${try_paths[*]}" >&2
    exit 2
}

# ----------------------------- Source all helpers ----------------------------
gh_source toolbox_name.sh
gh_source log_execution.sh
gh_source missing_run.sh
gh_source folder_existence_function.sh
gh_source func_parameters_extraction.sh
gh_source bias_field_correction.sh
gh_source check_spikes.sh
gh_source coregistration.sh
gh_source data_conversion.sh
gh_source motion_correction.sh
gh_source quality_check.sh
gh_source signal_change_map.sh
gh_source smoothing_using_fsl.sh
gh_source temporal_snr_using_afni.sh
gh_source temporal_snr_using_fsl.sh
gh_source scm_visual.sh
gh_source print_function.sh

# ------------------------------- Identity/path -------------------------------
identity="$(whoami)@$(hostname)"

# Path map file lives in the repo (access via working copy or your run dir)
datafile="individual_project_based_scripts/path_definition.txt"

if [[ ! -f "$datafile" ]]; then
    echo "WARNING: $datafile not found relative to current dir: $(pwd)"
    echo "         If you run via process substitution, ensure CWD contains this file."
fi

# Extract the matching path (3rd..Nth columns) for current user@host
matched_path=$(
    awk -v id="$identity" -F',' '
        $2 == id {
            path = $3
            for (i=4; i<=NF; i++) path = path "," $i
            print path
        }' "$datafile" 2>/dev/null || true
)

echo "Running data analysis for $(whoami) on $(hostname)."
if [[ -n "${matched_path:-}" ]]; then
    echo "Matched root path from map: $matched_path"
fi

# ------------------------------- CONFIG --------------------------------------
# Choose root_location: prefer matched_path if it exists; else fallback (edit!)
if [[ -n "${matched_path:-}" && -d "$matched_path" ]]; then
    root_location="$matched_path"
else
    # Fallback used in your example
    root_location="/Volumes/Extreme_Pro/fMRI"
fi

csv="Animal_Experiments_Sequences.csv"
header_lines=2
DELIM=$'\x1f'   # ASCII Unit Separator

# ------------------------------- SETUP ---------------------------------------
cd "$root_location/RawData" || { echo "ERROR: cannot cd to $root_location/RawData"; exit 1; }
[[ -f "$csv" ]] || { echo "ERROR: $csv not found in $(pwd)"; exit 1; }

# ------------------------------- PREVIEW -------------------------------------
echo
echo "== Available rows (CSV line | first 8 columns) =="

python3 - "$csv" "$header_lines" <<'PY' | { column -t -s'|' || cat; }
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

print("Line|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8")

with open(path, 'rb') as f:
    for lineno, raw in enumerate(f, start=1):
        if lineno <= skip_n:
            continue
        if not raw.strip():
            continue
        row = next(csv.reader([try_decode(raw)]))
        row = [clean(c) for c in row]
        first8 = (row + [""]*8)[:8]
        print("|".join([str(lineno)] + first8))
PY

# ------------------------------ SELECTION ------------------------------------
echo
read -rp "Enter the CSV LINE NUMBER to run (e.g., 5, 7, 10), or q to quit: " sel
case "$sel" in q|Q) echo "Aborted."; exit 0 ;; esac
if ! [[ "$sel" =~ ^[0-9]+$ ]]; then
  echo "Invalid selection."; exit 1
fi
if (( sel <= header_lines )); then
  echo "Line $sel is within the header region (1..$header_lines). Pick a data line."; exit 1
fi

# ------------------------- PARSE SELECTED LINE -------------------------------
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

# -------------------------- Helpers for parsed fields ------------------------
get_n_field() { echo "$parsed" | cut -d"$DELIM" -f"$1"; }
trim() { printf '%s' "$1" | xargs; }

# Map CSV schema
project_name=$(trim "$(get_n_field 3)")
sub_project_name=$(trim "$(get_n_field 4)")
dataset_name=$(trim "$(get_n_field 2)")
structural_name=$(trim "$(get_n_field 5)")
functional_name=$(trim "$(get_n_field 6)")
struc_coregistration=$(trim "$(get_n_field 7)")
baseline_duration=$(trim "$(get_n_field 13)")
injection_duration=$(trim "$(get_n_field 14)")

# Export for downstream
export Project_Name="$project_name"
export Sub_project_Name="$sub_project_name"
export Dataset_Name="$dataset_name"
export structural_run="$structural_name"
export run_number="$functional_name"
export str_for_coreg="$struc_coregistration"
export baseline_duration_in_min="$baseline_duration"
export injection_duration_in_min="$injection_duration"

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

# ------------------------------ Path setup -----------------------------------
Path_Raw_Data="$root_location/RawData/$project_name/$sub_project_name"
Path_Analysed_Data="$root_location/AnalysedData/$project_name/$sub_project_name/$Dataset_Name"

# Find raw dataset dir; take first match
datapath="$(find "$Path_Raw_Data" -type d -name "*${Dataset_Name}*" 2>/dev/null | head -n1 || true)"
if [[ -z "${datapath:-}" ]]; then
    echo "ERROR: Could not find raw data dir for '$Dataset_Name' under: $Path_Raw_Data" >&2
    exit 1
fi
echo "Raw dataset path: $datapath"

echo
echo "Dataset Currently Being Analysed: $Dataset_Name"
echo "Project: $Project_Name | Subproject: $Sub_project_Name"
echo "Structural run: $structural_run | Functional run: $run_number"
echo

mkdir -p "$Path_Analysed_Data"

cd "$Path_Analysed_Data"

# ----------------------- Structural: convert & prepare -----------------------
# NOTE: the function name was misspelled previously as FUNC_PARAM_EXTARCT.
# If your helper defines a different name, adjust here accordingly.
FUNC_PARAM_EXTARCT "$datapath/$structural_run"

CHECK_FILE_EXISTENCE "$Path_Analysed_Data/${structural_run}${SequenceName}"
cd "$Path_Analysed_Data/${structural_run}${SequenceName}"

run_if_missing "anatomy.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$structural_run" "$datapath/$structural_run/method"
cp -f G1_cp.nii.gz anatomy.nii.gz || echo "WARNING: G1_cp.nii.gz not found for structural; continuing."

# ----------------------- Functional: convert & prepare -----------------------
FUNC_PARAM_EXTARCT "$datapath/$run_number"

CHECK_FILE_EXISTENCE "$Path_Analysed_Data/${run_number}${SequenceName}"
cd "$Path_Analysed_Data/${run_number}${SequenceName}"

run_if_missing "G1_cp.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$run_number" "$datapath/$run_number/method"

# ----------------------- Step 1: Motion Correction ---------------------------
PRINT_YELLOW "Performing Step 1: Motion Correction"
run_if_missing "mc_func.nii.gz" "mc_func+orig.HEAD" "mc_func+orig.BRIK" -- \
    MOTION_CORRECTION "$MiddleVolume" G1_cp.nii.gz mc_func

# ----------------------- Step 2: tSNR (AFNI) --------------------------------
PRINT_YELLOW "Performing Step 2: Obtaining Mean func, Std func and tSNR Maps"
run_if_missing "tSNR_mc_func.nii.gz" "tSNR_mc_func+orig.HEAD" "tSNR_mc_func+orig.BRIK" -- \
    TEMPORAL_SNR_using_AFNI mc_func+orig

# ----------------------- Step 3: N4 Bias Field -------------------------------
PRINT_YELLOW "Performing Step 3: Performing N4 Bias Field Correction of mean_mc_func"
run_if_missing "cleaned_mc_func.nii.gz" -- \
    BIAS_CORRECTED_IMAGE mean_mc_func.nii.gz 100 mc_func.nii.gz
# -b (input #2) [54,3] example: start with 32-pt scale ... adjust in function if needed.

# ----------------------- Epoch-wise means (AFNI) -----------------------------
if [[ -f mean_image_min_1_to_10.nii.gz ]]; then
   echo "Mean files exist. Skipping 3dTstat blocks."
else
   echo "Computing mean images for 0–40 min windows (10-min bins)."
   # Adjust frame ranges to your TR and total timepoints
   3dTstat -mean -prefix mean_image_min_1_to_10.nii.gz  cleaned_mc_func.nii.gz"[0..599]"
   3dTstat -mean -prefix mean_image_min_11_to_20.nii.gz cleaned_mc_func.nii.gz"[600..1199]"
   3dTstat -mean -prefix mean_image_min_21_to_30.nii.gz cleaned_mc_func.nii.gz"[1200..1799]"
   3dTstat -mean -prefix mean_image_min_31_to_40.nii.gz cleaned_mc_func.nii.gz"[1800..2399]"
fi

echo "Pipeline completed up to mean-image generation."
