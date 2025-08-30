#!/usr/bin/env bash
set -Eeuo pipefail
set -o errtrace
PS4='+ ${BASH_SOURCE}:${LINENO}:${FUNCNAME[0]}: '
trap 'code=$?; echo "ERROR: command \"${BASH_COMMAND}\" exited $code at ${BASH_SOURCE[0]}:${LINENO}"; exit $code' ERR

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH." >&2; exit 127; }; }
need curl; need awk; need python3; command -v column >/dev/null 2>&1 || true

# --- Robust GitHub sourcing helper ---
gh_source() {
  local fname="$1"
  local repo_base="https://raw.githubusercontent.com/njainmpi/amplify/main"
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
  echo "ERROR: Could not fetch '$fname' from any known path. Tried: ${try_paths[*]}" >&2
  exit 2
}

# --- Source all helpers (your list) ---
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

identity="$(whoami)@$(hostname)"
datafile="individual_project_based_scripts/path_definition.txt"

matched_path=$(
  awk -v id="$identity" -F',' '
    $2 == id {
      path = $3
      for (i=4; i<=NF; i++) path = path "," $i
      print path
    }' "$datafile" 2>/dev/null || true
)

[[ -n "${matched_path:-}" && -d "$matched_path" ]] && root_location="$matched_path" || root_location="/Volumes/Extreme_Pro/fMRI"
csv="Animal_Experiments_Sequences.csv"
header_lines=2
DELIM=$'\x1f'

cd "$root_location/RawData"
[[ -f "$csv" ]] || { echo "ERROR: $csv not found in $(pwd)"; exit 1; }

echo -e "\n== Available rows (CSV line | first 8 columns) =="
python3 - "$csv" "$header_lines" <<'PY' | { column -t -s'|' || cat; }
import csv, sys
path = sys.argv[1]; skip_n = int(sys.argv[2])
def try_decode(b):
    for enc in ('utf-8-sig','cp1252','latin-1'):
        try: return b.decode(enc)
        except UnicodeDecodeError: continue
    return b.decode('utf-8', errors='replace')
def clean(s): return (s or '').replace('|','¦').strip()
print("Line|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8")
with open(path,'rb') as f:
    for lineno, raw in enumerate(f, start=1):
        if lineno <= skip_n: continue
        if not raw.strip(): continue
        row = next(csv.reader([try_decode(raw)]))
        row = [clean(c) for c in row]
        first8 = (row + [""]*8)[:8]
        print("|".join([str(lineno)] + first8))
PY

echo
read -rp "Enter the CSV LINE NUMBER to run (e.g., 5, 7, 10), or q to quit: " sel
case "$sel" in q|Q) echo "Aborted."; exit 0 ;; esac
[[ "$sel" =~ ^[0-9]+$ ]] || { echo "Invalid selection."; exit 1; }
(( sel > header_lines )) || { echo "Line $sel is within header (1..$header_lines)."; exit 1; }

parsed="$(
python3 - "$csv" "$sel" <<'PY'
import csv, sys
US = '\x1f'; path = sys.argv[1]; target = int(sys.argv[2])
def try_decode(b):
    for enc in ('utf-8-sig','cp1252','latin-1'):
        try: return b.decode(enc)
        except UnicodeDecodeError: continue
    return b.decode('utf-8', errors='replace')
with open(path,'rb') as f:
    for lineno, raw in enumerate(f, start=1):
        if lineno == target:
            txt = try_decode(raw)
            if not raw.strip() or txt.strip(', \t\r\n') == '':
                sys.stdout.write(''); sys.exit(0)
            row = next(csv.reader([txt]))
            row = [(c or '').replace(US, ' ') for c in row]
            sys.stdout.write(US.join(row)); break
PY
)"
[[ -n "$parsed" ]] || { echo "Selected line is blank."; exit 1; }

get_n_field(){ echo "$parsed" | cut -d"$DELIM" -f"$1"; }
trim(){ printf '%s' "$1" | xargs; }

project_name=$(trim "$(get_n_field 3)")
sub_project_name=$(trim "$(get_n_field 4)")
dataset_name=$(trim "$(get_n_field 2)")
structural_name=$(trim "$(get_n_field 5)")
functional_name=$(trim "$(get_n_field 6)")
struc_coregistration=$(trim "$(get_n_field 7)")
baseline_duration=$(trim "$(get_n_field 13)")
injection_duration=$(trim "$(get_n_field 14)")

export Project_Name="$project_name"
export Sub_project_Name="$sub_project_name"
export Dataset_Name="$dataset_name"
export structural_run="$structural_name"
export run_number="$functional_name"
export str_for_coreg="$struc_coregistration"
export baseline_duration_in_min="$baseline_duration"
export injection_duration_in_min="$injection_duration"

echo -e "\n== Selection summary =="
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

Path_Raw_Data="$root_location/RawData/$project_name/$sub_project_name"
Path_Analysed_Data="$root_location/AnalysedData/$project_name/$sub_project_name/$Dataset_Name"

datapath="$(find "$Path_Raw_Data" -type d -name "*${Dataset_Name}*" 2>/dev/null | head -n1 || true)"
[[ -n "${datapath:-}" ]] || { echo "ERROR: raw data dir for '$Dataset_Name' not found under $Path_Raw_Data"; exit 1; }
echo "Raw dataset path: $datapath"

echo
echo "Analysed Data folder check..."
if [[ -d "$Path_Analysed_Data" ]]; then
  echo "Analysed Data folder exists, Proceeding to Analyse the data"
else
  echo "Creating: $Path_Analysed_Data"
  mkdir -p "$Path_Analysed_Data"
fi
cd "$Path_Analysed_Data"

# ---------------- STRUCTURAL ----------------
# NOTE: If your helper defines FUNC_PARAM_EXTARCT (typo), change next line back.
FUNC_PARAM_EXTRACT "$datapath/$structural_run"

# Assert the helper exported these:
: "${SequenceName:?FUNC_PARAM_EXTRACT did not set SequenceName}"
# MiddleVolume might be set by motion function helper; assert later too if needed.
echo "DEBUG: SequenceName='$SequenceName' (structural)"

struct_dir="$Path_Analysed_Data/${structural_run}${SequenceName}"
mkdir -p "$struct_dir"
cd "$struct_dir"

# If your CHECK_FILE_EXISTENCE was creating dirs AND returns non-zero when missing, it would trigger ERR.
# We already mkdir -p + cd, so this is optional (won’t fail the script):
CHECK_FILE_EXISTENCE "$struct_dir" || true

run_if_missing "anatomy.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$structural_run" "$datapath/$structural_run/method"
cp -f G1_cp.nii.gz anatomy.nii.gz || echo "WARNING: G1_cp.nii.gz not found for structural; continuing."

# ---------------- FUNCTIONAL ----------------
FUNC_PARAM_EXTRACT "$datapath/$run_number"
: "${SequenceName:?FUNC_PARAM_EXTRACT did not set SequenceName}"
echo "DEBUG: SequenceName='$SequenceName' (functional)"

func_dir="$Path_Analysed_Data/${run_number}${SequenceName}"
mkdir -p "$func_dir"
cd "$func_dir"
CHECK_FILE_EXISTENCE "$func_dir" || true

run_if_missing "G1_cp.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$run_number" "$datapath/$run_number/method"

# ---------------- MOTION ----------------
: "${MiddleVolume:?FUNC_PARAM_EXTRACT (or motion helper) did not set MiddleVolume}"
PRINT_YELLOW "Performing Step 1: Motion Correction"
run_if_missing "mc_func.nii.gz" "mc_func+orig.HEAD" "mc_func+orig.BRIK" -- \
  MOTION_CORRECTION "$MiddleVolume" G1_cp.nii.gz mc_func

# ---------------- tSNR (AFNI) -------------
PRINT_YELLOW "Performing Step 2: Obtaining Mean func, Std func and tSNR Maps"
run_if_missing "tSNR_mc_func.nii.gz" "tSNR_mc_func+orig.HEAD" "tSNR_mc_func+orig.BRIK" -- \
  TEMPORAL_SNR_using_AFNI mc_func+orig

# ---------------- N4 Bias ------------------
PRINT_YELLOW "Performing Step 3: Performing N4 Bias Field Correction of mean_mc_func"
run_if_missing "cleaned_mc_func.nii.gz" -- \
  BIAS_CORRECTED_IMAGE mean_mc_func.nii.gz 100 mc_func.nii.gz

# ---------------- Means --------------------
if [[ -f mean_image_min_1_to_10.nii.gz ]]; then
  echo "Mean files exist. Skipping 3dTstat blocks."
else
  echo "Computing mean images for 0–40 min windows (10-min bins)."
  3dTstat -mean -prefix mean_image_min_1_to_10.nii.gz  cleaned_mc_func.nii.gz"[0..599]"
  3dTstat -mean -prefix mean_image_min_11_to_20.nii.gz cleaned_mc_func.nii.gz"[600..1199]"
  3dTstat -mean -prefix mean_image_min_21_to_30.nii.gz cleaned_mc_func.nii.gz"[1200..1799]"
  3dTstat -mean -prefix mean_image_min_31_to_40.nii.gz cleaned_mc_func.nii.gz"[1800..2399]"
fi

echo "Pipeline completed up to mean-image generation."
