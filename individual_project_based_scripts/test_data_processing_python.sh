#!/usr/bin/env bash
set -Eeuo pipefail
set -o errtrace
PS4='+ ${BASH_SOURCE}:${LINENO}:${FUNCNAME[0]}: '
trap 'code=$?; echo "ERROR: command \"${BASH_COMMAND}\" exited $code at ${BASH_SOURCE[0]}:${LINENO}"; exit $code' ERR

# --- ensure we're using bash, not sh ---
if [[ -z "${BASH_VERSION:-}" ]]; then
  echo "Please run this script with bash (not sh)." >&2
  exit 1
fi


# ===== Amplify: ensure Python deps (Py 3.7–3.12) =====
echo "[deps] Checking Python environment..."

if command -v python3 >/dev/null 2>&1; then
  PYTHON=python3
elif command -v python >/dev/null 2>&1; then
  PYTHON=python
else
  echo "ERROR: python3/python not found in PATH." >&2
  exit 1
fi

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  VENV_DIR=".venv"
  if [[ ! -d "$VENV_DIR" ]]; then
    echo "[deps] Creating venv at $VENV_DIR"
    "$PYTHON" -m venv "$VENV_DIR"
  fi
  # shellcheck disable=SC1090
  source "$VENV_DIR/bin/activate"
  PYTHON=python
else
  echo "[deps] Conda detected ($CONDA_PREFIX); using current env."
fi

$PYTHON -m pip install --upgrade pip setuptools wheel >/dev/null
$PYTHON -m pip uninstall -y argparse >/dev/null 2>&1 || true

REQ_URL="https://raw.githubusercontent.com/njainmpi/amplify/main/requirements.txt"
REQ_FILE="$(mktemp)"
curl -fsSL "$REQ_URL" -o "$REQ_FILE" || { echo "ERROR: unable to download requirements.txt"; exit 1; }

echo "[deps] Installing Python dependencies..."
$PYTHON -m pip install -r "$REQ_FILE"

if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "[deps] FFmpeg not found. Install manually (apt-get or brew)."
fi
echo "[deps] Python deps ready."
# ===== /Amplify deps =====



need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH." >&2; exit 127; }; }
need curl; need awk; need python3; command -v column >/dev/null 2>&1 || true
need mktemp

# Nice glob behaviour for optional files like Static_SCM*.nii.gz
shopt -s nullglob

# --- Expand user paths like ~, $HOME, relative → absolute (safe; no eval) ---
expand_path() {
  python3 - "$1" <<'PY_EXPAND'
import os, sys
p = sys.argv[1] if len(sys.argv)>1 else ''
print(os.path.abspath(os.path.expanduser(os.path.expandvars(p or '.'))))
PY_EXPAND
}

# --- Portable readline prompt with default (works on macOS Bash 3.2) ---
# Usage: read_default "Prompt text" "DEFAULT" varname
read_default() {
  local _prompt="$1" _default="$2" _outvar="$3"
  if help read 2>/dev/null | grep -q ' -i '; then
    read -e -r -p "${_prompt} [${_default}]: " -i "${_default}" REPLY || true
  else
    read -e -r -p "${_prompt} [${_default}]: " REPLY || true
  fi
  printf -v "${_outvar}" '%s' "${REPLY:-${_default}}"
}

# --- Robust GitHub sourcing helper (Bash files only) ---
gh_source() {
  local fname="$1"
  local repo_base="https://raw.githubusercontent.com/njainmpi/amplify/main"
  local try_paths=(
    "$fname"
    "individual_project_based_scripts/$fname"
    "toolbox/$fname"
    "scripts/$fname"
  )
  local url status tmp="/tmp/ghsrc.$$"
  for p in "${try_paths[@]}"; do
    url="$repo_base/$p"
    status="$(curl -sS -L -w '%{http_code}' -o "$tmp" "$url" || echo 000)"
    if [[ "$status" == "200" ]] && [[ -s "$tmp" ]]; then
      # shellcheck source=/dev/null
      source "$tmp"
      rm -f "$tmp"
      echo "Sourced: $p"
      return 0
    fi
  done
  echo "ERROR: Could not fetch '$fname' from any known path. Tried: ${try_paths[*]}" >&2
  exit 2
}

# --- Run a Python file from local or GitHub ---
# Usage: gh_py_exec make_static_maps_and_movie.py --args...
gh_py_exec() {
  local fname="$1"; shift || true
  local local_py="./$fname"
  if [[ -f "$local_py" ]]; then
    echo "Running local Python: $local_py $*"
    python3 "$local_py" "$@"
    return $?
  fi

  local repo_base="https://raw.githubusercontent.com/njainmpi/amplify/main"
  local try_paths=(
    "$fname"
    "individual_project_based_scripts/$fname"
    "toolbox/$fname"
    "scripts/$fname"
  )

  local url tmp status
  tmp="$(mktemp)"; mv "$tmp" "${tmp}.py"; tmp="${tmp}.py"

  for p in "${try_paths[@]}"; do
    url="$repo_base/$p"
    status="$(curl -sS -L -w '%{http_code}' -o "$tmp" "$url" || echo 000)"
    if [[ "$status" == "200" && -s "$tmp" ]]; then
      echo "Running Python from GitHub: $p $*"
      python3 "$tmp" "$@"; local rc=$?
      rm -f "$tmp"
      return "$rc"
    fi
  done

  rm -f "$tmp"
  echo "ERROR: Could not fetch Python '$fname' from known paths." >&2
  return 2
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
gh_source static_map.sh
# Python is executed via gh_py_exec (not sourced)

# --- Color print fallbacks (if helper didn't provide them) ---
if ! declare -F PRINT_CYAN   >/dev/null; then PRINT_CYAN()   { printf "\033[36m%s\033[0m\n" "$*"; }; fi
if ! declare -F PRINT_YELLOW >/dev/null; then PRINT_YELLOW() { printf "\033[33m%s\033[0m\n" "$*"; }; fi
if ! declare -F PRINT_RED    >/dev/null; then PRINT_RED()    { printf "\033[31m%s\033[0m\n" "$*"; }; fi

# --- resolve script dir for local files like path_definition.txt ---
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
identity="$(whoami)@$(hostname)"

default_root="/Volumes/Extreme_Pro/fMRI"

# ---- prompt/CLI for root_location ----
root_location="${1:-}"
if [[ "${root_location:-}" == "--root" ]]; then
  shift
  root_location="${1:-}"
  shift || true
fi

if [[ -z "${root_location:-}" ]]; then
  echo
  read_default "Root location" "${default_root}" root_location_input
  root_location="${root_location_input}"
fi

# Expand ~, $VARS, and relative paths → absolute
root_location="$(expand_path "$root_location")"

if [[ ! -d "$root_location" ]]; then
  echo "ERROR: root location '$root_location' does not exist." >&2
  exit 1
fi
echo "Using root location: $root_location"

# ---- CSV config (use absolute path!) ----
csv="Animal_Experiments_Sequences.csv"
csv_path="$root_location/RawData/$csv"
header_lines=2
DELIM=$'\x1f'

if [[ ! -f "$csv_path" ]]; then
  echo "ERROR: $csv not found at $csv_path"
  exit 1
fi

# ---- show available rows (first 8 columns) ----
echo -e "\n== Available rows (CSV line | first 8 columns) =="
python3 - "$csv_path" "$header_lines" <<'PY_AVAIL' | { column -t -s'|' || cat; }
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
PY_AVAIL

# ---- ask for multiple line numbers / ranges ----
echo
echo "You can choose multiple CSV lines, e.g.: 5,7,10-12"
read -rp "Enter CSV LINE NUMBERS (or q to quit): " sel
case "$sel" in q|Q) echo "Aborted."; exit 0 ;; esac
[[ -n "$sel" ]] || { echo "No selection."; exit 1; }

# ---- expand comma/range list to a unique, sorted list of integers > header_lines ----
expand_lines() {
  local input="$1" part a b n
  IFS=',' read -r -a parts <<< "$input"
  for part in "${parts[@]}"; do
    part="${part//[[:space:]]/}"
    if [[ "$part" =~ ^[0-9]+-[0-9]+$ ]]; then
      a="${part%-*}"; b="${part#*-}"
      if (( a > b )); then n="$a"; a="$b"; b="$n"; fi
      for (( n=a; n<=b; n++ )); do echo "$n"; done
    elif [[ "$part" =~ ^[0-9]+$ ]]; then
      echo "$part"
    else
      echo "ERROR: invalid token '$part' in selection." >&2
      exit 1
    fi
  done | awk -v hdr="$header_lines" '($1>hdr){print $1}' | sort -n | uniq
}

# ---- macOS-friendly read into array (no mapfile on Bash 3.2) ----
LINES=()
if [[ -n "${BASH_VERSINFO:-}" && "${BASH_VERSINFO[0]}" -ge 4 ]]; then
  mapfile -t LINES < <(expand_lines "$sel")
else
  while IFS= read -r ln; do
    [[ -n "$ln" ]] && LINES+=("$ln")
  done < <(expand_lines "$sel")
fi
((${#LINES[@]})) || { echo "No valid data lines selected (remember header lines are 1..$header_lines)."; exit 1; }

# ---- helpers ----
is_int() { [[ "$1" =~ ^-?[0-9]+$ ]]; }
prompt_if_unset() {
  local __var="$1"; shift
  local __prompt="$1"; shift || true
  local __def="${1:-}"
  if [[ -z "${!__var:-}" ]]; then
    if [[ -n "$__def" ]]; then
      read -rp "$__prompt [$__def]: " __ans
      printf -v "$__var" "%s" "${__ans:-$__def}"
    else
      read -rp "$__prompt: " __ans
      printf -v "$__var" "%s" "$__ans"
    fi
  fi
}

# ---- function: process one CSV line number (SANDBOXED in subshell) ----
process_csv_line() {
  local line_no="$1"

  (
    PRINT_CYAN "=== Processing CSV line $line_no ==="

    local parsed
parsed="$(
python3 - "$csv_path" "$line_no" <<'PY_PARSE'
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
PY_PARSE
)"
    [[ -n "$parsed" ]] || { echo "Selected line $line_no is blank. Skipping."; exit 0; }

    local project_name sub_project_name dataset_name structural_name functional_name struc_coregistration baseline_duration injection_duration
    local injection_on_left_side injection_on_right_side
    get_n_field(){ echo "$parsed" | cut -d"$DELIM" -f"$1"; }
    trim(){ printf '%s' "$1" | xargs; }

    project_name=$(trim "$(get_n_field 3)")
    sub_project_name=$(trim "$(get_n_field 4)")
    dataset_name=$(trim "$(get_n_field 2)")
    structural_name=$(trim "$(get_n_field 5)")
    functional_name=$(trim "$(get_n_field 6)")
    struc_coregistration=$(trim "$(get_n_field 7)")
    baseline_duration=$(trim "$(get_n_field 8)")
    injection_duration=$(trim "$(get_n_field 9)")
    injection_on_left_side=$(trim "$(get_n_field 22)")
    injection_on_right_side=$(trim "$(get_n_field 23)")

    export Project_Name="$project_name"
    export Sub_project_Name="$sub_project_name"
    export Dataset_Name="$dataset_name"
    export structural_run="$structural_name"
    export run_number="$functional_name"
    export str_for_coreg="$struc_coregistration"
    export baseline_duration_in_min="$baseline_duration"
    export injection_duration_in_min="$injection_duration"
    export injected_liquid_on_left_side="$injection_on_left_side"
    export injected_liquid_on_right_side="$injection_on_right_side"

    echo -e "\n== Selection summary =="
    echo "CSV line:                                 $line_no"
    echo "Project_Name:                             $Project_Name"
    echo "Sub_project_Name:                         $Sub_project_Name"
    echo "Dataset_Name:                             $Dataset_Name"
    echo "First Structural Run:                     $structural_run"
    echo "Functional Run Number:                    $run_number"
    echo "Structural Data used for Coregistration:  $str_for_coreg"
    echo "Baseline Duration (in min):               $baseline_duration_in_min"
    echo "Injection Duration (in min):              $injection_duration_in_min"
    echo "Liquid injected on Left Side:             $injected_liquid_on_left_side"
    echo "Liquid injected on Right Side:            $injected_liquid_on_right_side"

    local Path_Raw_Data="$root_location/RawData/$project_name/$sub_project_name"
    local Path_Analysed_Data="$root_location/AnalysedData/$project_name/$sub_project_name/$Dataset_Name"

    local datapath
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

    # ---------------- STRUCTURAL ----------------
    FUNC_PARAM_EXTRACT "$datapath/$structural_run"
    : "${SequenceName:?FUNC_PARAM_EXTRACT did not set SequenceName}"

    local struct_dir="$Path_Analysed_Data/${structural_run}${SequenceName}"
    mkdir -p "$struct_dir"
    cd "$struct_dir"
    CHECK_FILE_EXISTENCE "$struct_dir" || true

    run_if_missing "anatomy.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$structural_run" "$datapath/$structural_run/method"
    cp -f G1_cp.nii.gz anatomy.nii.gz || echo "WARNING: G1_cp.nii.gz not found for structural; continuing."

    # ---------------- STRUCTURAL FOR COREGISTRATION ----------------
    FUNC_PARAM_EXTRACT "$datapath/$str_for_coreg"
    : "${SequenceName:?FUNC_PARAM_EXTRACT did not set SequenceName}"

    local struct_coreg_dir="$Path_Analysed_Data/${str_for_coreg}${SequenceName}"
    mkdir -p "$struct_coreg_dir"
    cd "$struct_coreg_dir"
    CHECK_FILE_EXISTENCE "$struct_coreg_dir" || true

    run_if_missing "anatomy.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$str_for_coreg" "$datapath/$str_for_coreg/method"
    cp -f G1_cp.nii.gz anatomy.nii.gz || echo "WARNING: G1_cp.nii.gz not found for coreg structural; continuing."

    # ---------------- FUNCTIONAL ----------------
    FUNC_PARAM_EXTRACT "$datapath/$run_number"
    : "${SequenceName:?FUNC_PARAM_EXTRACT did not set SequenceName}"

    local func_dir="$Path_Analysed_Data/${run_number}${SequenceName}"
    mkdir -p "$func_dir"
    cd "$func_dir"
    CHECK_FILE_EXISTENCE "$func_dir" || true

    run_if_missing "G1_cp.nii.gz" -- BRUKER_to_NIFTI "$datapath" "$run_number" "$datapath/$run_number/method"

    # ---------------- Motion Correction (Using AFNI) ----------------
    : "${MiddleVolume:?FUNC_PARAM_EXTRACT (or motion helper) did not set MiddleVolume}"
    PRINT_YELLOW "Performing Step 1: Motion Correction"
    run_if_missing "mc_func.nii.gz" "mc_func+orig.HEAD" "mc_func+orig.BRIK" -- \
      MOTION_CORRECTION "$MiddleVolume" G1_cp.nii.gz mc_func

    # ---------------- tSNR Estimation (Using AFNI) -----------------
    PRINT_YELLOW "Performing Step 2: Obtaining Mean func, Std func and tSNR Maps"
    run_if_missing "tSNR_mc_func.nii.gz" "tSNR_mc_func+orig.HEAD" "tSNR_mc_func+orig.BRIK" -- \
      TEMPORAL_SNR_using_FSL mc_func.nii.gz

    # ---------------- N4 Bias Field Correction ---------------------
    PRINT_YELLOW "Performing Step 3: Performing N4 Bias Field Correction of mean_mc_func"
    run_if_missing "cleaned_mc_func.nii.gz" -- \
      BIAS_CORRECTED_IMAGE mean_mc_func.nii.gz 100 mc_func.nii.gz

  # ---------------- Time Series Inspection ---------------------
    echo "Inspect Time Series of your data to decide baseline and signal periods."
    # Open viewer non-blocking, then pause explicitly
    fsleyes cleaned_mc_func.nii.gz &  # do not block
    read -rp "After reviewing in FSLeyes, press ENTER to continue..." _

    # Optional sanity check if run non-interactively
    if [[ ! -t 0 ]]; then
      echo "WARNING: stdin is not a TTY; interactive prompts may not appear."
    fi

    # ---------------- Static Map (Generation) ----------------------
    PRINT_YELLOW "Performing Step 4: Generating Static Map (SCM)"
    local map_candidates=(Static_SCM*.nii.gz)
    local map_file=""

    # Show what we see (helps debug)
    if ((${#map_candidates[@]})); then
      printf 'DEBUG: %s\n' "${map_candidates[@]}"
    fi

    # Always ask for indices; if a map exists, let the user choose reuse vs regenerate
    local base_start base_end sig_start sig_end

    prompt_if_unset base_start "Enter baseline start index" "100"
    prompt_if_unset base_end   "Enter baseline end index"   "300"
    prompt_if_unset sig_start  "Enter signal start index"   "301"
    prompt_if_unset sig_end    "Enter signal end index"     "500"

    for v in base_start base_end sig_start sig_end; do
      is_int "${!v}" || { echo "ERROR: $v must be an integer (got '${!v}')"; exit 1; }
    done

    if ((${#map_candidates[@]})); then
      read -rp "Static map exists (${map_candidates[0]}). Reuse it? [Y/n]: " reuse
      reuse="${reuse:-Y}"
      if [[ "$reuse" =~ ^[Yy]$ ]]; then
        map_file="${map_candidates[0]}"
        echo "Reusing existing Static Map: $map_file"
      else
        echo "Regenerating Static Map with your indices..."
        Static_Map cleaned_mc_func.nii.gz "$base_start" "$base_end" "$sig_start" "$sig_end"
      fi
    else
      echo "No Static_SCM*.nii.gz found. Generating now..."
      Static_Map cleaned_mc_func.nii.gz "$base_start" "$base_end" "$sig_start" "$sig_end"
    fi

    # Refresh candidates after (re)generation
    map_candidates=(Static_SCM*.nii.gz)
    if ((${#map_candidates[@]})); then
      map_file="${map_candidates[0]}"
      echo "Using Static Map: $map_file"
    else
      echo "WARNING: Static Map generation did not produce Static_SCM*.nii.gz"
    fi


    # ---------------- Coregistration (Using AFNI) ------------------
    PRINT_YELLOW "Performing Step 5: Coregistration of functional/static map to structural"
    local base_anat="$struct_coreg_dir/anatomy.nii.gz"

    3dAllineate \
      -base "$base_anat" \
      -input mean_mc_func.nii.gz \
      -1Dmatrix_save mean_func_struct_aligned.aff12.1D \
      -cost lpa \
      -prefix mean_func_struct_aligned.nii.gz \
      -1Dparam_save params.1D \
      -twopass

    if [[ -n "$map_file" ]]; then
      3dAllineate \
        -base "$base_anat" \
        -input "$map_file" \
        -1Dmatrix_apply mean_func_struct_aligned.aff12.1D \
        -master "$base_anat" \
        -final linear \
        -prefix Static_Map_coreg.nii.gz
    else
      echo "Skipping Static Map coregistration (no map_file)."
    fi

    # ---------------- Sliding-window Movie (Python) ----------------
    PRINT_YELLOW "Performing Step 6: Sliding-window static-map movie"
    # local tr_val="1.0"
    prompt_if_unset rotation_degree "By what degree do you want to rotate the brain from display?"
    # local out_prefix="Static_Fast"

    gh_py_exec make_static_maps_and_movie.py \
      cleaned_mc_func.nii.gz \
      --baseline-start 100 --baseline-end 300 \
      --underlay-first-n 600 --window 200 \
      --mode pos --vmax 15 --cmap blackbody \
      --fps 10 --tr 1.0 --rotate $rotation_degree \
      --out-prefix SCM_Map_movie_positive_only

    echo "✔ Completed pipeline for CSV line $line_no."
  )
}

# ---- iterate over all selected lines (each in its own subshell) ----
for ln in "${LINES[@]}"; do
  process_csv_line "$ln"
done

echo
echo "All requested CSV lines processed."
