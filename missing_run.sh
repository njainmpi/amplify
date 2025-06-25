#!/bin/sh

run_if_missing() {
    local command_to_run
    local all_missing=true
    local files_to_check=()

    # Read arguments until we hit "--", the rest is the command
    while [[ "$1" != "--" && "$#" -gt 0 ]]; do
        files_to_check+=("$1")
        shift
    done

    # Skip the "--"
    shift

    command_to_run="$*"

    for file in "${files_to_check[@]}"; do
        if [ -f "$file" ]; then
            echo "Skipping: $file already exists."
            all_missing=false
            break
        fi
    done

    if $all_missing; then
        echo "None of ${files_to_check[*]} found. Running: $command_to_run"
        eval "$command_to_run"
    fi
}
