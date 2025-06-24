#!/bin/bash

# Function to log script execution
log_execution() {
    local log_dir="$1"
    local script_name; script_name=$(basename "$0")

    if [[ -z "$log_dir" ]]; then
        printf "Error: Log directory path is not provided.\n" >&2
        return 1
    fi

    mkdir -p "$log_dir" || {
        printf "Error: Failed to create log directory '%s'.\n" "$log_dir" >&2
        return 1
    }

    local log_file="$log_dir/${script_name%.sh}.log"
    printf "[%s] Executed script: %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$script_name" >> "$log_file" || {
        printf "Error: Failed to write to log file '%s'.\n" "$log_file" >&2
        return 1
    }
}

# Function to log function execution
log_function_execution() {
    local log_dir="$1"
    local function_name="$2"
    local script_name; script_name=$(basename "$0")
    local user; user=$(whoami)

    if [[ -z "$log_dir" || -z "$function_name" ]]; then
        printf "Error: Log directory or function name is not provided.\n" >&2
        return 1
    fi

    local log_file="$log_dir/${script_name%.sh}.log"
    printf "[%s] Executed function: %s by user: %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$function_name" "$user" >> "$log_file" || {
        printf "Error: Failed to write to log file '%s'.\n" "$log_file" >&2
        return 1
    }
}

# # Main execution block
# main() {
#     log_execution "$LOG_DIR" || exit 1

#     # Example function calls directly in main
#     log_function_execution "$LOG_DIR" "Function_1" || exit 1
#     log_function_execution "$LOG_DIR" "Function_2" || exit 1
# }

# # Define the log directory
# LOG_DIR="$HOME/bash_script_logs"

# main "$@"
