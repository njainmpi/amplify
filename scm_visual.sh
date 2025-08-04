#!/bin/bash

# Terminal colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
GRAY='\033[1;37m'
BOLD='\033[1m'
NC='\033[0m'

draw_case () {
    local label=$1
    local total_min=$2
    local discard_min=$3
    local baseline_min=$4
    local injection_min=$5
    local TR_sec=$6
    local Sequence=$7


    local post_min=$((total_min - baseline_min - injection_min))
    # local total_vol=$(( (total_min * 60) / TR_sec ))
    baseline_vol=$(( (baseline_min * 60) / TR_sec ))
    local injection_vol=$(( (injection_min * 60) / TR_sec ))
    local post_vol=$(( (post_min * 60) / TR_sec ))

    # Table Header
    echo ""
    echo "Data Acquired using $Sequence with TR of $TR_sec sec"
    echo ""
    printf "%-25s | %-20s %-20s %-20s\n" " " "Baseline" "Injection" "Post-Injection"
    echo "----------------------------------------------------------------------------------------"
    printf "%-25s | %-20s %-20s %-20s\n" "Time (in minutes)" "$baseline_min" "$injection_min" "$post_min"
    printf "%-25s | %-20s %-20s %-20s\n" "Volumes/Repetitions" "$baseline_vol" "$injection_vol" "$post_vol"
    echo "----------------------------------------------------------------------------------------"

    # Visual Bar
    echo -e "\nVisual Timeline for Data Acquisition (1 block = 1 min):"
    echo -n "        "

    for ((i=0; i<total_min; i++)); do
        if (( i >= 0 && i < baseline_min )); then
            echo -ne "${BLUE}■${NC}" # Entire baseline is retained
        elif (( i >= baseline_min && i < baseline_min + injection_min )); then
            echo -ne "${YELLOW}${BOLD}■${NC}" # Injection
        else
            echo -ne "${GREEN}■${NC}" # Post-injection
        fi
    done

    echo ""
    echo ""


    # Visual Bar
    echo -e "\nVisual Timeline for Data Analysis (1 block = 1 min):"
    echo -n "        "

    for ((i=0; i<total_min; i++)); do
        if (( i < discard_min )) || (( i >= baseline_min - discard_min && i < baseline_min )); then
            echo -ne "${RED}■${NC}"   # Discarded baseline
        elif (( i >= discard_min && i < baseline_min - discard_min )); then
            echo -ne "${GRAY}■${NC}" # Retained baseline
        elif (( i >= baseline_min && i < baseline_min + injection_min )); then
            echo -ne "${YELLOW}${BOLD}■${NC}" # Injection
        else
            echo -ne "${GREEN}■${NC}" # Post-injection
        fi
    done


    echo -e "\n"
    echo -e "Legend:" 
    echo ""
    echo -e "${RED}■${NC}=Discarded"
    echo -e "${BLUE}■${NC}=Baseline Acquired"
    echo -e "${GREEN}■${NC}=Retained"
    echo -e "${YELLOW}${BOLD}■${NC}=Injection"
    echo -e "${GRAY}${BOLD}■${NC}=Baseline for Analysis"
    echo ""
}

# CASE 1: TR = 17s, 180 volumes → 51 minutes
draw_case "CASE 1 (TR = 17s)" 51 1 20 10 1
