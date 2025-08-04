#!/bin/bash

# Define terminal color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
GRAY='\033[1;37m'
BOLD='\033[1m'
NC='\033[0m'  # No Color



PRINT_RED () {
    label=$1
    echo ""
    echo ""
    echo -e "${RED}${label}${NC}"
    echo ""
    echo ""
}


PRINT_YELLOW () {
    label=$1
    echo ""
    echo ""
    echo -e "${YELLOW}${label}${NC}"
    echo ""
    echo ""
}



PRINT_GREEN () {
    label=$1
    echo ""
    echo ""
    echo -e "${GREEN}${label}${NC}"
    echo ""
    echo ""
}