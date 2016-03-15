#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: make.sh"
    echo "make script for this program"
    exit
fi
selfdir=$(cd $(dirname $0);pwd)
make clean
make main >/dev/null 2>&1