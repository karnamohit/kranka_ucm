#!/bin/bash

locate gauss_hf.py | cat > junk

IFS=$'\n'; file_array=($(<junk)); declare -p file_array | for i in "${file_array[@]}"; do cp "${file_array[-1]}" "$i"; done

diff "/mnt/d/Box Sync/Project Isborn/research_projects/tddft_ml/ml_tests/lih/tdci/gauss_hf.py" /mnt/d/GitHub/kranka_ucm/scripts/gauss_hf.py
