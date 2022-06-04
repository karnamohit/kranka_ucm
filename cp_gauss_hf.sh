#!/bin/bash

locate gauss_hf.py | cat > junk

IFS=$'\n'; file_array=($(<junk)); declare -p file_array | for i in "${file_array[@]}"; do cp "${file_array[-1]}" "$i"; done
