#!/bin/bash
dir=$1
for file in `ls $dir`; do 
    sample=$(echo $file | sed 's/\..*//')
    size=$(du -s $dir/$file | cut -f1)
    walltime=$(echo "((($size / 150000)*10)/60)" | bc)
    echo "${sample},${walltime}"
done
