#!/bin/bash

module load parallel/20210322

# check/find out mode
# first argument should be either wes or wgs
if [[ $# == 0 ]]; then
    echo "are sample \"wes\" or \"wgs\"?"
    echo -n "enter here: "
    read -r mode
else
    if [[ "$1" = *"wes"* ]]; then
        mode="wes"
    elif [[ "$1" = *"wgs"* ]]; then
        mode="wgs"
    else
        echo "mode not recognized..."
        echo -n "enter \"wes\" or \"wgs\" here: "
        read -r mode
    fi
fi

# check if mode var is correct
if [[ "${mode}" = *"wes"* || "${mode}" = *"wgs"* ]]; then
    echo "mode is ${mode}"
else
    echo "mode is still not recognized..."
    echo -n "enter \"wes\" or \"wgs\" here: "
    read -r mode
    if [[ "${mode}" != "wes" || "${mode}" != "wgs" ]]; then
        bash 1_submit_jobs.sh
    fi
fi

# extract samples from csv file
samples=$(cat file_list.csv | cut -d, -f1 | sort -u)
export mode

# load all paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# print date
date

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi 

echo "submitting command:"
# submit jobs in parallel
cat file_list.csv | parallel --dry-run --colsep="," 'wt=$(get_walltime {4} {5}); qsub -l walltime=${wt}:00:00 -v sample={1},lane={2},lib={3},forward={4},reverse={5},mode=${mode} ${pipeline_dir}/02_align_and_sort_bam_to_ref.bwa.sh' | tee start.log
cat file_list.csv | parallel --colsep="," 'wt=$(get_walltime {4} {5}); qsub -l walltime=${wt}:00:00 -v sample={1},lane={2},lib={3},forward={4},reverse={5},mode=${mode} ${pipeline_dir}/02_align_and_sort_bam_to_ref.bwa.sh'

mv start.log all_logfiles


