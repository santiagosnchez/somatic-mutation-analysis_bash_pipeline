#!/bin/bash
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

echo "submitting command:"
# submit jobs in parallel
parallel --dry-run --colsep="," qsub -v sample={},mode=${mode} ${pipeline_dir}/02_fastq_to_ubam.picard.Fastq2Sam.samtools.merge.sh ::: ${samples} 
parallel --colsep="," qsub -v sample={},mode=${mode} ${pipeline_dir}/02_fastq_to_ubam.picard.Fastq2Sam.samtools.merge.sh ::: ${samples} 


