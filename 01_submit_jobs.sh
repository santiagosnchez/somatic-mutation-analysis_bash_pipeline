#!/bin/bash

#####################
# example:
# bash ~/pipline/01_submit_jobs.sh wes
# or
# bash ~/pipline/01_submit_jobs.sh wes added_file_list.csv
# where "added_file_list.csv" is a list with additional samples
# to add
#####################

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

# export mode
export mode

# submit default default file list or specific list
if [[ ! -z $2 ]]; then
    file_list=$2
else
    file_list="file_list.csv"
fi

# append to existing run
append=0

# as if user wants to overwrite results
if [[ -e main.log ]]; then
    echo "It looks like the pipeline has already started..."
    echo -n "Do you want to rerun (r), append (a), or quit (q)? (rerun will overwrite results) [r|a|n]:"
    read -r response
    if [[ "${response}" == "r"* ]]; then
        # quit running jobs
        jobids=$(head -1 $(ls *.log | grep -v "main") | grep -o "^[1-9]*$")
        qdel ${jobids}
        # delete logfiles
        rm -rf all_logfiles *.log
    elif [[ "${response}" == "a"* ]]; then
        append=1
    else
        echo "Exiting..."
        exit 0
    fi
fi


# load all paths
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# print date
if [[ "$append" == 0 ]]; then
    date > main.log
fi

# export file list var
export file_list

echo "01: submitting command:" | tee -a main.log
# submit jobs in parallel
# first dry run
cat ${file_list} | parallel --tmpdir ./tmp --dry-run --colsep="," '
01: wt=$(get_walltime {2} {3});
qsub -l walltime="${wt}":00:00 -v file_list=${file_list},index={#},sample={1},forward={2},reverse={3},mode=${mode} ${pipeline_dir}/02a_check_pairs.sh' | tee -a main.log
# then submit
echo "submitting ..." | tee -a main.log
cat ${file_list} | parallel --tmpdir ./tmp --colsep="," '
wt=$(get_walltime {2} {3});
rg=`get_read_group_info {2} {1}`;
qsub -l walltime="${wt}":00:00 -v file_list=${file_list},index={#},sample={1},forward={2},reverse={3},mode=${mode} ${pipeline_dir}/02a_check_pairs.sh' | tee -a main.log
