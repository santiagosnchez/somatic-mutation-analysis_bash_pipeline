#!/bin/bash

# print message and exit
die(){
    echo -e "$@"
    return 0
}
export -f die

help_message="

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Somatic and Germline Mutation Discovery Pipeline
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

©Santiago Sánchez-Ramírez, SickKids
Contact: santiago.sanchezramirez@sickkids.ca

Required arguments:

--mode, -m              STR       Data mode. Options are: wes, wgs

Optional arguments:

--file_list, -f         STR       3-column CSV file with paths to FASTQ files.
                                  First column is the sample name. Second is the path
                                  to the forward reads (R1.fastq.gz). Third is the path
                                  to the reverse reads (R2.fastq.gz).
                                  Default: file_list.csv

--organism, -o          STR       Organismal data source. Needed to adjust the genome
                                  reference files. Options are: human, mouse
                                  Default: human

--append, -a            BOOL      Flag to append to current run.
                                  Default: false

--reference, -r         STR       Provide a the version of the reference to use.
                                  Default: hg38
                                  Options: hg38, hs37d5, b37, (need to add those from other organisms)

--skip-alignment, -s    BOOL      Flag to skip the alignment steps. Go directly to
                                  variant calling.
                                  Default: false

--pipeline, -p          STR       Specify a different source location for pipeline scripts.
                                  Useful for testing new or old versions.
                                  Default: /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery

--help, -h              BOOL      Print this help message.

Example commands:

# run whole exome data on human genome hg38
01_submit_jobs -m wes

# run whole exome data on mouse genome
01_submit_jobs -m wes -o mouse

# run a whole genome sequencing analysis on human data
01_submit_jobs -m wgs

# add/append more samples to already started WES analysis
01_submit_jobs -m wes -f file_list.more_samples.csv -a

"
export help_message

# function that reads and then exports arguments
read_and_export_arguments(){
    args=($@)
    # default arguments
    export organism="human"
    export file_list="file_list.csv"
    export append=0
    export genome="hg38"
    export pipeline_dir="/hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery"
    export skip_aln=0
    # required
    export mode=""
    # loop through arguments if there are any
    if [[ ${#args[@]} > 0 ]]; then
        for i in `seq 1 ${#args[@]}`; do
            i=$(( i - 1 ))
            if [[ "${args[$i]}" == "-h" || "${args[$i]}" == "--help" || "${args[$i]}" == "-help" ]]; then
                die "$help_message" && return 1
            elif [[ "${args[$i]}" == "-m" || "${args[$i]}" == "--mode" ]]; then
                export mode=${args[$(( i + 1 ))]}
            elif [[ "${args[$i]}" == "-o" || "${args[$i]}" == "--organism" ]]; then
                export organism=${args[$(( i + 1 ))]}
            elif [[ "${args[$i]}" == "-f" || "${args[$i]}" == "--file_list" ]]; then
                export file_list=${args[$(( i + 1 ))]}
            elif [[ "${args[$i]}" == "-r" || "${args[$i]}" == "--reference" ]]; then
                export genome=${args[$(( i + 1 ))]}
            elif [[ "${args[$i]}" == "-p" || "${args[$i]}" == "--pipeline" ]]; then
                export pipeline_dir=${args[$(( i + 1 ))]}
            elif [[ "${args[$i]}" == "-a" || "${args[$i]}" == "--append" ]]; then
                export append=1
            elif [[ "${args[$i]}" == "-s" || "${args[$i]}" == "--skip-alignment" ]]; then
                export skip_aln=1
            fi
        done
        return 0
    else
        die "$help_message" && return 1
    fi
}

# load modules
module load parallel/20210322
module load samtools/1.10

# read all arguments
read_and_export_arguments $* || exit 1

# if user wants to overwrite results
if [[ -e main.log && ${append} == 0 ]]; then
    echo -e "$help_message"
    echo -e "Error: It looks like the pipeline has already started. Run command with the -a flag to\nappend new samples (see examples)."
    exit 1
fi

# print date
if [[ "$append" == 0 ]]; then
    # add date to first line
    date | tee main.log
fi

# test required mode
if [[ ${#mode} == "" ]]; then
    echo -e "$help_message"
    echo "Error: -m/--mode is required. Select wes or wgs."
    exit 1
elif [[ ${mode} == "wes" || ${mode} == "wgs" ]]; then
    echo "Running as mode: ${mode}"
else
    echo -e "$help_message"
    echo -e "Error: -m/--mode can only be \"wes\" or \"wgs\" (all lowercase)."
    exit 1
fi

# check tumors_and_normals.csv
if [[ ! -e tumors_and_normals.csv ]]; then
    current=`pwd -P`
    echo "tumors_and_normals.csv file not found in working directory: $current"
    exit 1
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/export_paths_to_reference_files.sh ${organism} ${genome} ${mode} || exit 1

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# add last commit version of pipeline
echo -e "\nLast git commit version:" | tee -a main.log
cd $pipeline_dir && (git log | head -3) && cd - | tee -a main.log

# submit jobs in parallel
echo -e "\n01: Running arguments:" | tee -a main.log
echo "organism: ${organism}" | tee -a main.log
echo "reference: ${genome}" | tee -a main.log
echo "pipeline path: ${pipeline_dir}" | tee -a main.log
echo "file list: ${file_list}" | tee -a main.log

if [[ ${skip_aln} == 0 ]]; then
    # first record arguments
    cat ${file_list} | parallel --tmpdir ./tmp --colsep="," '
  if [[ -e {2} ]]; then
    wt=$(get_walltime {2} {3});
    echo "sample: {1}";
    echo "R1: {2}";
    echo "R2: {3}";
    echo "walltime: ${wt}";
    echo "index: {#}";
  else
    echo "File not found:";
    echo {2};
  fi
' | tee -a main.log

    njobs=$(cat ${file_list} | wc -l)
    # then submit
    echo -e "\n01: Submitting ${njobs} jobs now ..." | tee -a main.log
    cat ${file_list} | parallel --tmpdir ./tmp --colsep="," '
if [[ -e {2} ]]; then
  wt=$(get_walltime {2} {3});
  rg=`get_read_group_info {2} {1}`;
  qsub -l walltime="${wt}":00:00 -v \
wt="${wt}",\
file_list=${file_list},\
index={#},\
sample={1},\
forward={2},\
reverse={3},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/02a_check_pairs.sh
else
  echo "File not found:"
  echo {2}
fi
' | tee -a main.log
else
    echo -e "\n01: Skipping alignment. Moving directly to variant calling."
    # first check that file_list includes bam files
    # first record arguments
    # make bam dir
    if [[ ! -e bam ]]; then
        mkdir bam
    fi
    # make symlinks for all bams
    cat ${file_list} | parallel --tmpdir ./tmp --colsep="," '
    if [[ {2} == *".bam" ]]; then
      if [[ $(samtools quickcheck {2}) && echo 1) == 1 ]]; then
        if [[ $(samtools view -H {2} | grep "SO:coordinate" &> /dev/null && echo 1) == 1 ]]; then
          echo "bam {2} is sorted. Linking..."
          ln -s {2} bam/{1}.bqsr.bam ;
          if [[ -e {2.}.bai ]]; then
            ln -s {2.}.bai bam/{1}.bqsr.bai ;
          elif [[ -e {2}.bai ]];
            ln -s {2}.bai bam/{1}.bqsr.bai ;
          else
            echo "bam index not found for {2}"
          fi
          wt=$(get_walltime {2});
          echo "sample: {1}";
          echo "bam: {2}";
          echo "walltime: ${wt}";
        else
          echo "bam {2} is unsorted. Omitting..."
        fi
      else
        echo "samtools quickcheck failed for bam {2}"
      fi
      ' | tee -a main.log

    njobs=$(cat ${file_list} | wc -l)
    # then submit
    echo -e "\n01: Submitting ${njobs} jobs now ..." | tee -a main.log
    ready=$(ls bam | grep ".bam$" | sed 's/\..*//')
    if [[ ${#ready} == 0 ]]; then
      echo "bam files are not ready to run. Check working paths and files."
      exit 1
    else
      for sample in "$ready"; do
        TN=$(grep "${sample}," tumors_and_normals.csv)
        if [[ ${#TN} -gt 0 ]]; then
          tumor=$(echo $TN | sed 's/,.*//')
          wt=$(get_walltime $(readlink -f bam/${tumor}.bqsr.bam))
          qsub -l walltime=${wt}:00:00 -v \
sample=${sample},\
wt=${wt},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh | tee -a main.log
        fi
      done
    fi
fi

if [[ "$?" != 0 ]]; then
    how_far_away=$(cat main.log | grep -o "^[0-9][0-9]: " | sort -u | wc -l)
    if [[ ${how_far_away} == 1 ]]; then
        echo "01: Deleting log and reverting back to starting conditions"
        rm -rf main.log ./tmp
    fi
fi
