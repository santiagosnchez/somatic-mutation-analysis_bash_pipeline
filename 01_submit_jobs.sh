#!/bin/bash

# print message and exit
die(){
    echo -e "$@"
    return 0
}
export -f die

# random 10 alphanumeric
rand(){
    alphanum=$(head /dev/urandom | tr -dc '[:alnum:]' | fold -10 | head -1)
    echo $alphanum
}
export -f rand

# find pipeline directory
default_pipe_dir=$(dirname $(readlink -f "${0}"))

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
                                  Default: hg38 (human), mm10 (mouse)
                                  Options: hg38, mm10, b37 (need to add those from other organisms)

--skip-alignment, -s    BOOL      Flag to skip the alignment steps. Go directly to
                                  variant calling.
                                  Default: false

--alignment-only, -l    BOOL      Just do alignment up to BQSR. Skip variant calling and onwards.
                                  Default: false

--pipeline, -p          STR       Specify a different source location for pipeline scripts.
                                  Useful for testing new or old versions.
                                  Default: ${default_pipe_dir}

--dry-run, -d           BOOL      Indicates if starting variables should be printed to screen only.
                                  Without submitting any jobs.

--fresh-start, -e       BOOL      Clears the working directory completely except for csv/text files.
                                  <<TAKE CAUTION>> using this argument.

--help, -h              BOOL      Print this help message.

Example commands:

# run whole exome data on human genome hg38
01_submit_jobs -m wes

# run whole exome data on mouse genome (default to mm10)
01_submit_jobs -m wes -o mouse

# run a whole genome sequencing analysis on human data
01_submit_jobs -m wgs

# add/append more samples to already started WES analysis
01_submit_jobs -m wes -f file_list.more_samples.csv -a

# run WES in \"alignment-only\" mode
01_submit_jobs -m wes --alignment-only

# run WES in \"variant-calling-only\" mode, i.e, --skip-alignment
# for a new analysis (i.e., bams are specified in the file_list.csv file)
01_submit_jobs -m wes --skip-alignment

# run WES in \"variant-calling-only\" mode, i.e, --skip-alignment
# from a previous run (i.e, bams are in either \"bam\" or \"BQSR\" directories)
01_submit_jobs -m wes -a --skip-alignment

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
    export pipeline_dir=${args[0]}
    #unset args[0]
    export skip_aln=0
    export aln_only=0
    export dry_run=0
    export fresh_start=0
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
                # export defaults
                if [[ "${organism}" == "mouse" ]]; then
                    export genome="mm10"
                fi
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
            elif [[ "${args[$i]}" == "-l" || "${args[$i]}" == "--alignment-only" ]]; then
                export aln_only=1
            elif [[ "${args[$i]}" == "-d" || "${args[$i]}" == "--dry-run" ]]; then
                export dry_run=1
            elif [[ "${args[$i]}" == "-e" || "${args[$i]}" == "--fresh-start" ]]; then
                export fresh_start=1
            elif [[ "${args[$i]}" == "-"* ]]; then
                echo "$help_message"
                die "Error: argument ${args[$i]} was not recognized." && return 1
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
read_and_export_arguments $default_pipe_dir $* || exit 1

# check if user wants a fresh start
if [[ ${fresh_start} == 1 ]]; then
    echo "Are you sure you want to delete these files? : "
    file * | grep -v "\.csv:\|\.txt:"
    echo -n "[y|n]? : "
    read -r response
    if [[ "${response}" == "y" || "${response}" == "Y" || "${response}" == "yes" || "${response}" == "YES" ]]; then
        if [[ $(qstat -u `whoami` | wc -l) > 0 ]]; then
            echo "stoping pipeline first..."
            ${pipeline_dir}/stop_pipeline
        fi
        rm -rf $(ls | grep -v "\.csv\|\.txt")
        echo "01: Cleared."
    else
        exit 0
    fi
fi

# if user wants to overwrite results
if [[ -e main.log && ${append} == 0 ]]; then
    echo -e "$help_message"
    echo -e "Error: It looks like the pipeline has already started. Run command with the -a flag to\nappend new samples (see examples)."
    exit 1
fi

# test required mode
if [[ ${#mode} == "" ]]; then
    echo -e "$help_message"
    echo "Error: -m/--mode is required. Select wes or wgs."
    exit 1
elif [[ ${mode} == "wes" || ${mode} == "wgs" ]]; then
    if [[ "$dry_run" == 0 ]]; then
        if [[ "$append" == 0 ]]; then
            # add date and mode to main.log
            date | tee main.log
        fi
        echo -e "\n01: Running as mode: ${mode}" | tee -a main.log
    else
        # print date and mode
        date
        echo -e "\n01: Running as mode: ${mode}"
    fi
else
    echo -e "$help_message"
    echo -e "Error: -m/--mode can only be \"wes\" or \"wgs\" (all lowercase)."
    exit 1
fi

# check for incompatible arguments
if [[ "${aln_only}" == 1 && ${skip_aln} == 1 ]]; then
    echo -e "$help_message"
    echo -e "Error: --skip-alignment [-s] and --alignment-only [-l] are incompatible arguments."
    exit 1
elif [[ "${aln_only}" == 1 ]]; then
    if [[ "$dry_run" == 0 ]]; then
        echo -e "\n01: Doing alignment only up to BQSR." | tee -a main.log
    else
        echo -e "\n01: Doing alignment only up to BQSR."
    fi
elif [[ "${skip_aln}" == 1 ]]; then
    if [[ "$dry_run" == 0 ]]; then
        echo -e "\n01: Skipping alignment step." | tee -a main.log
    else
        echo -e "\n01: Skipping alignment step."
    fi
fi

# check tumors_and_normals.csv
current=`pwd -P`
if [[ ! -e tumors_and_normals.csv ]]; then
    echo "tumors_and_normals.csv file not found in working directory: $current"
    exit 1
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode} || exit 1

# create .tmp dir (hidden)
if [[ ! -e ./.tmp ]]; then
    mkdir ./.tmp
fi

# add last commit version of pipeline
if [[ "$dry_run" == 0 ]]; then
    echo -e "\nLast git commit version:" | tee -a main.log
    cd $pipeline_dir && (git log | head -3) | tee -a $current/main.log
    cd -
else
    echo -e "\nLast git commit version:"
    cd $pipeline_dir && (git log | head -3)
    cd -
fi

# submit jobs in parallel
if [[ "$dry_run" == 0 ]]; then
    echo -e "\n01: Running arguments:" | tee -a main.log
    echo "organism: ${organism}" | tee -a main.log
    echo "reference: ${genome}" | tee -a main.log
    echo "pipeline path: ${pipeline_dir}" | tee -a main.log
    echo "file list: ${file_list}" | tee -a main.log
else
    echo -e "\n01: Running arguments:"
    echo "organism: ${organism}"
    echo "reference: ${genome}"
    echo "pipeline path: ${pipeline_dir}"
    echo "file list: ${file_list}"
fi

if [[ ${skip_aln} == 0 ]]; then
    # first record arguments
    if [[ "$dry_run" == 0 ]]; then
        cat ${file_list} | parallel --tmpdir ./.tmp --colsep="," '
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
    else
        cat ${file_list} | parallel --tmpdir ./.tmp --colsep="," '
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
    '
    fi

    njobs=$(cat ${file_list} | wc -l)
    # then submit
    if [[ "$dry_run" == 0 ]]; then
        echo -e "\n01: Submitting ${njobs} jobs now ..." | tee -a main.log
        cat ${file_list} | parallel --tmpdir ./.tmp --colsep="," '
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
genome=${genome},\
aln_only=${aln_only} \
${pipeline_dir}/02a_check_pairs.sh
    else
      echo "File not found:"
      echo {2}
    fi
    ' | tee -a main.log
    else
        echo -e "\n01: Submitting ${njobs} jobs now ..."
    fi
else
    # first check that file_list includes bam files
    # first record arguments
    # make bam dir
    if [[ -e BQSR ]]; then
        export dir=BQSR
    elif [[ -e bam ]]; then
        export dir=bam
    else
        mkdir bam
        export dir=bam
    fi
    bams=$(ls ${dir}/*.bqsr.bam 2> /dev/null || echo "")
    if [[ ${#bams[@]} -gt 0 ]]; then
        if [[ "$dry_run" == 0 ]]; then
            echo -e "\n01: bam dir is not empty. Overriding ${file_list}." | tee -a main.log
        else
            echo -e "\n01: bam dir is not empty. Overriding ${file_list}."
        fi
        file_list=`rand`
    fi
    # make symlinks for all bams if bam dir not empty
    if [[ ${dry_run} == 0 ]]; then
        cat ${file_list} 2> /dev/null || for bam in $bams; do
            echo $bam | perl -ne '$_ =~ m/\/(.*?)\.bqsr\.bam/; print "$1,$_"';
        done | parallel --plus --tmpdir ./.tmp --colsep="," '
        if [[ {2} == *".bam" && {2} != ${dir}/{1}.bqsr.bam ]]; then
          if [[ $(samtools quickcheck {2} && echo 1) == 1 ]]; then
            if [[ $(samtools view -H {2} | grep "SO:coordinate" &> /dev/null && echo 1) == 1 ]]; then
              echo -e -n "bam {2/} is sorted.\nLinking as "
              ln -s {2} bam/{1}.bqsr.bam ;
              echo "bam/{1}.bqsr.bam";
              if [[ -e {2.}.bai ]]; then
                ln -s {2.}.bai bam/{1}.bqsr.bai ;
                echo "bam/{1}.bqsr.bam";
              elif [[ -e {2}.bai ]]; then
                ln -s {2}.bai bam/{1}.bqsr.bai ;
                echo "bam/{1}.bqsr.bam";
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
        else
            if [[ {2} == ${dir}/{1}.bqsr.bam ]]; then
                if [[ -e {2.}.bai ]]; then
                    echo "bam {2} is indexed."
                elif [[ -e {2}.bai ]]; then
                    echo "bam {2} is indexed."
                else
                    echo "bam {2} is not indexed."
                fi
            else
                echo "{2} not a bam file."
            fi
        fi
        ' | tee -a main.log

        njobs=$(ls $dir/*.bam | wc -l)
        # then submit
        echo -e "\n01: Submitting ${njobs} jobs now ..." | tee -a main.log
        ready=$(ls $dir | grep ".bam$" | sed 's/\..*//')
        if [[ ${#ready} == 0 ]]; then
          echo "bam files are not ready to run. Check working paths and files."
          exit 1
        else
          for sample in $ready; do
            TN=$(grep "^${sample},\|,${sample}$" tumors_and_normals.csv)
            if [[ ${#TN} -gt 0 ]]; then
              wt=$(get_walltime $(readlink -f $dir/${sample}.bqsr.bam))
              qsub -l walltime=${wt}:00:00 -v \
sample=${sample},\
wt=${wt},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome},\
aln_only=0 \
${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh | tee -a main.log
            fi
          done
        fi
    else
      # first check that file_list includes bam files
      # first record arguments
      # make symlinks for all bams
      cat ${file_list} 2> /dev/null || for bam in $bams; do
          echo $bam | perl -ne '$_ =~ m/\/(.*?)\.bqsr\.bam/; print "$1,$_"';
      done | parallel --plus --tmpdir ./.tmp --colsep="," '
      if [[ {2} == *".bam" && {2} != ${dir}/{1}.bqsr.bam ]]; then
        if [[ $(samtools quickcheck {2} && echo 1) == 1 ]]; then
          if [[ $(samtools view -H {2} | grep "SO:coordinate" &> /dev/null && echo 1) == 1 ]]; then
            echo "bam {2/} is sorted. Linking..."
            if [[ -e {2.}.bai ]]; then
              echo "bam {2/} is indexed as {2./}.bai"
            elif [[ -e {2}.bai ]]; then
              echo "bam {2/} is indexed as {2/}.bai"
            else
              echo "bam index not found for {2/}"
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
      else
        if [[ {2} == ${dir}/{1}.bqsr.bam ]]; then
            if [[ -e {2.}.bai ]]; then
                echo "bam {2} is indexed."
            elif [[ -e {2}.bai ]]; then
                echo "bam {2} is indexed."
            else
                echo "bam {2} is not indexed."
            fi
        else
            echo "{2} not a bam file."
        fi
      fi
      '

      njobs=$(ls $dir/*.bam | wc -l)
      # then submit
      echo -e "\n01: Submitting ${njobs} jobs now ..."
    fi
fi

if [[ "$?" != 0 ]]; then
    how_far_away=$(cat main.log | grep -o "^[0-9][0-9]: " | sort -u | wc -l)
    if [[ ${how_far_away} == 1 ]]; then
        echo "01: Deleting log and reverting back to starting conditions"
        rm -rf main.log ./.tmp
    fi
fi
