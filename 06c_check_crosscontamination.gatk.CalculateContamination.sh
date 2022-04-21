#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.CalculateContamination.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create dir for contamination
if [[ ! -e contamination ]]; then
    mkdir contamination
fi
# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# create tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/export_paths_to_reference_files.sh ${organism} ${genome} ${mode}

if [[ "${normal}" == "PON" ]]; then
    echo "06: No CalculateContamination. Tumor-only mode." | tee -a main.log
    check_finish=0
else
  if [[ ! -e contamination/${tumor}__${normal}.calculatecontamination.table ]]; then
# run gatk's CalculateContamination
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" CalculateContamination \
 -I contamination/${tumor}.getpileupsummaries.table \
 -matched contamination/${normal}.getpileupsummaries.table \
 -O contamination/${tumor}__${normal}.calculatecontamination.table \
 --tumor-segmentation contamination/${tumor}__${normal}.tumorsegmentation.table
  else
    echo "06: CalculateContamination table found for ${tumor}__${normal}" | tee -a main.log
    check_finish=0
  fi
fi

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "06: ${tumor}__${normal} CalculateContamination completed." | tee -a main.log
    # move logfile
    mv ${tumor}__${normal}.CalculateContamination.log all_logfiles
fi
