#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=5:00:00
#PBS -e ${sample}.ReadOrientation.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load python/3.7.7

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create dir for contamination
if [[ ! -e orientation ]]; then
    mkdir orientation
fi

# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

# check if python libs can be found
python -c '
import sys
try:
  import pysam
except:
  print("pysam lib cannot be found.")
  sys.exit(1)
else:
  print("pysam lib found.")
'

if [[ "$?" == 0 ]]; then
    # run code
    python ${pipeline_dir}/scripts/bam_read_orientation_in_bed.py ${dir}/${sample}.bqsr.bam ${intervals_bed} orientation/${sample}
    # store sig
    check_finish=$?
else
    check_finish=1
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "06: calculated read orientation counts for ${sample}." | tee -a main.log
    # calc runtime
    runtime=$( how_long "${start}" h )
    echo "06: Step ${sample}.ReadOrientation.log took ${runtime} hours" | tee -a main.log
    # move logfile
    mv ${sample}.ReadOrientation.log all_logfiles
fi
