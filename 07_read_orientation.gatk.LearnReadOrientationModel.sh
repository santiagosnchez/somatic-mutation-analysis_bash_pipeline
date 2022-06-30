#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.gatk-learn-read-orientation.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
#module load gatk/4.2.2.0

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

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
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

# check if file exists
if [[ ! -e mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz ]]; then
# prepare input from scatter runs
all_f1_r2_input=$(ls mutect2/f1r2/${tumor}__${normal}.[1-9]*.f1r2.tar.gz | sed 's/^/-I /')

# run gatk's read orientation model
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" LearnReadOrientationModel $all_f1_r2_input -O mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz

else
# so that $? is 0
ls &> /dev/null
fi

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "07: ${tumor}__${normal} read-orientation analysis completed." | tee -a main.log
    # submit next step
    qsub -v \
tumor=${tumor},\
normal=${normal},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/08_filter_somatic_var.gatk.FilterMutectCalls.sh
    # calc runtime
    runtime=$( how_long "${start}" h )
    echo "07: Step ${tumor}__${normal}.gatk-learn-read-orientation.log took ${runtime} hours" | tee -a main.log
    # move logfiles if found
    if [[ -e ${tumor}__${normal}.gatk-learn-read-orientation.log ]]; then
        mv ${tumor}__${normal}.gatk-learn-read-orientation.log all_logfiles
        rm mutect2/f1r2/${tumor}__${normal}.[1-9]*.f1r2.tar.gz
    fi
fi
