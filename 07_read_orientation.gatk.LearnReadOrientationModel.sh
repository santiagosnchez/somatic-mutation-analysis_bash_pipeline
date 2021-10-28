#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.read-orientation.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.0.1.2

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load all paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# check if file exists
if [[ ! -e mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz ]]; then
# prepare input from scatter runs
all_f1_r2_input=$(ls mutect2/f1r2/${tumor}__${normal}.[1-9]*.f1r2.tar.gz | sed 's/^/-I /')

# run gatk's read orientation model
gatk-4.2.2.0 LearnReadOrientationModel $all_f1_r2_input -O mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz

else
# so that $? is 0
ls &> /dev/null
fi

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # mv logfiles if found
    if [[ ! -e all_logfiles/${tumor}__${normal}.read-orientation.log ]]; then
        mv ${tumor}__${normal}.read-orientation.log all_logfiles
        rm mutect2/f1r2/${tumor}__${normal}.[1-9]*.f1r2.tar.gz
    fi
    # submit next step
    qsub -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/08_filter_somatic_var.gatk.FilterMutectCalls.sh  
fi


