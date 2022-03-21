#PBS -l nodes=1:ppn=1,vmem=2g,mem=2g,walltime=6:00:00
#PBS -e ${sample}.waitforfile.log
#PBS -j eo
# scheduler settings

############### INFO #################
#
# Integrated Bash Pipeline for
# Somatic Mutation Discovery
#
# Author: Santiago Sanchez-Ramirez
# Year: 2021
# Email: santiago.snchez@gmail.com
#
# More info on README.md
#
# Notes:
# (1) Uses nextgenmap and samtools
#
#####################################

########################
# wait for file script #
########################

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/export_paths_to_reference_files.sh ${organism} ${genome} ${mode}

# set a timeout for the the file lookup function
# timeout is set for 5.5 hours (in seconds)
# file_lookup is a function in export_paths_to_reference_files
timeout 19800 bash -c "file_lookup $file"

# check if commands completes
if [[ "$?" == 0 ]]; then
  if [[ -z $tumor ]]; then
    qsub -v \
sample=${sample},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/${script}
  else
    qsub -v \
tumor=${tumor},\
normal=${normal},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/${script}
  fi
  mv ${sample}.waitforfile.log all_logfiles
else
  if [[ -z $tumor ]]; then
    qsub -v \
file=${file},\
sample=${sample},\
script=${script},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/wait_for_file.sh
  else
    qsub -v \
file=${file},\
sample=${sample},\
tumor=${tumor},\
normal=${normal},\
script=${script},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/wait_for_file.sh
  fi
fi
