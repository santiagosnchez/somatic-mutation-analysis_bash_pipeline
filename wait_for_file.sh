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

# load pipeline path
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh

# set a timeout for the the file lookup function
# timeout is set for 5.5 hours (in seconds)
# file_lookup is a function in export_paths_to_reference_files
timeout 19800 bash -c "file_lookup $file"

# check if commands completes
if [[ "$?" == 0 ]]; then
  if [[ -z $sample ]]; then
    qsub -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/${script}
  else
    qsub -v sample=${sample},mode=${mode} ${pipeline_dir}/${script}
  fi
else
  if [[ -z $sample ]]; then
    qsub -v file=${file},tumor=${tumor},normal=${normal},mode=${mode},script=${script} ${pipeline_dir}/wait_for_file.sh
  else
    qsub -v file=${file},sample=${sample},mode=${mode},script=${script} ${pipeline_dir}/wait_for_file.sh
  fi
fi
