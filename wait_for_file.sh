#PBS -l nodes=1:ppn=1,vmem=2g,mem=2g,walltime=24:00:00
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

# function to look for file
file_lookup(){
    until [[ -e $1 ]]; do
        # checks for the file every minute
        sleep 60
    done
    return 0
}
export -f file_lookup

# set a timeout for the the file lookup function
# timeout is set for 23 hours (in seconds)
timeout 82800 file_lookup $file

# check if commands completes
if [[ "$?" == 0 ]]; then
    qsub -v sample=${sample},mode=${mode} ${pipeline_dir}/${script}
else
    qsub -v file=${file},sample=${sample},mode=${mode},script=${script} ${pipeline_dir}/wait_for_file.sh
fi
