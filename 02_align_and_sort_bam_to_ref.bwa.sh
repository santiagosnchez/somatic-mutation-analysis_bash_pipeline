#PBS -l nodes=1:ppn=10,vmem=60g,mem=60g
#PBS -e ${sample}.${lane}.bwa.log
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

#############
# script 3  #
#############

# load modules
module load sambamba/0.7.0
module load samtools/1.10
module load bwa/0.7.17

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID
# get walltime
wt=$(qstat -f $PBS_JOBID | sed -rn 's/.*Resource_List.walltime = (.*)/\1/p' | sed 's/:.*//')

# write job details to log
qstat -f $PBS_JOBID > ${sample}.${lane}.bwa.log

# if forward and reverse not defined
# search in file
# if [[ -z "${forward}" ]]; then
#    forward=$(cat file_list.csv | grep "${sample}," | cut -d, -f4)
#    reverse=$(cat file_list.csv | grep "${sample},${lane}," | cut -d, -f5)
#    lib=$(cat file_list.csv | grep "${sample},${lane}," | cut -d, -f3)
#    lane=$(cat file_list.csv | grep "${sample},${lane}," | cut -d, -f3)
# fi

# load all paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# create dir for ubam
if [[ ! -e aligned_bam ]]; then
    mkdir aligned_bam
fi

# check if prev step finished correctly
if [[ -e "aligned_bam/${sample}.${lane}.bam" ]]; then
    if [[ $(samtools quickcheck aligned_bam/${sample}.${lane}.bam && echo 1) != 1 ]]; then
       echo "resubmitting step and increase time by 2 hrs"
       wt=$(( wt + 2 ))
       rm aligned_bam/${sample}.${lane}.bam 
       qsub -l walltime=${wt}:00:00 -v sample={1},lane={2},lib={3},forward={4},reverse={5},mode=${mode} ${pipeline_dir}/02_align_and_sort_bam_to_ref.bwa.sh
       exit 0
    else
       # check if bam is sorted
       if [[ $(samtools quickcheck aligned_bam/${sample}.${lane}.sorted.bam && echo 1) != 1 ]]; then
           sambamba sort \
            -m 5GB \
            -t 10 \
            -o aligned_bam/${sample}.${lane}.sorted.bam \
               aligned_bam/${sample}.${lane}.bam
       else
           ls &> /dev/null
       fi
    fi
else  
   # check if input files exist and run aligner
   bwa mem \
    -t 10 \
    -R "@RG\tID:${lane}\tSM:${sample}\tLB:${lib}\tPL:Illumina" \
    $reference \
    ${forward} ${reverse} \
   | sambamba view \
     -S -f bam \
     -o aligned_bam/${sample}.${lane}.bam /dev/stdin        
   # sort file
   if [[ "$?" == 0 ]]; then
       sambamba sort \
         -m 5GB \
         -t 10 \
         -o aligned_bam/${sample}.${lane}.sorted.bam \
         aligned_bam/${sample}.${lane}.bam
   fi
   # delete unsorted
   if [[ "$?" == 0 ]]; then
       rm aligned_bam/${sample}.${lane}.bam
   fi
fi

# check finish
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# if finished successfuly, submit next job
if [[ "$check_finish" == 0 ]]; then
    # calculate the number of bam files to merge
    expected_bams=$(cat file_list.csv | grep -c "^${sample},")
    found_sorted_bams=$(ls aligned_bam/${sample}.*.sorted.bam | wc -l)
    # check integrity of all files
    samtools quickcheck aligned_bam/${sample}.*.sorted.bam
    # check if all files are done
    if [[ "$?" == 0 && "${expected_bams}" == "${found_sorted_bams}" ]]; then 
        # submit next job
        # can switch this to picards MarkDuplicate method
        qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/03_merge_bams.sambamba.sh
    fi
    # move logfile
    mv ${sample}.${lane}.bwa.log all_logfiles    
fi

