#PBS -l nodes=1:ppn=10,vmem=30g,mem=30g,walltime=5:00:00
#PBS -e ${sample}.sambamba-merge.log
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
module load java/1.8
#module load gatk/4.2.2.0

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

# check all
all_check=0

# create .tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi
# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check for file integrity
# submit previous step if currupt
for bam in aligned_bam/${sample}.*.sorted.bam; do
    if [[ $(samtools quickcheck $bam && echo 1) != 1 ]]; then
       echo "resubmitting previous step and increase time by 2hrs"
       wt=$(( wt + 2 ))
       qsub -l walltime=${wt}:00:00 -v \
wt=${wt},\
sample=${sample},\
forward=${forward},\
reverse=${reverse},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome},\
aln_only=${aln_only} \
${pipeline_dir}/02b_align_bam_to_ref.bwa.sh
       all_check=1
    fi
done

# check if one if the files did
# not pass the check
if [[ "${all_check}" == 1 ]]; then
   echo "03: Resubmitting 02 with +2 walltime... "
   exit 0
else
    # check if input files exist and merge
   if [[ ! -e aligned_bam/${sample}.merged.bam ]]; then
       # check how many files
       if [[ $(ls aligned_bam/${sample}.*.sorted.bam | wc -l) == 1 ]]; then
           mv aligned_bam/${sample}.*.bam aligned_bam/${sample}.merged.bam
       else
           sambamba merge -t 10 /dev/stdout aligned_bam/${sample}.*.sorted.bam \
           | $gatk_path/gatk SetNmMdAndUqTags \
            -I /dev/stdin \
            -O aligned_bam/${sample}.merged.bam \
            -R $reference \
            --CREATE_INDEX
       fi
   else
      if [[ $(samtools quickcheck aligned_bam/${sample}.merged.bam && echo 1) != 1 ]]; then
          echo "03: Resubmitting 03 ..."
          wt=$(( wt + 2 ))
          qsub -l walltime=${wt}:00:00 -v \
wt=${wt},\
sample=${sample},\
forward=${forward},\
reverse=${reverse},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome},\
aln_only=${aln_only} \
${pipeline_dir}/03_merge_bams.sambamba.sh
      else
         ls &> /dev/null
      fi
   fi
fi

# check finish
check_finish=$?

# if finished successfuly, submit next job
if [[ "$check_finish" == 0 ]]; then
    # delete scattered bams
    rm aligned_bam/${sample}.*.sorted.bam*
    # delete temporary files
    if [[ -e .tmp/${sample}_file_list.csv ]]; then
        rm .tmp/${sample}_file_list.csv
    fi
    # get new walltime
    wt=$(get_walltime aligned_bam/${sample}.merged.bam)
    # submit next job
    # can switch this to picards MarkDuplicate method
    qsub -l walltime=${wt}:00:00 -v \
sample=${sample},\
wt=${wt},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome},\
aln_only=${aln_only} \
${pipeline_dir}/04_markduplicates.sambamba.markdup.sh
    # log to main
    echo "03: bam file(s) for sample ${sample} have been merged." | tee -a main.log
    # mv logfile
    mv ${sample}.sambamba-merge.log all_logfiles
fi
