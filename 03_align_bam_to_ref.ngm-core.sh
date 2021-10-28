#PBS -l nodes=1:ppn=20,vmem=30g,mem=30g
#PBS -e ${sample}.ngm-core.log
#PBS -j eo
# scheduler settings

# load modules
module load opencl/18.1.0.015
module load nextgenmap/0.5.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load all paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# check if prev step finished correctly
if [[ $(samtools quickcheck -u unmapped_bam/${sample}.unaligned.merged.bam && echo 1) != 1 ]]; then
     echo "resubmitting previous step and increase time by 2hrs"
     wt=$(( wt + 2 ))
     qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt} ${pipeline_dir}/02_fastq_to_ubam.picard.Fastq2Sam.samtools.merge.sh
     exit 0
fi

# create dir for ubam
if [[ ! -e aligned_bam ]]; then
    mkdir aligned_bam
fi

# check if input files exist and run aligner
if [[ ! -e aligned_bam/${sample}.bam ]]; then 
ngm-core -r $reference \
 -q unmapped_bam/${sample}.unaligned.merged.bam \
 -o aligned_bam/${sample}.bam \
 --bam \
 --keep-tags \
 --no-unal \
 -t 20 \
 --very-sensitive \
 --min-identity 0.7 \
 --min-residues 0.5 
else
# check integrity of bam file
samtools quickcheck aligned_bam/${sample}.bam
    if [[ "$?" != 0 ]]; then
ngm-core -r $reference \
 -q unmapped_bam/${sample}.unaligned.merged.bam \
 -o aligned_bam/${sample}.bam \
 --bam \
 --keep-tags \
 --no-unal \
 -t 20 \
 --very-sensitive \
 --min-identity 0.7 \
 --min-residues 0.5 
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
    # add read group to aligned bam
    Nreadgroups=$(samtools view -H unmapped_bam/${sample}.unaligned.merged.bam | grep -c "^@RG")
    for i in $( seq 1 $Nreadgroups ); do
        rg=$(samtools view -H unmapped_bam/${sample}.unaligned.merged.bam | grep "^@RG" | sed -n "${i}p")
        samtools addreplacerg -m orphan_only -r "$(echo -n $rg | tr ' ' '\t')" -o aligned_bam/${sample}.rg.bam aligned_bam/${sample}.bam
        mv aligned_bam/${sample}.rg.bam aligned_bam/${sample}.bam
    done
    echo "added read groups"
    # move logfiles
    mv ${sample}.ngm-core.log all_logfiles
    # submit next job
    # can switch this to picards MarkDuplicate method
    qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/04a_markduplicates.samtools.markdup.sh
fi

