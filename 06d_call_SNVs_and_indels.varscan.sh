#!/bin/bash
#PBS -l nodes=1:ppn=10,vmem=30g,mem=30g,walltime=12:00:00
#PBS -e ${tumor}__${normal}.VarScan.log
#PBS -j eo
# scheduler settings

# load modules
module load parallel/20210322
module load java/1.8
module load varscan/2.3.8
module load bcftools/1.11
module load sambamba/0.7.0
module load samtools/1.10
module load tabix

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e varscan ]]; then
    mkdir -p varscan/pileups
fi
# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

if [[ ! -e varscan/${tumor}__${normal}.varscan.all.Somatic.hc.${mode}.vcf.gz ]]; then
# first run mpileup in parallel
if [[ -e varscan/pileups/${normal}.pileups ]]; then
    sambamba mpileup -L $intervals_bed -t 10 ${dir}/${tumor}.bqsr.bam > varscan/pileups/${tumor}.pileup
else
    export dir
    parallel "sambamba mpileup -L $intervals_bed -t 5 ${dir}/{}.bqsr.bam > varscan/pileups/{}.pileup" ::: ${tumor} ${normal}
fi
# run varscan
if [[ "$?" == 0 ]]; then
    java -jar ${varscan_jar} somatic \
     varscan/pileups/${normal}.pileup \
     varscan/pileups/${tumor}.pileup \
     varscan/${tumor}__${normal}.varscan \
     --output-vcf 1 \
     --strand-filter

    if [[ "$?" == 0 ]]; then
        # first SNVs
        java -jar ${varscan_jar} processSomatic \
          varscan/${tumor}__${normal}.varscan.snp.vcf \
          --min-tumor-freq 0.10 \
          --max-normal-freq 0.05 \
          --p-value 0.05
        # then indels
        java -jar ${varscan_jar} processSomatic \
          varscan/${tumor}__${normal}.varscan.indel.vcf \
          --min-tumor-freq 0.10 \
          --max-normal-freq 0.05 \
          --p-value 0.05
        # compress and index
        index-vcf varscan/${tumor}__${normal}.varscan.snp.Somatic.hc.vcf
        index-vcf varscan/${tumor}__${normal}.varscan.indel.Somatic.hc.vcf
        # concat
        bcftools concat \
         -a \
         -Oz \
         varscan/${tumor}__${normal}.varscan.snp.Somatic.hc.vcf.gz \
         varscan/${tumor}__${normal}.varscan.indel.Somatic.hc.vcf.gz \
         > varscan/${tumor}__${normal}.varscan.all.Somatic.hc.vcf.gz
        # rename
        mv varscan/${tumor}__${normal}.varscan.all.Somatic.hc.vcf.gz varscan/${tumor}__${normal}.varscan.all.Somatic.hc.${mode}.vcf.gz
        # index
        tabix varscan/${tumor}__${normal}.varscan.all.Somatic.hc.${mode}.vcf.gz
    fi
fi

else
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
    # log to main
    echo "${tumor}__${normal} VarScan2 variant calling completed." | tee -a main.log
    # move logfile
    mv ${tumor}__${normal}.VarScan.log all_logfiles
fi
