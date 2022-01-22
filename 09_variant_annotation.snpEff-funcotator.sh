#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.snpEff.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load snpEff/4.11
module load bcftools/1.11
module load tabix

# set working dir
cd $PBS_O_WORKDIR

# create output dirs
if [[ ! -e vcf/snpEff ]]; then
    mkdir -p vcf/snpEff
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# make temporary header file with pedigree
echo "##PEDIGREE=<Derived=${tumor},Original=${normal}>" > ${tumor}__${normal}.tmp.vcf.header.txt

# normalize variants, reduce complex alleles and add header line
bcftools annotate \
 -h ${tumor}__${normal}.tmp.vcf.header.txt \
 -Oz -o mutect2/${tumor}__${normal}.mutect2.normalized_head.${mode}.vcf.gz \
 mutect2/${tumor}__${normal}.mutect2.normalized.${mode}.vcf.gz

# delete tmp header file
if [[ "$?" == 0 ]]; then
    rm ${tumor}__${normal}.tmp.vcf.header.txt
fi

# run snpEff
java -jar $snpeff_jar \
 -dataDir $snpeff_datadir \
 hg38 \
 -v \
 -cancer \
 -stats vcf/snpEff/${tumor}__${normal}.snpEff_summary.html \
 -csvStats vcf/snpEff/${tumor}__${normal}.snpEff_summary.csv \
 mutect2/${tumor}__${normal}.mutect2.normalized_head.${mode}.vcf.gz > \
 vcf/${tumor}__${normal}.mutect2.annotated-snpeff.${mode}.vcf

# run gatk's funcotator
$gatk_path/gatk Funcotator \
 --variant mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz \
 --reference $reference \
 --ref-version hg38 \
 --data-sources-path $funcotator_databases \
 --transcript-selection-mode CANONICAL \
 --output vcf/${tumor}__${normal}.mutect2.annotated-funcotator.${mode}.vcf \
 --output-file-format VCF

# same but output maf
 $gatk_path/gatk Funcotator \
  --variant mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz \
  --reference $reference \
  --ref-version hg38 \
  --data-sources-path $funcotator_databases \
  --transcript-selection-mode CANONICAL \
  --output vcf/${tumor}__${normal}.mutect2.annotated-funcotator.${mode}.maf \
  --output-file-format MAF

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # bgzip and tabix vcf
    ls vcf/${tumor}__${normal}.mutect2.annotated* | parallel index-vcf {}
    # log to main
    echo "Annotation with SnpEff completed for ${tumor}__${normal}." | tee -a main.log
    # run analyses
    qsub -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/10_run_analyses.signatures_and_TBM.sh
     # prepare for cleanup
#    echo "${tumor},${normal}" >> finished.csv
#    finished=$( cat finished.csv | wc -l )
#    started=$( cat tumors_and_normals.csv | wc -l )
#    if [[ "$finished" -eq "$started" ]]; then
#        # rename dirs
#        mv BQSR bam
#        # move stats to all_logfiles
#        mutect2/${tumor}__${normal}.mutect2.filtered.wes.vcf.filteringStats.tsv all_logfiles
#        mutect2/${tumor}__${normal}.mutect2.unfiltered.wes.merged.vcf.stats all_logfiles
#        # delete all other vcf files
#        rm -rf mutect2/${tumor}__${normal}*.vcf*
#        # delete directories with bam data
#        rm -rf preprocessed_bam aligned_bam
#        if [[ -e unmapped_bam ]]; then
#            rm -rf unmapped_bam
#        fi
#    fi
    # move logfile
    mv ${tumor}__${normal}.snpEff.log all_logfiles
fi
