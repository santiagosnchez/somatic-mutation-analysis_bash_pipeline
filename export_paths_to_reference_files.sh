#!/bin/bash
# pipeline folder
export pipeline_dir=/hpf/largeprojects/tabori/santiago/pipeline
# path to human reference genome assembly v38
export reference=/hpf/largeprojects/tabori/reference/gatk_bundle/hg38/Homo_sapiens_assembly38.fasta
# path to WES target intervals
export intervals=/hpf/largeprojects/tabori/reference/dap/AgilentSureSelectV5_LiftOver-hg19TOhg38_CoveredRegions.interval_list
# path to WES tergets in bed format
export intervals_bed=/hpf/largeprojects/tabori/reference/dap/AgilentSureSelectV5_LiftOver-hg19TOhg38_CoveredRegions.bed
# path to vcf file with known SNPs from the 1000 genomes project
export knownsites=/hpf/largeprojects/tabori/reference/gatk_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf
# path to WES intervals for running MuTect2
export bed30intervals=/hpf/largeprojects/tabori/reference/dap/AgilentSureSelectV5-exome-interval-files-gatk-split-30/
# path to gnomad resource
export gnomad_resource=/hpf/largeprojects/tabori/gnomad_resources/gnomad.biallelic.AF0.0001.PASS.vcf
# path to snpEff jar file
export snpeff_jar=/hpf/tools/centos6/snpEff/4.11/snpEff.jar
# path to snpEff data dir
export snpeff_datadir=/hpf/largeprojects/tabori/software/snpEff_data/4.11/data

# functions

# estimate walltime length
get_walltime(){
    size=$(du -sc $* | tail -1 | cut -f1)
    walltime=$(echo "scale=0; (${size} * 4)/10000000" | bc) 
    if [[ "${walltime}" == 0 ]]; then
        walltime=2
    fi
    echo $walltime
}
export -f get_walltime


