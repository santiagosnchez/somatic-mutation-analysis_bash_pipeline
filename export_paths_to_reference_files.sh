#!/bin/bash
# pipeline folder
export pipeline_dir=/hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery
# path to human reference genome assembly v38
export reference=/hpf/largeprojects/tabori/shared/resources/hg38/gatk_bundle/Homo_sapiens_assembly38.fasta
# path to WES target intervals
export intervals=/hpf/largeprojects/tabori/shared/resources/hg38/AgilentSureSelectV5/AgilentSureSelectV5_LiftOver-hg19TOhg38_CoveredRegions.interval_list
# path to WES tergets in bed format
export intervals_bed=/hpf/largeprojects/tabori/shared/resources/hg38/AgilentSureSelectV5/AgilentSureSelectV5_LiftOver-hg19TOhg38_CoveredRegions.bed
# path to vcf file with known SNPs from the 1000 genomes project
export knownsites_snps=/hpf/largeprojects/tabori/shared/resources/hg38/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
# path to vcf file with known indels from the 1000 genomes project
export knownsites_indels=/hpf/largeprojects/tabori/shared/resources/hg38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
# path to WES intervals for running MuTect2
export bed30intervals=/hpf/largeprojects/tabori/shared/resources/hg38/AgilentSureSelectV5/AgilentSureSelectV5-exome-interval-files-gatk-split-30/
# path to gnomad resource
export gnomad_resource=/hpf/largeprojects/tabori/shared/resources/hg38/gnomad/gnomad.biallelic.AF0.0001.PASS.vcf
# path to snpEff jar file
export snpeff_jar=/hpf/tools/centos6/snpEff/4.11/snpEff.jar
# path to snpEff data dir
export snpeff_datadir=/hpf/largeprojects/tabori/shared/resources/snpEff_data/4.11/data
# path to varscan jar file
export varscan_jar=/hpf/tools/centos6/varscan/2.3.8/VarScan.v2.3.8.jar
# point to recent version of gatk
export gatk_path=/hpf/largeprojects/tabori/shared/software/gatk/gatk-4.2.3.0
# funcotator data resources
export funcotator_databases_s=/hpf/largeprojects/tabori/shared/resources/funcotator_dataSources.v1.7.20200521s
export funcotator_databases_g=/hpf/largeprojects/tabori/shared/resources/funcotator_dataSources.v1.7.20200521g

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

# get read groups from illumina header
get_read_group_info(){
  # get the first line
  file $1 | grep "gzip" &> /dev/null
  if [[ "$?" == 0 ]]; then # gzipped
  # get the first line
    head=$(zcat $1 2> /dev/null | head -1 | sed 's/^@//')
  else
    head=$(cat $1 2> /dev/null | head -1 | sed 's/^@//')
  fi
  head_split=(`echo $head | tr ':' '\n'`)
  # default assume illumina
  PL=ILLUMINA
  # sample second arg
  SM=$2
  if [[ "${#head_split[@]}" == 11 ]]; then
    PM=${head_split[0]} # instrument id
    ID="${head_split[1]}-${head_split[3]}" # merge run id with lane id
    PU=${head_split[2]} # flowcell id
    BC=${head_split[10]} # barcode ID
    # build read group string
    #RG="@RG\\\tID:${ID}\\\tSM:${SM}\\\tLB:${BC}\\\tPL:${PL}\\\tBC:${BC}\\\tPU:${PU}\\\tPM:${PM}"
    RG="@RG\\tID:${ID}\\tSM:${SM}\\tLB:${BC}\\tPL:${PL}\\tBC:${BC}\\tPU:${PU}\\tPM:${PM}"
  else
    ID=1 # run id
    #RG="@RG\\\tID:${ID}\\\tSM:${SM}\\\tPL:${PL}"
    RG="@RG\\tID:${ID}\\tSM:${SM}\\tPL:${PL}"
  fi
  echo "$RG"
}
export -f get_read_group_info

# compresses and generates a tabix index
index-vcf(){
    if [[ -e $1.gz ]]; then
        rm $1.gz $1.gz.tbi
    fi
    bgzip $1 && tabix $1.gz
}
export -f index-vcf

# function to look for file
file_lookup(){
    until [[ -e $1 ]]; do
        # checks for the file every minute
        sleep 60
    done
    echo "file $1 found."
    return 0
}
export -f file_lookup

# calculate how long it took to run
how_long(){
  start_date=$(head -1 $1)
  end_date=$(date)
  # calculate total running time
  sds=$(date -d "$start_date" +%s)
  eds=$(date -d "$end_date" +%s)
  total_time_in_days=$( echo "scale=5; ($eds - $sds) / 86400" | bc)
  # add 0 if less than 1
  if [[ $(echo "${total_time_in_days} > 1" | bc) == 0 ]]; then
    total_time_in_days="0${total_time_in_days}"
  fi
  echo $total_time_in_days
}
