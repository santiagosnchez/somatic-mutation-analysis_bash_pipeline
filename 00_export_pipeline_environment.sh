#!/bin/bash

# add a bin directory and additional default locations
export PATH="/hpf/largeprojects/tabori/shared/bin:${PATH}"
export PERL5LIB="/hpf/largeprojects/tabori/shared/software/perl5/lib/perl5:${PERL5LIB}"
export PYTHONPATH="/hpf/largeprojects/tabori/shared/software/lib/python3.7:/hpf/largeprojects/tabori/shared/software/lib/python3.7/site-packages:/hpf/largeprojects/tabori/shared/software/lib/python3.7/dyn-lib:${PYTHONPATH}"

# base dirs
resources_dir=/hpf/largeprojects/tabori/shared/resources
software_dir=/hpf/largeprojects/tabori/shared/software
genomes=${resources_dir}/reference_genomes

# set hg38 as default
if [[ -z $1 ]]; then
    # path to human reference genome assembly v38
    export reference=${genomes}/hg38/gatk_bundle/Homo_sapiens_assembly38.fasta
    # path to reference dict
    export reference_dict=${genomes}/hg38/gatk_bundle/Homo_sapiens_assembly38.dict
    # path to WES target intervals
    export intervals=${genomes}/hg38/AgilentSureSelectV5/S04380110_Covered.hg38.interval_list
    # path to WES tergets in bed format
    export intervals_bed=${genomes}/hg38/AgilentSureSelectV5/S04380110_Covered.hg38.bed
    # path to vcf file with known SNPs from the 1000 genomes project
    export knownsites_snps=${genomes}/hg38/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    # same but biallelic only
    export knownsites_snps_biallelic=${genomes}/hg38/gatk_bundle/1000G_phase1.snps.high_confidence.biallelic.hg38.vcf.gz
    # path to vcf file with known indels from the 1000 genomes project
    export knownsites_indels=${genomes}/hg38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    # path to WES intervals for running MuTect2
    export bed30intervals=${genomes}/hg38/AgilentSureSelectV5/S04380110_Covered.hg38.30-bed-files/
    # path to gnomad resource
    export gnomad_resource=${genomes}/hg38/gatk_bundle/af-only-gnomad.hg38.vcf.gz
    # path to gatk's panel of normals vcf
    export gatk_pon=${genomes}/hg38/gatk_bundle/1000g_pon.hg38.vcf.gz
    # reference independent locations:
    # path to snpEff jar file
    export snpeff_jar=/hpf/tools/centos6/snpEff/4.11/snpEff.jar
    # path to snpEff data dir
    export snpeff_datadir=${resources_dir}/snpEff_data/4.11/data
    # path to varscan jar file
    export varscan_jar=/hpf/tools/centos6/varscan/2.3.8/VarScan.v2.3.8.jar
    # point to recent version of gatk
    export gatk_path=${software_dir}/gatk/gatk-4.2.3.0
    # funcotator data resources
    export funcotator_databases_s=${resources_dir}/funcotator_dataSources.v1.7.20200521s
    export funcotator_databases_g=${resources_dir}/funcotator_dataSources.v1.7.20200521g
    # path to annovar database
    export annovar_db=${resources_dir}/humandb
else
    # test organism
    if [[ ${1} == "human" ]]; then
        # test reference version
        if [[ ${2} == "hg38" ]]; then
            # path to human reference genome assembly v38
            export reference=${genomes}/${2}/gatk_bundle/Homo_sapiens_assembly38.fasta
            # path to reference dict
            export reference_dict=${genomes}/${2}/gatk_bundle/Homo_sapiens_assembly38.dict
            # path to vcf file with known SNPs from the 1000 genomes project
            export knownsites_snps=${genomes}/${2}/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
            # same but biallelic only
            export knownsites_snps_biallelic=${genomes}/hg38/gatk_bundle/1000G_phase1.snps.high_confidence.biallelic.hg38.vcf.gz
            # path to vcf file with known indels from the 1000 genomes project
            export knownsites_indels=${genomes}/${2}/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
            # path to gnomad resource
            export gnomad_resource=${genomes}/${2}/gatk_bundle/af-only-gnomad.hg38.vcf.gz
            # path to gatk's panel of normals vcf
            export gatk_pon=${genomes}/${2}/gatk_bundle/1000g_pon.hg38.vcf.gz
            # test mode
            if [[ ${3} == "wes" ]]; then
                # path to WES target intervals
                export intervals=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.hg38.interval_list
                # path to WES tergets in bed format
                export intervals_bed=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.hg38.bed
                # path to WES intervals for running MuTect2
                export bed30intervals=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.hg38.30-bed-files/
            else
                # path to WES target intervals
                export intervals=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.interval_list
                # path to WES tergets in bed format
                export intervals_bed=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.bed
                # path to WES intervals for running MuTect2
                export bed30intervals=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.30-bed-files/
            fi
            # reference independent locations:
            # path to snpEff jar file
            export snpeff_jar=/hpf/tools/centos6/snpEff/4.11/snpEff.jar
            # path to snpEff data dir
            export snpeff_datadir=/hpf/largeprojects/tabori/shared/resources/snpEff_data/4.11/data
            # path to varscan jar file
            export varscan_jar=/hpf/tools/centos6/varscan/2.3.8/VarScan.v2.3.8.jar
            # point to recent version of gatk
            export gatk_path=/hpf/largeprojects/tabori/shared/software/gatk/gatk-4.2.3.0
            # funcotator data resources
            export funcotator_databases_s=${resources_dir}/funcotator_dataSources.v1.7.20200521s
            export funcotator_databases_g=${resources_dir}/funcotator_dataSources.v1.7.20200521g
            # path to annovar database
            export annovar_db=${resources_dir}/humandb
        elif [[ ${2} == "hs37d5" ]]; then
            # path to human reference genome assembly hs37d5
            export reference=${genomes}/${2}/gatk_bundle/Homo_sapiens_assembly38.fasta
            # path to vcf file with known SNPs from the 1000 genomes project
            export knownsites_snps=${genomes}/${2}/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
            # path to vcf file with known indels from the 1000 genomes project
            export knownsites_indels=${genomes}/${2}/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
            # path to gnomad resource
            export gnomad_resource=${genomes}/${2}/gatk_bundle/af-only-gnomad.hg38.vcf.gz
            # path to gatk's panel of normals vcf
            export gatk_pon=${genomes}/${2}/gatk_bundle/1000g_pon.hg38.vcf.gz
            # test mode
            if [[ ${3} == "wes" ]]; then
                # path to WES target intervals
                export intervals=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.edited.LiftOverToHg38.interval_list
                # path to WES tergets in bed format
                export intervals_bed=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.edited.LiftOverToHg38.bed
                # path to WES intervals for running MuTect2
                export bed30intervals=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.edited.LiftOverToHg38.30-bed-files/
            else
                # path to WES target intervals
                export intervals=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.interval_list
                # path to WES tergets in bed format
                export intervals_bed=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.bed
                # path to WES intervals for running MuTect2
                export bed30intervals=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.30-bed-files/
            fi
        elif [[ ${2} == "b37" ]]; then
              # path to human reference genome assembly hs37d5
              export reference=${genomes}/${2}/gatk_bundle/Homo_sapiens_assembly38.fasta
              # path to vcf file with known SNPs from the 1000 genomes project
              export knownsites_snps=${genomes}/${2}/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
              # path to vcf file with known indels from the 1000 genomes project
              export knownsites_indels=${genomes}/${2}/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
              # path to gnomad resource
              export gnomad_resource=${genomes}/${2}/gatk_bundle/af-only-gnomad.hg38.vcf.gz
              # path to gatk's panel of normals vcf
              export gatk_pon=${genomes}/${2}/gatk_bundle/1000g_pon.hg38.vcf.gz
              # test mode
              if [[ ${3} == "wes" ]]; then
                  # path to WES target intervals
                  export intervals=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.edited.LiftOverToHg38.interval_list
                  # path to WES tergets in bed format
                  export intervals_bed=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.edited.LiftOverToHg38.bed
                  # path to WES intervals for running MuTect2
                  export bed30intervals=${genomes}/${2}/AgilentSureSelectV5/S04380110_Covered.edited.LiftOverToHg38.30-bed-files/
              else
                  # path to WES target intervals
                  export intervals=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.interval_list
                  # path to WES tergets in bed format
                  export intervals_bed=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.bed
                  # path to WES intervals for running MuTect2
                  export bed30intervals=${genomes}/${2}/gatk_bundle/wgs_calling_regions.hg38.30-bed-files/
              fi
              # reference independent locations:
              # path to snpEff jar file
              export snpeff_jar=/hpf/tools/centos6/snpEff/4.11/snpEff.jar
              # path to snpEff data dir
              export snpeff_datadir=${resources_dir}/snpEff_data/4.11/data
              # path to varscan jar file
              export varscan_jar=/hpf/tools/centos6/varscan/2.3.8/VarScan.v2.3.8.jar
              # point to recent version of gatk
              export gatk_path=${software_dir}/gatk/gatk-4.2.3.0
              # funcotator data resources
              export funcotator_databases_s=${resources_dir}/funcotator_dataSources.v1.7.20200521s
              export funcotator_databases_g=${resources_dir}/funcotator_dataSources.v1.7.20200521g
              # path to annovar database
              export annovar_db=${resources_dir}/humandb
        # more refs
        else
            echo "Unknown or reference no databased: ${2}"
            echo "See options using the -h flag."
            return 1
        fi
    elif [[ ${1} == "mouse" ]]; then
        if [[ ${2} == "mm10" ]]; then # need editing
            # path to mouse reference genome assembly mm10 ... with "chr"
            export reference=${genomes}/${2}/ucsc/Mus_musculus.mm10.fa
            # path to vcf file with known SNPs from dbsnp v150
            export knownsites_snps=${genomes}/${2}/ucsc/dbsnp_v150.snps-biallelic.mm10.vcf.gz
            # path to vcf file with known indels from the 1000 genomes project
            export knownsites_indels=${genomes}/${2}/ucsc/mgp.v5.indels.pass.mm10.vcf.gz
            # path to gnomad resource
            export gnomad_resource=null
            # path to gatk's panel of normals vcf
            export gatk_pon=${genomes}/${2}/ucsc/pon.sorted.mm10.vcf.gz
            # test mode
            if [[ ${3} == "wes" ]]; then
                # path to WES target intervals
                export intervals=${genomes}/${2}/AgilentSureSelectV1/S0276129_Covered.sorted_noMT.LiftOverToMm10.interval_list
                # path to WES tergets in bed format
                export intervals_bed=${genomes}/${2}/AgilentSureSelectV1/S0276129_Covered.sorted_noMT.LiftOverToMm10.bed
                # path to WES intervals for running MuTect2
                export bed30intervals=${genomes}/${2}/AgilentSureSelectV1/S0276129_Covered.sorted_noMT.LiftOverToMm10.30-bed-files/
            else
                # path to WES target intervals
                export intervals=${genomes}/${2}/ucsc/wgs_callable_nonN-LC_regions.interval_list
                # path to WES tergets in bed format
                export intervals_bed=${genomes}/${2}/ucsc/wgs_callable_nonN-LC_regions.bed
                # path to WES intervals for running MuTect2
                export bed30intervals=${genomes}/${2}/ucsc/wgs_callable_nonN-LC_regions.30-bed-files/
            fi
            # reference independent locations:
            # path to snpEff jar file
            export snpeff_jar=/hpf/tools/centos6/snpEff/4.11/snpEff.jar
            # path to snpEff data dir
            export snpeff_datadir=${resources_dir}/snpEff_data/4.11/data
            # path to varscan jar file
            export varscan_jar=/hpf/tools/centos6/varscan/2.3.8/VarScan.v2.3.8.jar
            # point to recent version of gatk
            export gatk_path=${software_dir}/gatk/gatk-4.2.3.0
            # funcotator data resources
            #export funcotator_databases_s=${resources_dir}/funcotator_dataSources.v1.7.20200521s
            #export funcotator_databases_g=${resources_dir}/funcotator_dataSources.v1.7.20200521g
            # path to annovar database
            #export annovar_db=${resources_dir}/mousedb
        else
            echo "Unknown or reference no databased: ${2}"
            echo "See options using the -h flag."
            return 1
        fi
    else
        echo "Unknown or organism not databased: ${1}"
        echo "See options using the -h flag."
        return 1
    fi
fi



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
# 1st arg: date
# 2nd arg: h|d  (in fraction of days or hours)
how_long(){
  if [[ -f ${1} ]]; then
    start_date=$(head -1 $1)
  else
    start_date="$1"
  fi
  end_date=$(date)
  # calculate total running time
  sds=$(date -d "$start_date" +%s)
  eds=$(date -d "$end_date" +%s)
  # in days
  if [[ "${2}" == "d" ]]; then
    # divided by seconds in a day
    total_time=$( echo "scale=5; ($eds - $sds) / 86400" | bc )
  elif [[ "${2}" == "h" ]]; then
    # divided by seconds in an hour
    total_time=$( echo "scale=5; ($eds - $sds) / 3600" | bc )
  fi
  # add 0 if less than 1
  if [[ $(echo "${total_time} > 1" | bc) == 0 ]]; then
    total_time="0${total_time}"
  fi
  echo $total_time
}
export -f how_long
