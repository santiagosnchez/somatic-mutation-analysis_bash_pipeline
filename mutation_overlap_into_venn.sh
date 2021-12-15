#!/bin/bash

# first make sure we have the right modules

module load bcftools/1.11
module load R/4.1.0
module load tabix

# get files

vcf_files=($@)

# initialize samples

declare -a samples=()

# check files and get tumor samples

for idx in `seq 1 ${#vcf_files[@]}`; do
    idx=$(( idx - 1 ))
    if [[ "${vcf_files[idx]}" != *".vcf.gz" ]]; then
        # get sample from name
        sample=$(echo ${vcf_files[idx]} | sed 's/__.*//')
        # compress and index
        bgzip "${vcf_files[idx]}"
        tabix "${vcf_files[idx]}".gz
        # add samples and rename vcf
        vcf_files[idx]="${vcf_files[idx]}".gz
        samples[idx]=${sample}
    else
        # get sample from name
        sample=$(echo ${vcf_files[idx]} | sed 's/__.*//')
        samples[idx]=${sample}
        if [[ ! -e "${vcf_files[idx]}.tbi" ]]; then
            tabix "${vcf_files[idx]}"
        fi
    fi
done

# load paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# output file
outfile=$( echo ${samples[@]} | sed 's/ /__/g' ).merged_vars.csv

# use bcftools to generate an overlap table
# Vars with with multiple alternate alleles are excluded

bcftools merge \
 --force-samples ${vcf_files} | \
 bcftools view -m 2 -s $( echo ${samples[@]} | sed 's/ /,/g' ) | \
 bcftools query -f "%CHROM,%POS,%REF,%ALT,[%GT,]\n" | sed 's/,$//' > $outfile

# check if finished
if [[ "$?" == 0 ]]; then
    # add header
    samples_joined=$( echo ${samples[@]} | sed 's/ /,/g' )
    sed -i "1 s/^/CHROM,POS,REF,ALT,$samples_joined\n/" $outfile
    # run R script
    Rscript ${pipeline_dir} $outfile
else
    echo "bcftools failed. Check files"
fi
