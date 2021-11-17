#!/bin/bash
#PBS -l nodes=1:ppn=8,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.analyses.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load gatk/4.0.1.2
module load samtools/1.10
module load bcftools/1.11
module load R/4.1.0
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# create output dirs
if [[ ! -e analyses ]]; then
    mkdir analyses
fi
# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi


# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# estimate coverage
coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 ${dir}/${tumor}.bqsr.bam ${dir}/${normal}.bqsr.bam | awk '$3 >= 4 && $4 >= 4' | wc -l)
# expected coverage
expected=$(cat $intervals_bed | awk '{ count = count + ($3 - ($2 + 1)) } END { print count }')

# estimate tumor mutation burden (TMB)
# use prev coverage estimate

# total snvs
total_snvs=$(bcftools view --types snps vcf/${tumor}__${normal}.mutect2.annotated.${mode}.vcf.gz | grep -v "^#" | wc -l)
# total indels
total_indels=$(bcftools view --types indels vcf/${tumor}__${normal}.mutect2.annotated.${mode}.vcf.gz | grep -v "^#" | wc -l)
# calc TMB
TMB_snvs=$( echo "scale=4; ${total_snvs}/(${coverage}/1000000)" | bc )
TMB_indels=$( echo "scale=4; ${total_indels}/(${coverage}/1000000)" | bc )

# output
echo "${tumor},${normal},${coverage},${expected},${total_snvs},${total_indels},${TMB_snvs},${TMB_indels}" >> analyses/coverage_and_tmb.csv
echo "tumor mutation burden done"

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # final logs and tidy up dir, ie. gargabe collection
    # run final COSMIC signature analysis
    echo "${tumor},${normal}" >> finished.csv
    finished=$( cat finished.csv | wc -l )
    started=$( cat tumors_and_normals.csv | grep -v "^#" | wc -l )
    if [[ "$finished" -eq "$started" ]]; then
      # add example of how to load signature data to R
      echo -e "
# library path to standard and required additional libraries
.libPaths('/hpf/largeprojects/tabori/software/R_libs/4.1.0/')

# load libraries
library(sigminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)

# load results from analysis
load(\"analyses/mutational_signatures_as_R_object.Rdata\")

# do stuff with data...
# you should find:
# (1) linear_decomp_mt_sig_legacy_30 (matrix with matched signatures, COSMIC v2)
# (2) linear_decomp_mt_sig_sbs_96 (matrix with matched signatures, COSMIC v3)
# (3) maf (vcf files inmaf format)
# (4) matched_mt_sig_legacy_30 (data frame with matched mutational signatures, COSMIC v2)
# (5) matched_mt_sig_sbs_96 (data frame with matched mutational signatures, COSMIC v3)
# (6) mt_sig_bayes_sbs_96 (Bayesian NMF fitted signatures)
# (7) mt_tally (tally of SBS, DBS, and indel mutational patterns)
# (8) tmb_data (tumor mutational burden data from VCFs)

" > analyses/revisit_signature_data.R
        # add header to TMB csv
        sed -i '1 s/^/tumor,normal,obs_coverage,exp_coverage,snvs,indels,tmb_snvs,tmb_indels\n/' analyses/coverage_and_tmb.csv
        # run mutational signature analysis
        Rscript ${pipeline_dir}/cosmic_signature_analysis.R ${mode}
        # tidy and rename
        # rename dirs
        if [[ ! -e bam ]]; then
            mv BQSR bam
        else
            mv BQSR/* bam
            rm -rf BQSR
        fi
        # move stats to all_logfiles
        mv mutect2/*.mutect2.filtered.wes.vcf.filteringStats.tsv all_logfiles
        mv mutect2/*.mutect2.unfiltered.wes.merged.vcf.stats all_logfiles
        # delete all other vcf files except unfiltered VCFs
        rm $(ls mutect2/*vcf* | grep -v "unfiltered")
        # bgzip and tabix all vcf files
        ls mutect2/*.vcf | parallel --tmpdir ./tmp "bgzip {} && tabix {}.gz"
        # delete directories with bam data
        rm -rf preprocessed_bam aligned_bam tmp
        if [[ -e unmapped_bam ]]; then
            rm -rf unmapped_bam
        fi
        # final log
        echo "pipeline finished." | tee -a main.log
        # log final
        date | tee -a main.log
        # get date pipeline started
        start_date=$(head -1 main.log)
        end_date=$(tail -1 main.log)
        # calculate total running time
        sds=$(date -d "$start_date" +%s)
        eds=$(date -d "$end_date" +%s)
        total_time_in_days=$( echo "scale=5; ($eds - $sds) / 86400" | bc)
        # add 0 if less than 1
        if [[ $(echo "${total_time_in_days} > 1" | bc) == 0 ]]; then
          total_time_in_days="0${total_time_in_days}"
        fi
        # final log
        echo -e "\npipeline took ${total_time_in_days} days to complete" | tee -a main.log
        # change/adjust permisions
        # this configuration allows the main user and the users in the tabori group to
        # read/write/excecute
        # dirs first
        chmod 774 all_logfiles analyses bam contamination mutect2/f1r2 vcf/snpEff
        # files second
        chmod 664 all_logfiles/* analyses/* bam/* contamination/* mutect2/*vcf* mutect2/f1r2/* vcf/*vcf* vcf/snpEff/* ${tumor}__${normal}.analyses.log
    fi
    # last log and move logfile to dir
    echo "Done for ${tumor}__${normal}"
    mv ${tumor}__${normal}.analyses.log all_logfiles
fi
