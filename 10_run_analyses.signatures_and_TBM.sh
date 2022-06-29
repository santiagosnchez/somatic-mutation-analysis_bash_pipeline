#!/bin/bash
#PBS -l nodes=1:ppn=8,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.analyses.${tissue}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
#module load gatk/4.0.1.2
module load samtools/1.10
module load bcftools/1.11
#module load R/4.1.0
module load R/4.0.2
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e analyses ]]; then
    mkdir analyses
fi

# create tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi

# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

# list of MMR genes:
# MLH1 MSH2 MSH6 PMS2 POLD1 POLE IDH1 TP53 NF1
# grep -w "MLH1\|MSH2\|MSH6\|PMS2\|POLD1\|POLE\|IDH1\|TP53\|NF1"

# pull germline and somatic missense (nonsynonymous) mutations
# look for MMR genes
# germline on VarScan calls

# ${pipeline_dir}/get_gene_annotations_from_vcf-funcotator.sh \
#  vcf/${tumor}__${normal}.varscan.all.Germline.annotated-funcotator.${mode}.vcf.gz \
#  MLH1 \
#  MSH2 \
#  MSH6 \
#  PMS2 \
#  POLD1 \
#  POLE \
#  IDH1 \
#  TP53 \
#  NF1 \
#  > analyses/${tumor}__${normal}.germline_MMR_mutations.genes.csv
# add IDH4

# ${pipeline_dir}/get_gene_annotations_from_vcf.sh \
# vcf/${tumor}__${normal}.varscan.all.Germline.annotated-snpeff.${mode}.vcf.gz \
# POLD1 \
# POLD2 \
# POLD3 \
# POLD4 \
# POLE \
# POLE2 > analyses/${tumor}__${normal}.germline_POL_mutations.genes.csv
#
# # somatic on Mutect2
# ${pipeline_dir}/get_gene_annotations_from_vcf-funcotator.sh \
#   vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz \
#   MLH1 \
#   MSH2 \
#   MSH6 \
#   PMS2 \
#   POLD1 \
#   POLE \
#   IDH1 \
#   TP53 \
#   NF1 \
#   > analyses/${tumor}__${normal}.somatic_MMR_mutations.genes.csv
#
# ${pipeline_dir}/get_gene_annotations_from_vcf.sh \
# vcf/${tumor}__${normal}.mutect2.annotated-snpeff.${mode}.vcf.gz \
# POLD1 \
# POLD2 \
# POLD3 \
# POLD4 \
# POLE \
# POLE2 > analyses/${tumor}__${normal}.somatic_POL_mutations.genes.csv

echo "10: Fetching all variant annotations." | tee -a main.log

# get all annotations into csv
# funcotator Somatic
if [[ ${tissue} == "Somatic" ]]; then
    if [[ -e vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz ]]; then
        ${pipeline_dir}/scripts/funcotator-vcf2maf.sh \
        vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz \
        ${tumor} ${normal} Somatic \
        > analyses/${tumor}__${normal}.annotations_funcotator_somatic.csv
    fi
    if [[ -e vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-snpeff.${mode}.vcf.gz ]]; then
        # get all annotations into csv
        # snpeff Somatic
        ${pipeline_dir}/scripts/snpeff-vcf2tbl.sh \
        vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-snpeff.${mode}.vcf.gz \
        ${tumor} ${normal} Somatic \
        > analyses/${tumor}__${normal}.annotations_snpeff_somatic.csv
    fi
    echo "10: done extracting somatic variants from vcf file for ${tumor}__${normal}" | tee -a main.log
fi

# check if tumor-only mode or germline
if [[ "${normal}" != "PON" && ${tissue} == "Germline" ]]; then
    if [[ -e vcf/${tumor}__${normal}.varscan.all.Germline.annotated-funcotator.${mode}.vcf.gz ]]; then
        ${pipeline_dir}/scripts/funcotator-vcf2maf.sh \
        vcf/${tumor}__${normal}.varscan.all.Germline.annotated-funcotator.${mode}.vcf.gz \
        ${tumor} ${normal} Germline \
        > analyses/${tumor}__${normal}.annotations_funcotator_germline.csv
    fi
    if [[ -e vcf/${tumor}__${normal}.varscan.all.Germline.annotated-snpeff.${mode}.vcf.gz ]]; then
        ${pipeline_dir}/scripts/snpeff-vcf2tbl.sh \
        vcf/${tumor}__${normal}.varscan.all.Germline.annotated-snpeff.${mode}.vcf.gz \
        ${tumor} ${normal} Germline \
        > analyses/${tumor}__${normal}.annotations_snpeff_germline.csv
    fi
    echo "10: done extracting germline variants from vcf file for ${tumor}__${normal}" | tee -a main.log
fi

if [[ ${tissue} == "Somatic" ]]; then
    # log
    echo "10: calculating observed coverage, SNVs and indels (${tumor}__${normal})" | tee -a main.log

    # add header to analyses/coverage_and_tmb.csv
    if [[ ! -e analyses/coverage_and_tmb.csv ]]; then
        echo "tumor,normal,obs_coverage,exp_coverage,snvs,snvs_nob,indels,indels_nob,tmb_snvs,tmb_snvs_nob,tmb_indels,tmb_indels_nob" > analyses/coverage_and_tmb.csv
    fi

    # expected coverage
    expected=$(cat $intervals_bed | awk '{ count = count + ($3 - ($2 + 1)) } END { print count }')
    # find observed coverage
    if [[ ! -e ${dir}/${tumor}.bqsr.bam ]]; then
        # if no bam, just use expected
        coverage=${expected}
    else
        if [[ "${normal}" != "PON" ]]; then
            coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 ${dir}/${tumor}.bqsr.bam ${dir}/${normal}.bqsr.bam | awk '$3 >= 4 && $4 >= 4' | wc -l)
        else
            coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 ${dir}/${tumor}.bqsr.bam | awk '$3 >= 4' | wc -l)
        fi
    fi

    # estimate tumor mutation burden (TMB)
    # use prev coverage estimate
    # total snvs
    total_snvs=$(bcftools view -H --types snps mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz | wc -l)
    total_snvs_nob=$(bcftools view -H --types snps mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf.gz | wc -l)
    # total indels
    total_indels=$(bcftools view -H --types indels mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz | wc -l)
    total_indels_nob=$(bcftools view -H --types indels mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf.gz | wc -l)
    # calc TMB
    TMB_snvs=$( echo "scale=2; ${total_snvs}/(${coverage}/1000000)" | bc | sed 's/^\./0\./')
    TMB_indels=$( echo "scale=2; ${total_indels}/(${coverage}/1000000)" | bc | sed 's/^\./0\./' )
    TMB_snvs_nob=$( echo "scale=2; ${total_snvs_nob}/(${coverage}/1000000)" | bc | sed 's/^\./0\./')
    TMB_indels_nob=$( echo "scale=2; ${total_indels_nob}/(${coverage}/1000000)" | bc | sed 's/^\./0\./' )

    # output
    echo "${tumor},${normal},${coverage},${expected},${total_snvs},${total_snvs_nob},${total_indels},${total_indels_nob},${TMB_snvs},${TMB_snvs_nob},${TMB_indels},${TMB_indels_nob}" >> analyses/coverage_and_tmb.csv
    echo "10: Done calculating tumor mutation burden (TMB)." | tee -a main.log

    # look at differences in calls between varscan, mutect2 with ob-priors and without.
    # varscan_snvs=$(bcftools view -H -v snps varscan/${tumor}__${normal}.all.Somatic.hc.vcf.gz | wc -l)
    # mutect2_all_filters_snvs=$(bcftools view -H -v snps -f PASS mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf | wc -l)
    # mutect2_all_filters_snvs=$(bcftools view -H -v snps -f PASS mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf | wc -l)

    # log
    echo "10: Running signature analysis (${tumor}__${normal})" | tee -a main.log

    # run variant analysis
    Rscript ${pipeline_dir}/scripts/variant_analysis.nofigs2.R ${mode} ${tumor}__${normal} ${organism}
    # for no-ob
    mkdir analyses/no-obpriors
    Rscript ${pipeline_dir}/scripts/variant_analysis.nofigs2.R ${mode} ${tumor}__${normal} ${organism} "no-obpriors"
fi

# add to archive
# zip -ru ${tumor}__${normal}.analyses.zip analyses/${tumor}__${normal}.*

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    if [[ ${mode} == "wgs" || ${normal} == "PON" || (-e all_logfiles/${tumor}__${normal}.analyses.Germline.log && ${tissue} == "Somatic") || (-e all_logfiles/${tumor}__${normal}.analyses.Somatic.log && ${tissue} == "Germline") ]]; then
        # final logs and tidy up dir, ie. gargabe collection
        # run final COSMIC signature analysis
        echo "${tumor},${normal}" >> finished.csv
        # get time it took
        total_time_in_days=$( how_long main.log )
        echo "pipeline finished for ${tumor}__${normal} in ${total_time_in_days} days" | tee -a main.log
        # check if all samples finished
        finished=$( cat finished.csv | sort -u | wc -l )
        started=$( cat tumors_and_normals.csv | grep -v "^#" | sort -u | wc -l )
        if [[ "$finished" -eq "$started" ]]; then
            # log and fetch MMR genes in annotations
            if [[ ${normal} == "PON" ]]; then
                echo "10: Fetching somatic variants of interest (${tumor}__${normal})" | tee -a main.log
            else
                echo "10: Fetching germline and somatic variants of interest (${tumor}__${normal})" | tee -a main.log
            fi

            # function to extrac rows that match MMR/tumor genes
            fetch_mmr_ann(){
              skip_rows=2
              col=10
              if [[ $(echo $1 | grep "funcotator" &> /dev/null && echo 1) == 1 ]]; then
                skip_rows=3
                col=4
              fi
              awk -v FS="," -v SR=$skip_rows -v COL=${col} \
              '{
                if (NR >= SR){
                  if ($COL ~ /^(MLH1|MSH2|MSH6|PMS2|POLD1|POLE|IDH1|TP53|NF1)$/){
                    print $0
                  }
                } else {
                 print $0
               }
              }' $1
              }

            # fetch_mmr_ann analyses/all_annotations_snpeff_somatic.csv > analyses/mmr_annotations_snpeff_somatic.csv
            # fetch_mmr_ann analyses/all_annotations_funcotator_somatic.csv > analyses/mmr_annotations_funcotator_somatic.csv
            # if [[ "${normal}" != "PON" ]]; then
            #     fetch_mmr_ann analyses/all_annotations_snpeff_germline.csv > analyses/mmr_annotations_snpeff_germline.csv
            #     fetch_mmr_ann analyses/all_annotations_funcotator_germline.csv > analyses/mmr_annotations_funcotator_germline.csv
            # fi

            # export to zip file
            today=$(date -I)
            zip -r export_results.${today}.zip annovar/*_multianno.maf analyses/old_output*

            # add tmb_and_coverage to archive
            #zip -ru analyses.zip analyses/coverage_and_tmb.csv
            # tidyup and clean working dir
            if [[ ! -e bam ]]; then
                mv BQSR bam
            else
                mv BQSR/* bam
                rm -rf BQSR
            fi
            #
            if [[ -e varscan/pileups ]]; then
              rm -rf varscan/pileups
            fi
            # move stats to all_logfiles
            mv mutect2/*.filteringStats.tsv all_logfiles
            #mv mutect2/*.mutect2.unfiltered.wes.merged.vcf.stats all_logfiles
            # delete all other vcf files except filtered and unfiltered VCFs
            #rm $(ls mutect2/*vcf* | grep -v "filtered\|selected")
            # bgzip and tabix all vcf files
            ls mutect2/*.vcf | parallel --tmpdir ./.tmp "index-vcf {}"
            # delete directories with bam data
            rm -rf preprocessed_bam aligned_bam .tmp
            if [[ -e unmapped_bam ]]; then
                rm -rf unmapped_bam
            fi
            # write final logs and change permissions
            grep "pipeline took" main.log
            if [[ "$?" != 0 ]]; then
                # final log
                echo "10: pipeline finished for batch." | tee -a main.log
                # log final
                date | tee -a main.log
                # get date pipeline started
                total_time_in_days=$( how_long main.log )
                # final log
                echo -e "\n10: pipeline took ${total_time_in_days} days for batch to complete" | tee -a main.log
                # change/adjust permisions
                # this configuration allows the main user and the users in the tabori group to
                # read/write/excecute
                # dirs first
                echo "10: changing permissions" | tee -a main.log
                find . -type d -user `whoami` -exec chmod 775 {} \;
                find . -type f -user `whoami` -exec chmod 664 {} \;
                echo "@@@@ All done. See you next time! @@@@"
            fi
        fi
    fi
    # run mutational signature analysis
    # Rscript ${pipeline_dir}/cosmic_signature_analysis.R ${mode} ${tumor}__${normal}.

    # last move logfile to dir
    if [[ -e ${tumor}__${normal}.analyses.${tissue}.log ]]; then
        echo "10: Done for ${tumor}__${normal} in ${tissue}" | tee -a main.log
        mv ${tumor}__${normal}.analyses.${tissue}.log all_logfiles
    fi
fi
