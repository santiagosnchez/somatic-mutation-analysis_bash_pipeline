#!/bin/bash
#PBS -l nodes=1:ppn=8,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.analyses.log
#PBS -j eo
# scheduler settings

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

# create output dirs
if [[ ! -e analyses ]]; then
    mkdir analyses
fi

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
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

# log
echo "Fetching germline and somatic variants of interest"

# pull germline and somatic missense (nonsynonymous) mutations
# look for MMR genes
# germline on VarScan calls

${pipeline_dir}/get_gene_annotations_from_vcf-funcotator.sh \
 vcf/${tumor}__${normal}.varscan.all.Germline.annotated-funcotator.${mode}.vcf.gz \
 MLH1 \
 MSH2 \
 MSH6 \
 PMS2 \
 POLD1 \
 POLE \
 IDH1 \
 TP53 \
 NF1 \
 > analyses/${tumor}__${normal}.germline_MMR_mutations.genes.csv
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
${pipeline_dir}/get_gene_annotations_from_vcf-funcotator.sh \
  vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz \
  MLH1 \
  MSH2 \
  MSH6 \
  PMS2 \
  POLD1 \
  POLE \
  IDH1 \
  TP53 \
  NF1 \
  > analyses/${tumor}__${normal}.somatic_MMR_mutations.genes.csv
#
# ${pipeline_dir}/get_gene_annotations_from_vcf.sh \
# vcf/${tumor}__${normal}.mutect2.annotated-snpeff.${mode}.vcf.gz \
# POLD1 \
# POLD2 \
# POLD3 \
# POLD4 \
# POLE \
# POLE2 > analyses/${tumor}__${normal}.somatic_POL_mutations.genes.csv

# add header to analyses/coverage_and_tmb.csv
if [[ ! -e analyses/coverage_and_tmb.csv ]]; then
    echo "tumor,normal,obs_coverage,exp_coverage,snvs,indels,tmb_snvs,tmb_indels" > analyses/coverage_and_tmb.csv
fi

coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 ${dir}/${tumor}.bqsr.bam ${dir}/${normal}.bqsr.bam | awk '$3 >= 4 && $4 >= 4' | wc -l)
# expected coverage
expected=$(cat $intervals_bed | awk '{ count = count + ($3 - ($2 + 1)) } END { print count }')

# estimate tumor mutation burden (TMB)
# use prev coverage estimate
# total snvs
total_snvs=$(bcftools view --types snps vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz | grep -v "^#" | wc -l)
# total indels
total_indels=$(bcftools view --types indels vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz | grep -v "^#" | wc -l)
# calc TMB
TMB_snvs=$( echo "scale=2; ${total_snvs}/(${coverage}/1000000)" | bc | sed 's/^\./0\./')
TMB_indels=$( echo "scale=2; ${total_indels}/(${coverage}/1000000)" | bc | sed 's/^\./0\./' )

# output
echo "${tumor},${normal},${coverage},${expected},${total_snvs},${total_indels},${TMB_snvs},${TMB_indels}" >> analyses/coverage_and_tmb.csv
echo "tumor mutation burden done"

# look at differences in calls between varscan, mutect2 with ob-priors and without.
# varscan_snvs=$(bcftools view -H -v snps varscan/${tumor}__${normal}.all.Somatic.hc.vcf.gz | wc -l)
# mutect2_all_filters_snvs=$(bcftools view -H -v snps -f PASS mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf | wc -l)
# mutect2_all_filters_snvs=$(bcftools view -H -v snps -f PASS mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf | wc -l)

# run variant analysis
Rscript ${pipeline_dir}/variant_analysis.R ${mode} ${tumor}__${normal}

# add to archive
zip -ru analyses.zip analyses/${tumor}__${normal}.*

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
    # get time it took
    total_time_in_days=$( how_long main.log )
    echo "pipeline finished for ${tumor}__${normal} in ${total_time_in_days} days" | tee -a main.log
    # check if all samples finished
    finished=$( cat finished.csv | wc -l )
    started=$( cat tumors_and_normals.csv | grep -v "^#" | wc -l )
    if [[ "$finished" -eq "$started" ]]; then
        # add tmb_and_coverage to archive
        zip -ru analyses.zip analyses/coverage_and_tmb.csv
        # tidyup and clean working dir
        if [[ ! -e bam ]]; then
            mv BQSR bam
        else
            mv BQSR/* bam
            rm -rf BQSR
        fi
        # move stats to all_logfiles
        mv mutect2/*.filteringStats.tsv all_logfiles
        #mv mutect2/*.mutect2.unfiltered.wes.merged.vcf.stats all_logfiles
        # delete all other vcf files except filtered and unfiltered VCFs
        rm $(ls mutect2/*vcf* | grep -v "filtered")
        # bgzip and tabix all vcf files
        ls mutect2/*.vcf | parallel --tmpdir ./tmp "index-vcf {}"
        # delete directories with bam data
        rm -rf preprocessed_bam aligned_bam tmp
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
            echo -e "10: \npipeline took ${total_time_in_days} days for batch to complete" | tee -a main.log
            # change/adjust permisions
            # this configuration allows the main user and the users in the tabori group to
            # read/write/excecute
            # dirs first
            echo "10: changing permissions" | tee -a main.log
            chmod 774 all_logfiles \
                      analyses \
                      bam \
                      contamination \
                      mutect2/f1r2 \
                      vcf/snpEff
            # files second
            chmod 664 all_logfiles/* \
                      analyses/* \
                      bam/* \
                      contamination/* \
                      mutect2/*vcf* \
                      mutect2/f1r2/* \
                      vcf/*vcf* \
                      vcf/*maf* \
                      vcf/snpEff/* \
                      varscan/*vcf* \
                      ${tumor}__${normal}.analyses.log
            #
        fi
    fi

    # # run mutational signature analysis
    # Rscript ${pipeline_dir}/cosmic_signature_analysis.R ${mode} ${tumor}__${normal}.

    # last move logfile to dir
    if [[ -e ${tumor}__${normal}.analyses.log ]]; then
        mv ${tumor}__${normal}.analyses.log all_logfiles
    fi
fi
