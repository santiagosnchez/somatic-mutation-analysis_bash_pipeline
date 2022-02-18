# from sigminer
library(stringr)

scan_vcf_to_maf <- function(vcf_file){
  # get file
  vcf = scan(gzfile(vcf_file), what=character(), sep="\n", quiet=T)
  # get funcotator header
  func_head = vcf[ grep("^##INFO=<ID=FUNCOTATION", vcf) ]
  func_head = sub(".*Funcotation fields are: ","",func_head)
  func_head = sub("\\\">","",func_head)
  func_head = strsplit(func_head, "\\|")[[1]]
  main_maf_heads = c(Gencode_34_hugoSymbol="Hugo_Symbol",
                     Gencode_34_ncbiBuild="NCBI_Build",
                     Gencode_34_chromosome="Chromosome",
                     Gencode_34_start="Start_Position",
                     Gencode_34_end="End_Position",
                     Gencode_34_variantClassification="Variant_Classification",
                     Gencode_34_secondaryVariantClassification="Secondary_Variant_Classification",
                     Gencode_34_variantType="Variant_Type",
                     Gencode_34_refAllele="Reference_Allele",
                     Gencode_34_tumorSeqAllele1="Tumor_Seq_Allele1",
                     Gencode_34_tumorSeqAllele2="Tumor_Seq_Allele2",
                     Gencode_34_genomeChange="Genome_Change",
                     Gencode_34_annotationTranscript="Annotation_Transcript",
                     Gencode_34_transcriptStrand="Strand",
                     Gencode_34_transcriptExon="Transcript_Exon",
                     Gencode_34_transcriptPos="Transcript_position",
                     Gencode_34_cDnaChange="cDNA_change",
                     Gencode_34_codonChange="Codon_change",
                     Gencode_34_proteinChange="Protein_change",
                     Gencode_34_gcContent="GC_content",
                     Gencode_34_referenceContext="Reference_Context",
                     Gencode_34_otherTranscripts="Other_Transcripts")
  new_maf_heads=c(main_maf_heads[func_head[1:22]], func_head[22:length(func_head)])
  # process only variant calls in vcf
  vcf = vcf[ grep("^#", vcf, invert=T) ]
  # extract INFO l
}


read_vcf <- function (vcfs, samples = NULL, genome_build = c("hg19", "hg38", "mm10", "mm9"), keep_only_pass = FALSE, verbose = TRUE){
  genome_build <- match.arg(genome_build)
  # maf header
  vcf = scan(gzfile(vcf_with_filter[1]), what=character(), sep="\n", quiet=T)
  func_head = func_head[ grep("^##INFO=<ID=FUNCOTATION", func_head) ]
  func_head = sub(".*Funcotation fields are: ","",func_head)
  func_head = sub("\\\">","",func_head)
  func_head = strsplit(func_head, "\\|")[[1]]
  main_maf_heads = c(Gencode_34_hugoSymbol="Hugo_Symbol",
                     Gencode_34_ncbiBuild="NCBI_Build",
                     Gencode_34_chromosome="Chromosome",
                     Gencode_34_start="Start_Position",
                     Gencode_34_end="End_Position",
                     Gencode_34_variantClassification="Variant_Classification",
                     Gencode_34_secondaryVariantClassification="Secondary_Variant_Classification",
                     Gencode_34_variantType="Variant_Type",
                     Gencode_34_refAllele="Reference_Allele",
                     Gencode_34_tumorSeqAllele1="Tumor_Seq_Allele1",
                     Gencode_34_tumorSeqAllele2="Tumor_Seq_Allele2",
                     Gencode_34_genomeChange="Genome_Change",
                     Gencode_34_annotationTranscript="Annotation_Transcript",
                     Gencode_34_transcriptStrand="Strand",
                     Gencode_34_transcriptExon="Transcript_Exon",
                     Gencode_34_transcriptPos="Transcript_position",
                     Gencode_34_cDnaChange="cDNA_change",
                     Gencode_34_codonChange="Codon_change",
                     Gencode_34_proteinChange="Protein_change",
                     Gencode_34_gcContent="GC_content",
                     Gencode_34_referenceContext="Reference_Context",
                     Gencode_34_otherTranscripts="Other_Transcripts")
  new_maf_heads=c(main_maf_heads[func_head[1:22]], func_head[22:length(func_head)])

  # read data from vcfs
  vcfs_name <- vcfs
  if (verbose)
    message("Reading file(s): ", paste(vcfs, collapse = ", "))
  vcfs <- purrr::map(vcfs, ~tryCatch(data.table::fread(., select = 8, skip = "#CHROM"),
          error = function(e) {
          message("It seems ", ., " has no normal VCF header, try parsing without header.")
          data.table::fread(., select = 8)
          }))
