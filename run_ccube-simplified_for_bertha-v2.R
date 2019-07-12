#!/usr/bin/env Rscript

library(ccube)
library(dplyr) # required for Ccube to work
library(foreach) # required for Ccube to work
library(GenomicRanges)
library(stringr)
library(VariantAnnotation)

# CCUBE WRAPPER SCRIPT FOR RUNNING ON GENOMICS ENGLAND LSF
# THIS IS AN ADAPTED VERSION OF CODE BY ALEX CORNISH AT THE ICR
# WHICH WAS WRITTEN FOR THE GENOMICS ENGLAND RESEARCH ENVIRONMENT
# FOR RUNNING USING DIFFERENT INPUTS
# THIS CODE HAS BEEN TESTED FOR RENAL CANCER, COLORECTAL CANCER COHORTS
# AND TO A LESSER EXTENT TRACERx LUNG CANCER SAMPLES.

# SEE https://cnfl.extge.co.uk/pages/viewpage.action?spaceKey=BTS&title=Ccube+Tumour+Purity+Estimate
# FOR DOCUMENTATION OF SCIENTIFIC TESTING

# CCUBE USES AS INPUTS SNV AND CNV CALLS, GENERATED FOR EXAMPLE BY STRELKA AND CANVAS SOFTWARE
# FROM ILLUMINA. CCUBE USES THIS INFORMATION TO IDENTIFY AREAS IN THE GENOME WHERE IT WILL PERFORM
# CALCULATIONS I.E. WHERE THESE SNVS AND CNVS OVERLAP. IT ALSO TAKES FROM CANVAS A TUMOUR PURITY
# ESTIMATE AND PLOIDY VALUE. THIS IS STILL USEFUL, AS ALTHOUGH CCUBE TAKES IN A TUMOUR PURITY
# VALUE TO PRODUCE A TUMOUR PURITY ESTIMATE, WE ALREADY HAVE A CANVAS TUMOUR PURITY ESTIMATE,
# SO THERE IS NO COST IN OBTAINING THIS, BUT WE ARE NOT CONFIDENT IN THIS VALUE, IN PARTICULAR,
# VALUES OF <= 0.2 ARE INDICATED TO BE POOR PREDICTIONS. THIS IS WHY WE ARE USING CCUBE AS
# AN ALTERNATIVE. SEE THE CONFLUENCE WEB PAGE GIVEN ABOVE FOR FURTHER INFORMATION.

# WE FOUND THAT IT IS IMPORTANT TO USE FILTERED STRELKA INPUT FILES. RESULTS ARE NOT RELIABLE
# USING UNFILTERED STRELKA INPUTS. IN PARTICULAR WE USED THIS SCRIPT & INPUTS FOR FILTERING
# STRELKA SNV INPUT FILES PRIOR TO USE:

#///home/jambrose/jca_scripts/filter_variants_plus_gnomad.py
#GERMLINE_FREQ=0.01
#GNOMAD_FREQ=0.01
#GNOMAD_STUDY_TYPE=GNOMAD_GENOMES
#The above were inputs in the wrapper script which feeds them into the python script:
# /home/jambrose/jca_scripts/filtering_somatic_variants_gnomad_bsub.sh

# THE ORIGINAL CODE BY ALEX CORNISH IS HERE:
# /home/pcarter/bin/ccube/Alex_code/

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 8){
    SNV_file_name <- args[1]
    CNV_file_name <- args[2]
    tumour_sample_name <- args[3]
    normal_sample_name <- args[4]
    output_results_file_name <- args[5]
    output_PDF_file_name <- args[6]
    output_RData_file_name <- args[7]
    SV_vcf_bed_output_file <- args[8]
} else {
    cat(paste0("Please provide 10 arguments in this order: SNV_file CNV_file tumour_sample_name normal_sample_name output_file output_PDF output_RData SV_vcf_bed_output_file\n"))
    quit(status=1)
}

cat(paste0("SNV file (FILTERED Strelka) = ", SNV_file_name, "\n"))
cat(paste0("CNV file (Canvas)           = ", CNV_file_name, "\n"))
cat(paste0("Tumour sample name          = ", tumour_sample_name, "\n"))
cat(paste0("Normal sample name          = ", normal_sample_name, "\n"))
cat(paste0("Output values file name     = ", output_results_file_name, "\n"))
cat(paste0("Output PDF file name        = ", output_PDF_file_name, "\n"))
cat(paste0("Output RData file name      = ", output_RData_file_name, "\n"))
cat(paste0("SV vcf bed output file name = ", SV_vcf_bed_output_file, "\n"))

if(!file.exists(SNV_file_name)){
    cat(paste0("SNV (Strelka) file not found.\n"))
    quit(status=1)
}

if(!file.exists(CNV_file_name)){
    cat(paste0("CNV (Canvas) file not found.\n"))
    quit(status=1)
}

# THIS FUNCTION IS USED TO GET CLUSTER VARIANCE, WHICH THEORETCIALLY
# CAN BE USED AS A CONFIDENCE SCORE FOR THE PURITY ESTIMATE VALUE.
# HOWEVER, IN PRACTICE THE RESULTS DO NOT LOOK USEFUL, AND SO MAY
# NEED FURTHER CHECKING WITH THE ORIGINAL CCUBE AUTHOR (KE YUAN AT GLASGOW UNIVERSITY)
# FOR A BETTER IMPLEMENTATION.
# THIS FUNCTION IS MOSTLY KE YUAN'S CODE
GetPurityModified <- function (mydata, wgd = F, K = 6, th = 0.015){

    vtox <- function(v, nA, nB, tA, tB) {
        (nB - nB * v - nA * v)/(tA * v + tB * v + nB - tB - nA * v - nB * v)
    }
    tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn != 0)
    if (!wgd) {
        tmpdata <- dplyr::filter(mydata, major_cn == minor_cn & major_cn == 1)
    }
    if (nrow(tmpdata) == 0) {
        return(NA)
    }
    tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), ploidy = major_cn + minor_cn)
    tmpdata <- dplyr::mutate(tmpdata, vaf = var_counts/(var_counts + ref_counts))
    tmpdata <- dplyr::mutate(dplyr::rowwise(tmpdata), cp = vtox(vaf, 2, 0, ploidy/2, ploidy/2))
    tmpdata <- dplyr::filter(tmpdata, !is.infinite(cp) & !is.na(cp))
    if (K >= nrow(tmpdata)) {
        K = nrow(tmpdata) - 1
        if (K == 0) {
            K = 1
        }
    }
    if (nrow(tmpdata) == 1) {
        purity <- tmpdata$cp
	if(purity > 1){
            purity <- 1
        }
    }
    else {
        res <- vbsmm(tmpdata$cp, init = K, tol = 1e-05, verbose = F)

        uniqLabels <- sort(unique(res$label))
        ww <- uniqLabels[which((table(res$label)/length(res$label)) > th)]
        pool <- res$mu[ww]
        maxCp <- max(pool[pool <= 1])
        if(maxCp > 1){
            purity <- 1
        } else {
	    purity <- maxCp
        }
    }

    return(purity)
}

pastedot <- function(...) paste(..., sep=".")

# CHANGE TO 0 IF USING NEWER VERSION OF CANVAS:
#use_NSV4_Canvas <- 0
use_NSV4_Canvas <- 1

plateid.normal <- normal_sample_name
plateid.tumour <- tumour_sample_name
combined_IDs <- paste0(normal_sample_name, "_", tumour_sample_name)

strelka_snv_file <- SNV_file_name
canvas_cnv_file <- CNV_file_name

# READ IN THE CALLS
vcf_snvs <- readVcf(strelka_snv_file, genome="hg38")
vcf_cnvs <- readVcf(canvas_cnv_file, genome="hg38")
vcf_cnv_hdr <- header(vcf_cnvs)

purity_column_number = 0
ploidy_column_number = 0
purity_column_number <- which(rownames(vcf_cnv_hdr@header[[1]]) == "EstimatedTumorPurity")
ploidy_column_number <- which(rownames(vcf_cnv_hdr@header[[1]]) == "OverallPloidy")
if((purity_column_number == 0) || (ploidy_column_number == 0)){
    cat("Couldn't find a VCF column number.\n")
    quit(status=1)
}
# NEEDED AS INPUTS TO CCUBE
vcf_purity   <- as.numeric(vcf_cnv_hdr@header[[1]][purity_column_number,])
vcf_ploidy   <- as.numeric(vcf_cnv_hdr@header[[1]][ploidy_column_number,])

cnv_bed_info <- read.table(SV_vcf_bed_output_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
cnv_bed_info_with_MCC <- cnv_bed_info[which(!is.na(cnv_bed_info[,11])),]
cnv_bed_info_with_MCC2 <- cnv_bed_info_with_MCC[which(cnv_bed_info_with_MCC[,11] != "None"),]
cnv_bed_info <-  cnv_bed_info_with_MCC2

colnames(cnv_bed_info) <- c("Chr", "Start_pos", "Stop_pos", "Call_type", "Quality", "Sample1_genos", "Filter", "Alt", "Geno_RC", "Geno_CN", "Geno_MCC", "ID")
CNV_gain_or_loss_positions <- which(cnv_bed_info[,8] == "CNV")
copy_number_chromosomes     <- cnv_bed_info[,1]
copy_number_start_positions <- cnv_bed_info[,2]
copy_number_end_positions   <- cnv_bed_info[,3]
total_copy_numbers <- cnv_bed_info[,10]
major_copy_numbers <- cnv_bed_info[,11]

minor_copy_numbers <- total_copy_numbers - major_copy_numbers

canvas_copy_numbers_df <- data.frame(cnv_bed_info[,1], cnv_bed_info[,2], cnv_bed_info[,3], total_copy_numbers, major_copy_numbers, minor_copy_numbers)
colnames(canvas_copy_numbers_df) <- c("Chromosome", "Start_pos", "End_pos", "Total_copy_number", "Major_copy_number", "Minor_copy_number")

# SOME CCUBE PARAMETERS & SETTING UP
num.of.cluster.pool <- 1:6
num.of.repeat <- 1
maxiter <- 10000
set.seed(1234) # for reproducability
chrs <- paste0("chr", 1:22) # chromosomes to consider
nucleotides <- c("A", "C", "G", "T")
plateids.comb <- paste0("tumo", plateid.tumour, "_norm", plateid.normal)

filename.purityest     <- output_results_file_name
filename.ccube.results <- output_RData_file_name
filename.plot          <- output_PDF_file_name

wgd <- round(vcf_ploidy) >= 4 # ALEX TO DO: should probably match with a high-confidence list of WGD samples, but not currently available

# FORMAT SNV DATA
vcf_pass_pos <- (rowRanges(vcf_snvs)$FILTER == "PASS")
vcf_snvs_pass <- vcf_snvs[vcf_pass_pos,] # remove non-PASS variants
vcf_ref_no_indels_pos <- which(str_length(as.character(ref(vcf_snvs_pass))) == 1)

# FIXING FOR LIST OF LISTS INDEX > 1
vcf_alt_multi_lists_pos <- which(lengths(alt(vcf_snvs_pass)) > 1)

snvs_pass_limit <- length(alt(vcf_snvs_pass))
vcf_alt_indels_pos <- vector()

for(i in 1:snvs_pass_limit){
    this_str_length <- nchar(alt(vcf_snvs_pass)[[i]])
    if(this_str_length > 1){
        vcf_alt_indels_pos <- c(vcf_alt_indels_pos, i)
    }
}

vcf_all_indels_pos    <- unique(sort(c(vcf_alt_multi_lists_pos, vcf_alt_indels_pos)))
vcf_alt_no_indels_pos <- which(!1:max(vcf_all_indels_pos) %in% vcf_all_indels_pos)
vcf_ref_and_alt_no_indels_pos <- intersect(vcf_ref_no_indels_pos, vcf_alt_no_indels_pos)
vcf_snvs_no_indels <- vcf_snvs_pass[vcf_ref_and_alt_no_indels_pos,] # remove indels
canvas_copy_numbers_df_chr_filtered <- canvas_copy_numbers_df[canvas_copy_numbers_df$Chromosome %in% chrs,]

canvas_copy_numbers_df_GRanges <- GRanges(seqnames = canvas_copy_numbers_df$Chromosome,
                                          ranges = IRanges(start = canvas_copy_numbers_df$Start_pos, end = canvas_copy_numbers_df$End_pos),
                                          minor_cn = canvas_copy_numbers_df$Minor_copy_number, major_cn = canvas_copy_numbers_df$Major_copy_number)

vcf_cnv_and_snv_overlap <- findOverlaps(canvas_copy_numbers_df_GRanges, vcf_snvs_no_indels)
vcf_cnv_and_snv_overlap_no_dups <-
    vcf_cnv_and_snv_overlap[!duplicated(subjectHits(vcf_cnv_and_snv_overlap)), ] # remove SNVs overlapping more than one CNV segment
                                                                                 # (there shouldn't be any, but just to sure)
vcf_snvs_no_indels.overlap <- vcf_snvs_no_indels[subjectHits(vcf_cnv_and_snv_overlap_no_dups),]
canvas_copy_numbers_df_GRanges.overlap <- canvas_copy_numbers_df_GRanges[queryHits(vcf_cnv_and_snv_overlap_no_dups),]
snvs.overlap <- vcf_snvs_no_indels.overlap
segs.overlap <- canvas_copy_numbers_df_GRanges.overlap

ccube.input <- data.frame(
    mutation_id=names(snvs.overlap),
    minor_cn=segs.overlap$minor_cn,
    major_cn=segs.overlap$major_cn,
    total_cn=segs.overlap$minor_cn + segs.overlap$major_cn,
    stringsAsFactors=F
)

# ADD READ COUNTS
read.depths <- do.call("cbind", sapply(paste0(nucleotides, "U"), function(field) geno(snvs.overlap)[[field]][,, 1], simplify=F))
#read.depths.tumour <- read.depths[, colnames(read.depths) == plateid.tumour]
read.depths.tumour <- read.depths
colnames(read.depths.tumour) <- nucleotides
ccube.input$var_counts <- read.depths.tumour[cbind(rownames(read.depths.tumour), as.character(CharacterList(alt(snvs.overlap))))]
ccube.input$ref_counts <- read.depths.tumour[cbind(rownames(read.depths.tumour), as.character(ref(snvs.overlap)))]
ccube.input$total_counts <- ccube.input$var_counts + ccube.input$ref_counts

ccube.input$purity <- vcf_purity
ccube.input$normal_cn <- 2

# REMOVE ANY DUPLICATE MUTATIONS
duplicated_ID_rows <- which(duplicated(ccube.input[,1]))
number_of_duplicated_ID_rows <- length(duplicated_ID_rows)
if(number_of_duplicated_ID_rows > 0){
    cat(paste0(number_of_duplicated_ID_rows, " duplicated mutations were removed\n"))
    ccube.input <- ccube.input[-c(duplicated_ID_rows),]
}

cat(paste0("Preprocessing completed.\n"))

# RUN CCUBE
cat(paste0("Running Ccube.\n"))
# ccube input should be a data frame containing columns for: mutation_id, minor_cn, major_cn, total_cn, purity, normal_cn, total_counts, var_counts, ref_counts
ccube.results <- NULL
ccube.results <- RunCcubePipeline(
    sampleName = plateid.tumour,
    runAnalysis = T,
    runQC = T,
    purity = purity,
    numOfClusterPool = num.of.cluster.pool,
    numOfRepeat = num.of.repeat,
    maxiter = maxiter,
    ssm = ccube.input#,
    #returnAll = T # CHANGE FOR MORE RECENT VERSION OF CCUBE ON BERTHA
)
if(is.null(ccube.results)){
    cat("Ccube failed.\n")
    purityest_val <- 'FAIL'
    write.table(data.frame(purity=purityest_val), file=filename.purityest, quote=F, sep="\t", row.names=F, col.names=T)
    quit(status=1)
}
		
num.of.cluster.opt <- num.of.cluster.pool[order(ccube.results$lb, decreasing=T)][1] # number of clusters in optimal model solution
cat(paste0("Ccube completed.\n"))
save(ccube.results, file = filename.ccube.results)
cat(paste0("Ccube results saved.\n"))

# GET PURITY ESTIMATE USING OPTIMAL NUMBER OF CLUSTERS OBTAINED FROM INITIAL CCUBE RUN
cat(paste0("Getting purity estimate.\n"))
purityest_val <- GetPurityModified(
    mydata = ccube.input,
    wgd = wgd,
    K = num.of.cluster.opt
)
cat(paste0("Purity estimate complete.\n"))
options(scipen=999) # PREVENT SCIENTIFIC NOTATION FOR RESULT
show(as.numeric(purityest_val))
#write.table(data.frame(purity=purityest_vals[1], cluster_variance=purityest_vals[2]), file=filename.purityest, quote=F, sep="\t", row.names=F, col.names=T)
write.table(data.frame(purity=purityest_val), file=filename.purityest, quote=F, sep="\t", row.names=F, col.names=T)
cat(paste0("Purity estimate result saved.\n"))

# PLOT CCUBE RESULTS
cat(paste0("Plotting results.\n"))
MakeCcubeStdPlot(
    ssm=ccube.results$ssm,
    res=ccube.results$res,
    printPlot=T,
    fn=filename.plot
)
cat(paste0("Plotting completed.\n"))

cat(paste0("All tests completed.\n"))

