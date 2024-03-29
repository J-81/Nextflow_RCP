#! /usr/bin/env Rscript

# TODO: DISCUSS: new section
## Set options

# allow more flexibility in download time
# useful for slower connections where the default of 60 seconds might be exceeded
options(timeout=600)

## Import libraries (tximport, DESeq2, tidyverse, Risa)
library(tximport)
library(DESeq2)
library(tidyverse)

################################################################################
# Retrieve args from commandline
################################################################################
args = commandArgs(trailingOnly=TRUE)
print("ARGS")
print(args)

## Define which organism is used in the study - this should be consistent with the name in the organisms.csv file, which matches the abbreviations used in the Panther database for each organism
organism <- args[1]


## Create a directory for the metadata and in that directory, download the *ISA.zip file for the study being analyzed, which is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files'

#MOD#metadata_dir=arg[2]
runsheet_path=args[2]
#END_MOD#
work_dir=getwd()
#REMOVE_MOD#counts_dir=args[4]
norm_output=args[3]
DGE_output=args[4]
DEBUG_RUN=args[5]
DGE_output_ERCC=args[6]

if (DEBUG_RUN) {
	print('DEBUG MODE ENABLED')
} else {
	print('DEBUG MODE DISABLED')
}

if (!is.na(DGE_output_ERCC)) {
	print('ERCC MODE ENABLED')
	ERCC_MODE=TRUE
} else {
	print('ERCC MODE DISABLED')
	ERCC_MODE=FALSE
}

####### Pull all factors for each sample in the study from the metadata in the runsheet 
compare_csv_from_runsheet <- function(runsheet_path) {
	df = read.csv(runsheet_path)
	# get only Factor Value columns
	print(colnames(df))
	factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
	colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")

	result = data.frame(sample_id = df[,c("Sample.Name")], factors)	
	return(result)
}

compare_csv <- compare_csv_from_runsheet(runsheet_path)

#### Create data frame containing all samples and respective factors
study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]

### If pulling from a metadata file rather than ISA.zip use the command below instead
#study <- read.csv(Sys.glob(file.path(work_dir,"*metadata.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)


# ## Set your working directory to the directory containing the organisms.csv file
#
# setwd(file.path(work_dir))


##### Format groups and indicate the group that each sample belongs to #####
if (dim(study) >= 2){
	group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
	group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") # human readable group names
group <- make.names(group) # group naming compatible with R models
names(group) <- group_names
rm(group_names)


##### Format contrasts table, defining pairwise comparisons for all groups #####
contrasts <- combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))),2)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names)

##### Import RSEM raw (gene) count data #####
files <- list.files(path = "Rsem_gene_counts", pattern = ".genes.results", full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples

# Replace spaces in sample names from ISA with "_", consistent with runsheet generation
library(stringr)
samples = str_replace_all(rownames(study), " ", "_")
reordering <- sapply(samples, function(x)grep(paste0(x,".genes.results$"), files, value=FALSE))
files <- files[reordering]
names(files) <- samples
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)


## Add 1 to genes with lengths of zero - needed to make DESeqDataSet object
txi.rsem$length[txi.rsem$length == 0] <- 1

## DEBUG MODE: ensures truncated runs process by supplying random gene counts for each sample
if (DEBUG_RUN) {
	print('DEBUG MODE: Replacing all gene counts with random counts from 0-5000!!')
	set.seed(77077)
	txi.rsem$counts <- txi.rsem$counts + matrix(sample( 0:5000, NROW(txi.rsem$counts)*NCOL(txi.rsem$counts), replace=TRUE),nrow=NROW(txi.rsem$counts))
}
##### Make DESeqDataSet object #####

## Create data frame defining which group each sample belongs to
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
summary(dds)


##### Create ERCC unfiltered raw counts table #####

## Note: These data are used internally at GeneLab for QC

if (ERCC_MODE) {
	## Make a DESeqDataSet object using only unfiltered ERCC genes
	ercc_rows_all <- grep("ERCC-",rownames(dds))
	ercc_dds_all <- dds[ercc_rows_all,]
	
	## Print ERCC unfiltered raw counts table
	ERCC_rawCounts_all = as.data.frame(counts(ercc_dds_all))
	write.csv(ERCC_rawCounts_all,file='ERCC_rawCounts_unfiltered.csv')
}

##### Filter out genes with counts of less than 10 in all samples #####
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
summary(dds)
dim(dds)

##### Prepare filtered data to be normalized with and without considering ERCC genes #####

if (ERCC_MODE) {
	## Make a DESeqDataSet object using only filtered ERCC genes, which will be used to generate ERCC counts table
	ercc_rows <- grep("ERCC-",rownames(dds))
	ercc_dds <- dds[ercc_rows,]

	## Create list of rows containing ERCC group B genes to use for ERCC-normalization
	## Note: ERCC group B genes should be the same concentration in all samples
	ercc_rows_gpB <- grep("ERCC-00096|ERCC-00171|ERCC-00009|ERCC-00042|ERCC-00060|ERCC-00035|ERCC-00025|ERCC-00051|ERCC-00053|ERCC-00148|ERCC-00126|ERCC-00034|ERCC-00150|ERCC-00067|ERCC-00031|ERCC-00109|ERCC-00073|ERCC-00158|ERCC-00104|ERCC-00142|ERCC-00138|ERCC-00117|ERCC-00075",rownames(dds))
	ercc_dds_gpB <- dds[ercc_rows_gpB,]
	summary(ercc_dds_gpB)
	dim(ercc_dds_gpB)

	## Identify and list samples that do not contain counts for ERCC genes
	## This is to make sure all samples indeed contain ERCC spike-in
	cat("Samples that do not have detectable ERCC spike-ins: ", colnames(ercc_dds[,colSums(counts(ercc_dds))==0]), sep="\n")
	cat("Samples that do not have detectable ERCC group B spike-ins: ", colnames(ercc_dds[,colSums(counts(ercc_dds_gpB))==0]), sep="\n")


	# Create a new group object WITHOUT the samples that don't have detectable ERCC group B spike-ins
	remove <- colnames(dds[,colSums(counts(ercc_dds_gpB))==0])

	study_sub <- subset(study,!rownames(study) %in% remove) # new study object with non-ERCC-gpB samples removed

	## New group object with with non-ERCC-gpB samples removed
	if (dim(study_sub) >= 2){
		group_sub<-apply(study_sub,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
	} else{
		group_sub<-study_sub[,1]
	}
	group_names <- paste0("(",group_sub,")",sep = "") # human readable group names
	group_sub <- make.names(group_sub) # group naming compatible with R models
	names(group_sub) <- group_names
	rm(group_names)

	## New contrasts object that only reflects the groups in the subset group object
	contrasts_sub <- combn(levels(factor(group_sub)),2) # generate matrix of pairwise group combinations for comparison
	contrasts_sub.names <- combn(levels(factor(names(group_sub))),2)
	contrasts_sub.names <- c(paste(contrasts_sub.names[1,],contrasts_sub.names[2,],sep = "v"),paste(contrasts_sub.names[2,],contrasts_sub.names[1,],sep = "v")) # format combinations for output table files names
	contrasts_sub <- cbind(contrasts_sub,contrasts_sub[c(2,1),])
	colnames(contrasts_sub) <- contrasts_sub.names
	rm(contrasts_sub.names)

	## If no samples to remove, reassign group_sub, study_sub and contrasts_sub back to the original variable contents
	if (length(remove) == 0) {
		group_sub <- group
		study_sub <- study
		contrasts_sub <- contrasts
	}

	## Print ERCC filtered raw counts table
	## Note: These data are used internally at GeneLab for QC
	ERCC_rawCounts = as.data.frame(counts(ercc_dds))
	write.csv(ERCC_rawCounts,file='ERCC_rawCounts_filtered.csv')

	## Remove samples that do not contain ERCC counts
	## Note: All samples should contain ERCC spike-in and thus ERCC counts, if some samples do not contain ERCC counts, those samples should be removed and not used for downstream analysis
	ercc_dds <- ercc_dds[,colSums(counts(ercc_dds)) > 0]
}
##### Generate a DESeqDataSet object using only non-ERCC genes #####

## dds_1 will be used to generate data without considering ERCC genes
if (ERCC_MODE) {
	dds_1 <- dds[-c(ercc_rows),] # remove ERCCs from full counts table
	## dds_2 will be used to generate data with considering ERCC genes
	## First remove samples that do not contain ERCC counts
	## Note: All samples should contain ERCC spike-in and thus ERCC counts, if some samples do not contain ERCC counts, those samples should be removed and not used for downstream analysis
	dds_2 <- dds[,colSums(counts(ercc_dds_gpB)) > 0] 
	sampleTable_sub <- data.frame(condition=factor(group_sub)) # create a new sampleTable with only the kept samples
	dds_2$condition <- sampleTable_sub$condition # reassign the condition to the subset condition to keep the groups and levels consistent with the "keep" samples
	summary(dds_2)
	dim(dds_2)
} else {
	ercc_rows <- grep("ERCC-",rownames(dds))
	if ( length(ercc_rows) != 0 )
		{ stop("Found ERCC rows despite user indicating NO ERCC") }
	dds_1 <- dds
}


##### Perform DESeq analysis with and without considering ERCC genes #####
if (ERCC_MODE) {
	## Run DESeq analysis with ERCC-normalization by replacing size factor object with ERCC size factors for rescaling
	## Try first to use the default type="median", but if there is an error (usually due to zeros in genes), use type="poscounts"
	## From DESeq2 manual: "The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts."

	dds_2 <- tryCatch(
			  expr = { estimateSizeFactors(dds_2, controlGenes=ercc_rows_gpB) },
			  error = function(e) { estimateSizeFactors(dds_2, type="poscounts", controlGenes=ercc_rows_gpB)}
	)

	dds_2 <- dds_2[-c(ercc_rows),] # remove ERCCs from counts table after normalization
	dds_2 <- estimateDispersions(dds_2)
	dds_2 <- nbinomWaldTest(dds_2)
}
## Run DESeq analysis without considering ERCC genes
dds_1 <- DESeq(dds_1)

##### Export unnormalized, normalized, and ERCC-normalized counts as well as the sample table #####

normCounts = as.data.frame(counts(dds_1, normalized=TRUE))
setwd(file.path(norm_output))
write.csv(txi.rsem$counts,file='RSEM_Unnormalized_Counts.csv')
write.csv(normCounts,file='Normalized_Counts.csv')
write.csv(sampleTable,file='SampleTable.csv')

if (ERCC_MODE) {
	ERCCnormCounts = as.data.frame(counts(dds_2, normalized = TRUE))
	write.csv(ERCCnormCounts, file ='ERCC_Normalized_Counts.csv')
	write.csv(sampleTable_sub, file = 'ERCC_SampleTable.csv')
}
setwd(file.path(work_dir))


##### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design #####

## Add 1 to all counts to avoid issues with log transformation
normCounts <- normCounts +1
if (ERCC_MODE) {
	ERCCnormCounts <- ERCCnormCounts +1
}
## For non-ERCC normalized data
dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

if (ERCC_MODE) {
	## For ERCC-normalized data
	dds_2_lrt <- DESeq(dds_2, test = "LRT", reduced = ~ 1)
	res_2_lrt <- results(dds_2_lrt)
}

##### Generate annotated DGE tables #####

## Import table with organism db objects for annotation
organism_table <- read.csv(file.path(work_dir,"organisms.csv"))

## Load annotation libraries
## If working through VPN, internet via proxy may need to be defined
assign("has_internet_via_proxy", TRUE, environment(curl::has_internet))
library(STRINGdb) # for String database annotations
library(PANTHER.db) # for GOSLIM annotations

## Begin building anotation database
ann.dbi <- organism_table$annotations[organism_table$name == organism] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
  BiocManager::install(ann.dbi, ask = FALSE)
  library(ann.dbi, character.only=TRUE)
}


#### Generate annotated DGE table for non-ERCC normalized counts ####

## Start by creating output tables with (non-ERCC) normalized sample expression values

## reduced output table 1 will be used to generate human-readable DGE table
reduced_output_table_1 <- normCounts

## output table 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
output_table_1 <- normCounts
print(colnames(output_table_1))

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in 1:dim(contrasts)[2]){
  res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
  res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
  colnames(res_1)<-c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("Stat_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
  output_table_1<-cbind(output_table_1,res_1)
  reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
  rm(res_1)
}

## Determine the keytype to use for annotation database
## Order (proceed to next if not found)
## 'ENSEMBL' -> 'TAIR' -> Raise Error
if ( "ENSEMBL" %in% keytypes( eval( parse(text = ann.dbi) ) ) ) {
  keytype = "ENSEMBL"
} else if ( "TAIR" %in% keytypes( eval( parse(text = ann.dbi) ) ) ) {
  keytype = "TAIR"
} else {
  stop(
    sprintf("Neither 'ENSEMBL' nor 'TAIR' keytypes found in %s. The following keytypes were found '%s'",
            ann.dbi,
            paste( keytypes( eval( parse(text = ann.dbi) ) ), collapse="' '" )
          )
      )
}
sprintf("Using '%s' keytype for %s", keytype, ann.dbi )

## Create annotation table and add gene annotation columns
#keytype = "ENSEMBL" # will be different if primary annotations are not ENSEMBL # replaced by section above
annot <- data.frame(rownames(output_table_1), stringsAsFactors = FALSE)
colnames(annot)[1]<-keytype
if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
  annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "SYMBOL", multiVals = "first")
}
if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
  annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "GENENAME", multiVals = "first")
}
if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
  annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENSEMBL", multiVals = "first")
}
if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
  annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "REFSEQ", multiVals = "first")
}
if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
  annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENTREZID", multiVals = "first")
}

## Create and add string annotation columns to the annotation table
string_db <- STRINGdb$new( version="11", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
string_map <- string_map[!duplicated(string_map$SYMBOL),]
annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")

## Create and add columns containing GOSLIM ids using the panther annotation database to the annotation table
pthOrganisms(PANTHER.db) <- organism
panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
panther <- na.omit(panther)
annot$GOSLIM_IDS <- panther
rm(string_db,string_map,panther)

## Generate and add all sample mean column to the (non-ERCC) normalized counts table
output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the (non-ERCC) normalized counts table
output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the (non-ERCC) normalized counts table
output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

## Generate and add group mean and stdev columns to the (non-ERCC) normalized counts table
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- group
group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
group_means <- group_means[-c(1),]
colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))
group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
group_stdev <- group_stdev[-c(1),]
colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))
output_table_1 <- cbind(output_table_1,group_means)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means)

output_table_1 <- cbind(output_table_1,group_stdev)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_stdev)

rm(group_stdev,group_means,tcounts)


### Add columns needed to generate GeneLab visulaization plots to the (non-ERCC) normalized counts table
## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)
## Add column to indicate contrast significance with p <= 0.1
sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05
sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)
print(colnames(output_table_1))
## Add columns for the volcano plot with p-value and adjusted p-value
log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)
### Combine annotations table and the (non-ERCC) normalized counts table

output_table_1 <- cbind(annot,output_table_1)
reduced_output_table_1 <- cbind(annot,reduced_output_table_1)
rownames(output_table_1) <- NULL
rownames(reduced_output_table_1) <- NULL
output_table_1$GOSLIM_IDS <- vapply(output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))
reduced_output_table_1$GOSLIM_IDS <- vapply(reduced_output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))

### Export human- and computer/visualization- readable DGE tables

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))
write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)

### Generate and export PCA table for GeneLab visualization plots

exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)
rm(exp_raw,PCA_raw)

#### Generate annotated DGE table for ERCC-normalized counts ####
if (ERCC_MODE) {
	## Start by creating output tables with ERCC-normalized sample expression values

	## reduced output table 2 will be used to generate human-readable DGE table
	reduced_output_table_2 <- ERCCnormCounts

	## output table 2 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
	output_table_2 <- ERCCnormCounts

	## Iterate through Wald Tests to generate pairwise comparisons of all groups
	for (i in 1:dim(contrasts_sub)[2]){
		res_2 <- results(dds_2, contrast=c("condition",contrasts_sub[1,i],contrasts_sub[2,i]))
		res_2 <- as.data.frame(res_2@listData)[,c(2,4,5,6)]
		colnames(res_2)<-c(paste0("Log2fc_",colnames(contrasts_sub)[i]),paste0("Stat_",colnames(contrasts_sub)[i]),paste0("P.value_",colnames(contrasts_sub)[i]),paste0("Adj.p.value_",colnames(contrasts_sub)[i]))
		output_table_2<-cbind(output_table_2,res_2)
		reduced_output_table_2 <- cbind(reduced_output_table_2,res_2)
		rm(res_2)
	}

	## Create annotation table and add gene annotation columns
	# keytype = "ENSEMBL" # will be different if primary annotations are not ENSEMBL, now determined once during non-ercc and kept
	annot <- data.frame(rownames(output_table_2), stringsAsFactors = FALSE)
	colnames(annot)[1]<-keytype
	if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
		annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "SYMBOL", multiVals = "first")
	}
	if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
		annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "GENENAME", multiVals = "first")
	}
	if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
		annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "ENSEMBL", multiVals = "first")
	}
	if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
		annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "REFSEQ", multiVals = "first")
	}
	if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
		annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "ENTREZID", multiVals = "first")
	}

	## Create and add string annotation columns to the annotation table
	string_db <- STRINGdb$new( version="11", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
	string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
	string_map <- string_map[!duplicated(string_map$SYMBOL),]
	annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")

	## Create and add columns containing GOSLIM ids using the panther annotation database to the annotation table
	pthOrganisms(PANTHER.db) <- organism
	panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
	panther <- na.omit(panther)
	annot$GOSLIM_IDS <- panther
	rm(string_db,string_map,panther,keytype)

	## Generate and add all sample mean column to the ERCC-normalized counts table
	output_table_2$All.mean <- rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)
	reduced_output_table_2$All.mean <- rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)

	## Generate and add all sample stdev column to the ERCC-normalized counts table
	output_table_2$All.stdev <- rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)
	reduced_output_table_2$All.stdev <- rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)

	## Add F statistic p-value (similar to ANOVA p-value) column to the ERCC-normalized counts table
	output_table_2$LRT.p.value <- res_2_lrt@listData$padj
	reduced_output_table_2$LRT.p.value <- res_2_lrt@listData$padj

	## Generate and add group mean and stdev columns to the ERCC-normalized counts table
	tcounts <- as.data.frame(t(ERCCnormCounts))
	tcounts$group <- group
	group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
	group_means <- group_means[-c(1),]
	colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))
	group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
	group_stdev <- group_stdev[-c(1),]
	colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))

	output_table_2 <- cbind(output_table_2,group_means)
	reduced_output_table_2 <- cbind(reduced_output_table_2,group_means)

	output_table_2 <- cbind(output_table_2,group_stdev)
	reduced_output_table_2 <- cbind(reduced_output_table_2,group_stdev)

	rm(group_stdev,group_means,tcounts)

	### Add columns needed to generate GeneLab visulaization plots to the ERCC-normalized counts table

	## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
	updown_table <- sign(output_table_2[,grep("Log2fc_",colnames(output_table_2))])
	colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_2),value = TRUE))
	output_table_2 <- cbind(output_table_2,updown_table)
	rm(updown_table)

	## Add column to indicate contrast significance with p <= 0.1
	sig.1_table <- output_table_2[,grep("P.value_",colnames(output_table_2))]<=.1
	colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_2),value = TRUE))
	output_table_2 <- cbind(output_table_2,sig.1_table)
	rm(sig.1_table)

	## Add column to indicate contrast significance with p <= 0.05
	sig.05_table <- output_table_2[,grep("P.value_",colnames(output_table_2))]<=.05
	colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_2),value = TRUE))
	output_table_2 <- cbind(output_table_2,sig.05_table)
	rm(sig.05_table)

	## Add columns for the volcano plot with p-value and adjusted p-value
	log_pval_table <- log2(output_table_2[,grep("P.value_",colnames(output_table_2))])
	colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
	output_table_2 <- cbind(output_table_2,log_pval_table)
	rm(log_pval_table)
	log_adj_pval_table <- log2(output_table_2[,grep("Adj.p.value_",colnames(output_table_2))])
	colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
	output_table_2 <- cbind(output_table_2,log_adj_pval_table)
	rm(log_adj_pval_table)

	### Combine annotations table and the ERCC-normalized counts table

	output_table_2 <- cbind(annot,output_table_2)
	reduced_output_table_2 <- cbind(annot,reduced_output_table_2)
	rownames(output_table_2) <- NULL
	rownames(reduced_output_table_2) <- NULL

	output_table_2$GOSLIM_IDS <- vapply(output_table_2$GOSLIM_IDS, paste, collapse = ", ", character(1L))
	reduced_output_table_2$GOSLIM_IDS <- vapply(reduced_output_table_2$GOSLIM_IDS, paste, collapse = ", ", character(1L))

	### Export human- and computer/visualization- readable DGE tables

	write.csv(output_table_2,file.path(DGE_output_ERCC, "visualization_output_table_ERCCnorm.csv"), row.names = FALSE)
	write.csv(contrasts,file.path(DGE_output_ERCC, "ERCCnorm_contrasts.csv"))
	write.csv(reduced_output_table_2,file.path(DGE_output_ERCC, "ERCCnorm_differential_expression.csv"), row.names = FALSE)

	### Generate and export PCA table for GeneLab visualization plots

	exp_raw <- log2(ERCCnormCounts)
	PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
	write.csv(PCA_raw$x,file.path(DGE_output_ERCC, "visualization_PCA_table_ERCCnorm.csv"), row.names = TRUE)
	rm(exp_raw,PCA_raw)
}

## print session info ##
print("Session Info below: ")
sessionInfo()
## Log same info into versions.txt file """
writeLines(capture.output(sessionInfo()), "versions.txt")
write("BioC_version_associated_with_R_version", file="versions.txt",append=TRUE)
write(toString(tools:::.BioC_version_associated_with_R_version()), file="versions.txt",append=TRUE)
