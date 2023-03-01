#### Load libraries ####

library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")

features <- c("ensembl_gene_id", "chromosome_name", 
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs, 
               mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)
dim(annot)
head(annot)
annot[annot == ""] = NA
annot <- annot[!is.na(annot$HGNC_symbol), ]
dim(annot)
annot <- annot[!duplicated(annot$HGNC_symbol), ]
dim(annot)

#### Get TCGA data ####

qry.rna_COAD <- GDCquery(project = "TCGA-COAD",
                          data.category= "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")
GDCdownload(qry.rna_COAD)
COAD <- GDCprepare(qry.rna_COAD, summarizedExperiment = TRUE) 

#Get TPM
COAD
COADL_TPM <- assay(COAD, 4)

colnames(COADL_TPM) <- gsub("-", "\\.", colnames(COADL_TPM))
rownames(COADL_TPM) <- rowData(COAD)$gene_name
length(which(duplicated(rownames(COADL_TPM))))
dim(COADL_TPM)
COADL_TPM <- COADL_TPM[!duplicated(rownames(COADL_TPM)),]
dim(COADL_TPM)
COADL_TPM <- COADL_TPM[row.names(COADL_TPM) %in% annot$HGNC_symbol,]
dim(COADL_TPM)

COADL_TPM <- cbind(rownames(COADL_TPM), COADL_TPM)
colnames(COADL_TPM)[1] <- "Gene"
write.table(COADL_TPM, file = "rnas_COADL_TPM_PC.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

#Get clinical annot 
clinical_info <- as.data.frame(colData(COAD))

rownames(clinical_info) <- gsub("-", "\\.", rownames(clinical_info))
clinical_info$sample <- gsub("-", "\\.", clinical_info$sample)
clinical_info <- clinical_info %>% mutate(ID=sample)

colnames(clinical_info)

write.table(clinical_info, file = "clinical_info_COAD.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
