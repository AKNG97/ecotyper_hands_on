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

qry.rna_NBL <- GDCquery(project = "TARGET-NBL",
                         data.category= "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts")
GDCdownload(qry.rna_NBL)
NBL <- GDCprepare(qry.rna_NBL, summarizedExperiment = TRUE) 

#Get TPM
NBL
NBLL_TPM <- assay(NBL, 4)

colnames(NBLL_TPM) <- gsub("-", "\\.", colnames(NBLL_TPM))
rownames(NBLL_TPM) <- rowData(NBL)$gene_name
length(which(duplicated(rownames(NBLL_TPM))))
dim(NBLL_TPM)
NBLL_TPM <- NBLL_TPM[!duplicated(rownames(NBLL_TPM)),]
dim(NBLL_TPM)
NBLL_TPM <- NBLL_TPM[row.names(NBLL_TPM) %in% annot$HGNC_symbol,]
dim(NBLL_TPM)

NBLL_TPM <- cbind(rownames(NBLL_TPM), NBLL_TPM)
colnames(NBLL_TPM)[1] <- "Gene"
write.table(NBLL_TPM[,1:61], file = "rnas_NBLL_TPM_PC.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)

#Get clinical annot 
clinical_info <- as.data.frame(colData(NBL))

rownames(clinical_info) <- gsub("-", "\\.", rownames(clinical_info))
clinical_info$sample <- gsub("-", "\\.", clinical_info$sample)
clinical_info <- clinical_info %>% mutate(ID=sample)

colnames(clinical_info)

write.table(clinical_info, file = "clinical_info_NBL_tets.txt", row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
