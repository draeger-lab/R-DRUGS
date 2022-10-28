##This script loads a Seurat object containing scRNA-Seq count data, normalizes them by transforming it into TPMs and
##outputs TPMs as well as cell type information in a format that can be loaded by the script to derive core reactions

#Load Dependencies & Dataset
library(data.table)
library(vctrs)
library(Seurat)
#Load scRNA Data
load("PBMC.rda")
df_countMatrix<- GetAssayData(object = All_pre3, slot = "scale.data")

##########################################################################################
# Use human ENSEMBL-gene lengths to calculate FPKM values
df_hg38GeneLengths <- read.table(file = "scRNA/geneLength.v23", 
                                 header = F, sep = "\t", col.names = c("GeneID","Length"), stringsAsFactors = F)
df_hg38GeneLengths$GeneID <- substr(x = df_hg38GeneLengths$GeneID, start = 1, stop = 15)
##########################################################################################
# Set a function to convert count-matrix to TPM-values (columns = samples, rows = gene/transcript IDs)
df_count2TPM <- function(df_countMatrix, df_referenceTranscriptSizes = df_hg38GeneLengths) {
    df_outputMatrix <- df_countMatrix[rownames(df_countMatrix) %in% df_referenceTranscriptSizes[,1],]
  if(nrow(df_outputMatrix) < nrow(df_countMatrix)) {
    print( paste0("We lost ", nrow(df_countMatrix) - nrow(df_outputMatrix), " transcripts due to incompatible transcript/gene ID-matching.") )  
  }
  rownames(df_referenceTranscriptSizes) <- df_referenceTranscriptSizes[,1]
  df_outputMatrix <- df_outputMatrix / df_referenceTranscriptSizes[ rownames(df_outputMatrix), "Length" ]
  df_outputMatrix <- t( t( df_outputMatrix ) / colSums( df_outputMatrix ) )
  df_outputMatrix <- df_outputMatrix * 1000000
  return( df_outputMatrix )}

##########################################################################################
# HGNC ID to Ensembl conversion table
df_HGNC2ENSEMBL <- read.table(file = paste0("/scRNA/mapping_ensembl_hgnc.txt"), sep =",", header = T, stringsAsFactors = F)
colnames(df_HGNC2ENSEMBL) <- c("ENSEMBL","HGNC")
df_HGNC2ENSEMBL <- df_HGNC2ENSEMBL[!duplicated(df_HGNC2ENSEMBL$ENSEMBL) & !duplicated(df_HGNC2ENSEMBL$HGNC),]
# Load ensembl to genesymbol, for human version
df_hg38Ensembl2Genesymbol <- read.csv(file = paste0("scRNA/gencode.v23.annotation_hg38_CHR.gff3_crossref_protein.csv"),
                                      sep = "\t",
                                      header = F,
                                      stringsAsFactors = F)
colnames(df_hg38Ensembl2Genesymbol) <- c("TranscriptID","ProteinID","GeneID","GeneSymbol","GeneSymbolAlt")
df_hg38Ensembl2Genesymbol <- df_hg38Ensembl2Genesymbol[,c(3,4)]
df_hg38Ensembl2Genesymbol$GeneID <- substr(x = df_hg38Ensembl2Genesymbol$GeneID, start = 1, stop = 15)
df_hg38Ensembl2Genesymbol <- df_hg38Ensembl2Genesymbol[!duplicated(df_hg38Ensembl2Genesymbol$GeneSymbol),]
rownames(df_hg38Ensembl2Genesymbol) <- df_hg38Ensembl2Genesymbol$GeneSymbol
#Determine relevant gene sets (HGNC)
lst_str_metabolicGenes <- read.table(file = "scRNA/recon22_genes.csv", header = F, col.names = c("","HGNC"),
                                     skip = 1, sep = ",", row.names = 1, stringsAsFactors = F)[,1]
lst_str_metabolicGenes[!(lst_str_metabolicGenes %in% df_HGNC2ENSEMBL$HGNC)]
lst_str_metabolicGenes <- lst_str_metabolicGenes[lst_str_metabolicGenes != "HGNC:116"]
#convert to ENSEMBL ids
rownames(df_HGNC2ENSEMBL) <- df_HGNC2ENSEMBL$HGNC 
lst_str_metabolicGenes <- df_HGNC2ENSEMBL[lst_str_metabolicGenes, "ENSEMBL"]
# Select only relevant gene-symbols from scRNA data
lst_str_metabolicGenes[!(lst_str_metabolicGenes %in% df_hg38Ensembl2Genesymbol$GeneID)]
# make a list of GeneSymbols we're interested in, instead of a list of HNGC-IDs
lst_str_metabolicGeneSymbols <- unique( df_hg38Ensembl2Genesymbol[df_hg38Ensembl2Genesymbol$GeneID %in% lst_str_metabolicGenes, "GeneSymbol"] )
lst_str_metabolicGeneSymbols[!(lst_str_metabolicGeneSymbols %in% rownames(PBMC))]
# make a copy of Seurat object with reduced gene set
obj_tmp <- PBMC[rownames(PBMC) %in% lst_str_metabolicGeneSymbols,]
obj_tmp
#####An object of class Seurat 
#####1286 features across 8516 samples within 1 assay 
#####Active assay: RNA (1286 features, 163 variable features)
#####2 dimensional reductions calculated: pca, umap

# Convert count data to TPMs
df_tmp <- obj_tmp@assays$RNA@counts
# Make df_tmp into a dataframe and adjust coloumn markers
df_tmp <- data.frame(df_tmp)
colnames(df_tmp) <- gsub(".","-",colnames(df_tmp), fixed=TRUE)

# Check if we have geneLengths for ENSEMBLids?
table(df_hg38Ensembl2Genesymbol$GeneID %in% df_hg38GeneLengths$GeneID)
# Convert Genesymbols to ENSEMBL-IDs (required for TPM calculation)
rownames(df_tmp) <- df_hg38Ensembl2Genesymbol[rownames(df_tmp),"GeneID"]
# Check if the HGNC-ids match?
rownames(df_tmp)[rownames(df_tmp) %in% df_HGNC2ENSEMBL$ENSEMBL]
df_tmp <- df_tmp[rownames(df_tmp)[rownames(df_tmp) %in% df_HGNC2ENSEMBL$ENSEMBL],]
# TPM conversion (function "df_count2TPM()" is stored in setupEnv.R)
df_tmp <- df_count2TPM(df_tmp)
colSums( df_tmp ) # <- should all return 1 million
# Finally convert ENSEMBL gene IDs to HGNC ids
rownames(df_HGNC2ENSEMBL) <- df_HGNC2ENSEMBL$ENSEMBL
rownames(df_tmp) <- df_HGNC2ENSEMBL[rownames(df_tmp),"HGNC"]
# table(is.na(rownames(df_tmp)))
#################### Write files for StanDep ####################


# 1st: Expression data matrix
write.table(t(df_tmp ),"/scRNA/STANDEP/Input/PBMC_expression.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
## simple expression matrix with genes as rows and cells as columns (before transformation via 't()')
## Standep - Condition information
## Standep requires this information for clustering. 
## tissue/ cell type is minimally required and ideally also age of donor
## We convert this to a string and StanDep will group all data points with the same string together
# first check if cells are in similar order
all.equal( colnames(obj_tmp),colnames(df_tmp) )
# TRUE

str_cellData <- paste(obj_tmp@meta.data$CellType, sep=".")
write.table(str_cellData,"/scRNA/STANDEP/Input/PBMC_conditions.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

# Standep - Sample names (i.e. cell names)
# some cell identifiers (barcode or whatever unique string we like)
write.table(colnames(df_tmp),"/scRNA/STANDEP/Input/PBMC_cellnames.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)


## Standep - Gene names
write.table(rownames(df_tmp),"/scRNA/STANDEP/Input/PBMC_genes.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
# should be HGNC ids
## clean up and finish

rm(obj_tmp, df_tmp, str_cellData)

#########################################################
