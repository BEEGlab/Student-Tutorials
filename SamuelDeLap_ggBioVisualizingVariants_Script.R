#Installation instructions have been included with each library
#Installation only needs to be done once
#Accession needs to be done in every run

#Access necessary libraries
#install.packages(BiocManager)
library(BiocManager)
#install.packages(ggplot2)
library(ggplot2)
#BiocManager::install(ggbio)
library(ggbio)
#BiocManager::install(biomaRt)
library(biomaRt)
#BiocManager::install(bioVizBase)
library(biovizBase)
#BiocManager::install(Homo.sapiens)
library(Homo.sapiens)
#BiocManager::install(VariantAnnotation)
library(VariantAnnotation)

#Perform a biomaRt search to determine desired genes
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
GeneSelection <- getBM(attributes = "external_gene_name", filters = "phenotype_description", values = "Breast Cancer", mart = human)

#Make Annotation Graphs for Genes we want
data(genesymbol, package = "biovizBase")
wh <- genesymbol["BRCA1"]
wh <- range(wh, ignore.strand = TRUE)
p.BRCA1 <- autoplot(Homo.sapiens, which = wh)
p.BRCA1

wh <- genesymbol["ATM"]
wh <- range(wh, ignore.strand = TRUE)
p.ATM <- autoplot(Homo.sapiens, which = wh)
p.ATM

#Import the .vcf file containing variant information
fl.vcf <- system.file("extdata", "17-1409-CEU-brca1.vcf.bgz", package="biovizBase")
vcf <- readVcf(fl.vcf, "hg19")

#Coerce the vcf file into a VRanges option
vr <- as(vcf[, 1:3], "VRanges")

#Rename the columns of the VRanges object so that they match what we have in wh
vr <- renameSeqlevels(vr, value = c("17" = "chr17"))

#Use autoplot on the VRanges object
p.vr <- autoplot(vr, which = wh)

#Display the autoplotted function
p.vr

#Create a genomic range (GRanges) that we can use with autoplot to zoom in
gr17 <- GRanges("chr17", IRanges(41234400, 41234530))

#Zoom in on the desired genomic range
p.vr + xlim(gr17)

#Zoom in even more
p.vr + xlim(gr17) + zoom()

#Using a custom geometry
autoplot(vr, which = wh, geom = "rect", arrow = FALSE)