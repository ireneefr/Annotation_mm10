###########################
###   Annotation mm10   ###
###########################

### Download file

# Specify URL where file is stored
url <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/knownGene.txt.gz"
# Specify destination where file should be saved
destfile <- "~/Annotation_mm10/knownGene.txt.gz"
# Apply download.file function
download.file(url, destfile)


### Load file

data <- read.table(gzfile("~/Annotation_mm10/knownGene.txt.gz"), sep = "\t")
head(data)
# Assign colnames based on knownGene.sql file
cols <- c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
          "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID")
colnames(data) <- cols
dim(data) #142446     12
length(unique(data$proteinID)) #53246
length(unique(sub("\\..*", "", data$name))) 
head(data)
# Generate exon start and end columns for 1st exon
vapply(strsplit(data$exonStarts,","), `[`, 1, FUN.VALUE=character(1))

### Obtain gene names

#BiocManager::install("biomaRt")
library(biomaRt)
mart <- useMart("ensembl", "mmusculus_gene_ensembl")
ensemble2gene <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name", "ensembl_gene_id"),
                       filters = "ensembl_transcript_id",
                       values = sub("\\..*", "", data$name), mart = mart)
head(ensemble2gene)
dim(ensemble2gene) #139239      3
length(unique(ensemble2gene$ensembl_gene_id)) #52564
ensemble2gene[ensemble2gene$external_gene_name == "Ptpru",]


