########################################
###   Obtain mm10 annotation table   ###
########################################

### Read annotation transcript
ann_transcript <- read.csv("annotation_transcript.csv")
dim(ann_transcript) #139239     17
head(ann_transcript)
ann_transcript <- ann_transcript[order(ann_transcript$chrom),]

### Read annotation CpGs
ann_cpgs <- read.csv("MouseMethylation-12v1-0_A2.csv", skip = 7)
dim(ann_cpgs) #287693     38
head(ann_cpgs)
ann_cpgs <- ann_cpgs[order(ann_cpgs$CHR),]

### Function to annotate CpG given its position
annotate_cg <- function(pos_cg, transcript){
  group <- c()
  name <- c()
  gene <- c()
  for(i in 1:nrow(transcript)){
    if(transcript[i, "strand"] == "+"){
      if(pos_cg %in% seq(transcript[i, "txStart"]-200, transcript[i, "txStart"])){
        group <- append(group, "TSS200")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "txStart"]-1500, transcript[i, "txStart"]-201)){
        group <- append(group, "TSS1500")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "txStart"]+1, transcript[i, "cdsStart"]-1)){
        group <- append(group, "5'UTR")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "cdsEnd"]+1, transcript[i, "txEnd"]-1)){
        group <- append(group, "3'UTR")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "cdsStart"], transcript[i, "end_1stExon"])){
        group <- append(group, "1stExon")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "end_1stExon"]+1, transcript[i, "cdsEnd"])){
        group <- append(group, "Body")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      }
    } else if(transcript[i, "strand"] == "-"){
      if(pos_cg %in% seq(transcript[i, "txEnd"], transcript[i, "txEnd"]+200)){
        group <- append(group, "TSS200")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "txEnd"]+201, transcript[i, "txEnd"]+1500)){
        group <- append(group, "TSS1500")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "cdsEnd"]+1, transcript[i, "txEnd"]-1)){
        group <- append(group, "5'UTR")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "txStart"]+1, transcript[i, "cdsStart"]-1)){
        group <- append(group, "3'UTR")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "end_1stExon"], transcript[i, "cdsEnd"])){
        group <- append(group, "1stExon")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      } else if(pos_cg %in% seq(transcript[i, "cdsStart"], transcript[i, "end_1stExon"]-1)){
        group <- append(group, "Body")
        name <- append(name, transcript[i, "id_transcript"])
        gene <- append(gene, transcript[i, "external_gene_name"])
      }
    }
  }
  return(list(UCSC_Transcript = name, UCSC_Group = group, UCSC_Gene = gene))
}

#change MT chromosome name by M to match with annotation transcript
ann_cpgs[ann_cpgs$CHR == "MT", "CHR"] <- "M"

### Assign number of cores to parallelize
#n_cores <- detectCores()
#registerDoParallel(cores = n_cores-1)

dim(ann_cpgs)
ann_cpgs_sub <- ann_cpgs[!(ann_cpgs$CHR %in% c("", "0")),]
dim(ann_cpgs_sub)

### Run annotate_cg function for each CpG
start <- Sys.time()
for(i in 1:nrow(ann_cpgs_sub)){
    res <- annotate_cg(ann_cpgs_sub[i, "MAPINFO"], ann_transcript[ann_transcript$chrom == paste0("chr", ann_cpgs_sub[i, "CHR"]),])
    ann_cpgs_sub[i, "UCSC_Transcipt"] <- paste(res$UCSC_Transcript, collapse = ";")
    ann_cpgs_sub[i, "UCSC_Group"] <- paste(res$UCSC_Group, collapse = ";")
    ann_cpgs_sub[i, "UCSC_Gene"] <- paste(res$UCSC_Gene, collapse = ";")
}
end <- Sys.time()
time_diff <- end-start

### Save table
write.csv(ann_cpgs_sub, "~/Desktop/Annotation_mm10/mm10_cpg_annotation.csv", row.names = FALSE)
