# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#                                  GLIOTRAIN Data Analysis: Preprocessing
#
# Prepares the datasets as needed in the analyses
#
# Started on the 22th of March, 2021
# 
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


### Environment set-up ----
setwd("/path/to/working_directory/")

library(DESeq2)
library(lumi)
library(lumiHumanIDMapping)
#library(arrayQualityMetrics)
library(sva)
library(org.Hs.eg.db)
library(biomaRt)
# library(readxl)
# library(ade4)
# library(drc)
# library(Bolstad2)
# library(ggplot2)
# library(reshape)

# oooooooooooooooooooooooooooooooooooooooooooooooooo
# I - RNA-Seq ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  # +++ 1) Load data ----
  clinic_data = read.table("./raw/subset1_GTDATAreformatted_clinical_data.tsv", sep = '\t', h = T)
  rnaseq_data_raw = read.csv("../../../GBMdatabase/GTdata/workingFolder/tmp/rnaSeq_raw_resequencedUpdated_200110.csv", header = T, row.names = 1)
  rnaseq_annots = read.table("../../../GBMdatabase/GTdata/workingFolder/tmp/conversion_rnaSeqIDs.txt", h=T, sep='\t')
  EMC_annots = read.csv("../../../GBMAnalysis/EMCcollab/TMZ_Lomustine/Data/tmp/RNA_seq_annotation_Ensembl_BioMart_v89_protein_coding.csv", h=T)
  rnaseq_samplesMap = read.csv("../../../GBMdatabase/GTdata/workingFolder/tmp/samplesLabels_RNASeq.csv", h=T)
  rnaseq_samplesBatch = read.csv("../../../GBMdatabase/GTdata/workingFolder/tmp/samplesBatches_RNASeq.csv", h=T)
  
  
  # +++ 2) Prepare raw counts dataset ----
  colnames(rnaseq_data_raw) = rnaseq_samplesMap$Gliotrain_ID[match(colnames(rnaseq_data_raw), rnaseq_samplesMap$DILA_ID_RNA)]
  rna_data_complete = rnaseq_data_raw[, -grep("CP", colnames(rnaseq_data_raw))]
  colnames(rna_data_complete) = gsub("-\\d+-\\w+-\\d+-\\w+$", "", colnames(rna_data_complete))
  rna_data_complete = rna_data_complete[apply(rna_data_complete, 1, sd) != 0, ]
  
  rnaseq_samplesBatch$sampleid = rnaseq_samplesMap$Gliotrain_ID[match(rnaseq_samplesBatch$sampleid, rnaseq_samplesMap$DILA_ID_RNA)]
  rnaseq_samplesBatch = rnaseq_samplesBatch[, -grep("CP", rnaseq_samplesBatch$sampleid)]
  rnaseq_samplesBatch$sampleid = gsub("-\\d+-\\w+-\\d+-\\w+$", "", rnaseq_samplesBatch$sampleid)
  
  
  # +++ 3) Reduce dataset ----
  ### Remove non-protein-coding genes and combine transcripts that correspond to the same gene by adding them.
  rnaseq_annots$geneID[match(EMC_annots$Ensembl_ID, rnaseq_annots$id)] = EMC_annots$Gene_name # replace gene names in rnaseq_annots with those from EMC
  
  ## Make sure to only keep genes present in the data.
  rnaseq_annots = rnaseq_annots[rnaseq_annots$id %in% rownames(rna_data_complete),]
  EMC_annots = EMC_annots[EMC_annots$Ensembl_ID %in% rownames(rna_data_complete), ]
  rnaseq_annots = rnaseq_annots[match(rownames(rna_data_complete), rnaseq_annots$id),] # put them in the same order as in the dataset
  
  geneNames = unique(EMC_annots$Gene_name) # and determine how many unique coding-protein genes there are.
  rna_data_reduced = as.data.frame(matrix(data = NA, nrow = length(geneNames), ncol = ncol(rna_data_complete),
                                          dimnames = list(geneNames, colnames(rna_data_complete)))) # create the empty dataframe that will receive the "reduced" data
  
  for(g in 1:nrow(rna_data_reduced)){ # for each protein coding gene,
    if(g%%2000 == 0){  print(g)  }
    
    geneData = rna_data_complete[match(EMC_annots$Ensembl_ID[EMC_annots$Gene_name == rownames(rna_data_reduced)[g]], rownames(rna_data_complete)), ] # extract the subset of data corresponding to that gene
    
    # And insert that data in the reduced dataset, adding it per sample if there are more than one transcripts corresponding to that gene.
    if(nrow(geneData) == 1){
      rna_data_reduced[g, ] = geneData
    }else{
      combinedData = sapply(geneData, sum)
      cat(g, " combine\n")
      rna_data_reduced[g, ] = combinedData
    }
    
  }
  rna_data_reduced = rna_data_reduced[apply(rna_data_reduced, 1, sd) != 0, ]
  
  
  # +++ 4) Normalize with vst ----
  ## Assemble relevant cofactors for normalization: source institute, batches, tumour MGMT status.
  normFactors_rna <- data.frame(GTid = colnames(rna_data_complete),
                                center = gsub("-\\d+$", "", colnames(rna_data_complete)),
                                batch = rnaseq_samplesBatch$fcid[match(colnames(rna_data_complete), rnaseq_samplesBatch$sampleid)],
                                MGMT = clinic_data$MGMT_Methylation_Status[match(colnames(rna_data_complete), clinic_data$Subject_ID)])
  normFactors_rna$center = gsub('-', '\\.', normFactors_rna$center)
  
  
  ## Apply normalization.
  # For complete dataset.
  ddsMat <- DESeqDataSetFromMatrix(countData = rna_data_complete, colData = normFactors_rna, design = ~ batch + center + MGMT)
  dim(ddsMat)
  ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 10, ] # filter out transcripts that have less than 10 count total
  dim(ddsMat)
  vstMat <- DESeq2::vst(ddsMat, blind = FALSE)
  rna_norm_complete <- t(assay(vstMat))
  
  # For reduced dataset
  ddsMat <- DESeqDataSetFromMatrix(countData = rna_data_reduced, colData = normFactors_rna, design = ~ batch + center + MGMT)
  dim(ddsMat)
  ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 10, ] # filter out transcripts that have less than 10 count total
  dim(ddsMat)
  vstMat <- DESeq2::vst(ddsMat, blind = FALSE)
  rna_norm_reduced<- t(assay(vstMat))
  
  
  # +++ 5) Save results ----
  save(rna_data_complete, rna_data_reduced, rna_norm_complete, rna_norm_reduced, normFactors_rna, rnaseq_annots, EMC_annots,
       file = "./tmp/analysesDatasets_GT.RNASeq.RData")
  
} # end of I - RNA-Seq



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# II - CNV ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  cnv_data = read.table("./raw/cnv_data_focalEvents.tsv", sep = '\t', h = T)
  wgs2genes = read.table("./raw/logR2genes_map.tsv", sep = '\t', h = T)
  rownames(cnv_data) = cnv_data$region_name
  cnv_data = cnv_data[, -1]
  colnames(cnv_data) = gsub("\\.\\d+\\.\\w+\\.\\d+\\.\\w+\\.\\w+$", "", colnames(cnv_data))
  colnames(cnv_data) = gsub("\\.", "\\-", colnames(cnv_data))
  save(cnv_data, file = "./tmp/analysesDatasets_GT.CNV.RData")
  
} # end of II - CNV



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# III - DASL ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # +++ 1) Load data, transform it to expression data and preprocess ----
  ## Load raw files.
  batch1 <- lumiR.batch(c("./raw/BATCH1_A1222_DataExport_NoNormalization_BATCH1_FFonly.txt"),
                        columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")
  batch2 <- lumiR.batch(c("./raw/BATCH2_A1395_SampleProbeProfile_NoNormalization_BATCH2_FFonly.txt"), 
                        columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")
  batch3 <- lumiR.batch(c("./raw/BATCH3_A1733_Probe_Profile_No_Normalization_Xname_BATCH3_FFonly.txt"), 
                        columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")
  batch4 <- lumiR.batch(c("./raw/BATCH4_A2259_Sample_Probe_Profile_nonorm.txt"),
                        columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")
  
  ## Convert from intensities to expression values and quantile-normalize.
  qc1 = lumiExpresso(batch1, QC.evaluation = T)
  qc2 = lumiExpresso(batch2, QC.evaluation = T)
  qc3 = lumiExpresso(batch3, QC.evaluation = T)
  qc4 = lumiExpresso(batch4, QC.evaluation = T)
  expressions = cbind(exprs(qc1),exprs(qc2),exprs(qc3),exprs(qc4))

  
  ## Filter out probes based on high detection p-value.
  detectionPvalThres = 0.05 # threshold for considering a given expression value reliable
  detectionPvalFreq  = 0.75 # if a probe has p-value above the threshold beyond this frequency in the population, it will be removed
  
  # Identify probesets that should be removed due to high detection p-values in too many samples
  allPvals = cbind(detection(batch1), detection(batch2), detection(batch3), detection(batch4))
  probesToRemove_pvals = apply(allPvals, 1, FUN=function(p) (sum(p>detectionPvalThres)/ncol(allPvals))>=detectionPvalFreq )
  probesToRemove_pvals = names(probesToRemove_pvals[probesToRemove_pvals])
  length(probesToRemove_pvals)
  dasl_DetectPvalFiltered = expressions[ !(rownames(expressions) %in% probesToRemove_pvals), ]
  
  
  ## Correct for batch effect.
  batches = read.table("./raw/batchInfo.txt", h=T)
  batch = batches$Batch
  modcombat = model.matrix(~1, data = batches)
  dasl_batch = ComBat(dat = dasl_DetectPvalFiltered, batch = batch, mod = modcombat)
  
  
  # +++ 2) Reduce dataset ----

  ### Keep only protein-coding genes and combine transcripts that correspond to the same gene by averaging them.
  ## Determine which probes should be kept and averaged.
  # Get names of genes from the dataset.
  probesGenes = read.table("./raw/probesIDandSYMBOL.txt", h=T) # mapping between probesets and corresponding gene
  probesGenes = nuID2EntrezID(nuID=as.character(probesGenes[,2]), lib.mapping="lumiHumanIDMapping", returnAllInfo=TRUE) # mapping between probesets and corresponding gene
  HGNCsymbol = mapIds(org.Hs.eg.db, probesGenes[,2], "SYMBOL", "ENTREZID")
  probesGenes = (cbind(probesGenes, hgnc = as.character(HGNCsymbol)))
  probesGenesLumi = as.data.frame(nuID2IlluminaID(nuID=rownames(probesGenes), lib.mapping="lumiHumanIDMapping", 
                                                  species = "Human", idType = "All"))
  probesGenesMap = cbind(probesGenes, probe.id = probesGenesLumi[match(rownames(probesGenesLumi), rownames(probesGenes)), "Probe_Id"])
  
  # Identify which of these are protein-coding transcripts.
  ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  listAttributes(ensembl)[grep("ncbi", listAttributes(ensembl)[,2], ignore.case = TRUE), ]
  genesMap <- getBM( filters = "illumina_humanht_12_v4", values = probesGenesMap[, "probe.id"], mart = ensembl,
                     attributes = c("hgnc_symbol", "entrezgene_id", "external_gene_name", 
                                    "transcript_biotype", "illumina_humanht_12_v4"))
  genesMap = genesMap[genesMap$transcript_biotype == "protein_coding",]
  probesGenesMap2 = probesGenesMap[probesGenesMap[, "probe.id"] %in% genesMap$illumina_humanht_12_v4,]
  probesGenesMap2 = probesGenesMap2[rownames(probesGenesMap2) %in% rownames(dasl_batch),]
  probesGenesMap2[is.na(probesGenesMap2[, "hgnc"]), "hgnc"] = probesGenesMap2[is.na(probesGenesMap2[, "hgnc"]), "Symbol"]
  
  # Reduce the dataset.
  geneNames = unique(probesGenesMap2[, "hgnc"])
  reduced = as.data.frame(matrix(data = NA, nrow = length(geneNames), ncol = ncol(dasl_batch),
                                 dimnames = list(geneNames, colnames(dasl_batch)))) # create the empty dataframe that will receive the "reduced" data
  
  for(g in 1:nrow(reduced)){ # for each protein coding gene,
    if(g%%2000 == 0){  print(g)  }
    
    # Determine by how many probes the gene is represented.
    probesCount = sum(probesGenesMap2[, "hgnc"] == geneNames[g])
    
    # Insert the data in the reduced dataset, averaging it per sample if there are more than one transcripts corresponding to that gene.
    if(probesCount == 1){
      reduced[g, ] = dasl_batch[match(rownames(probesGenesMap2)[probesGenesMap2[,"hgnc"] == geneNames[g]], 
                                      rownames(dasl_batch)), ]
    }else{
      geneData = as.data.frame(dasl_batch[match(rownames(probesGenesMap2)[probesGenesMap2[,"hgnc"] == geneNames[g]], 
                                                rownames(dasl_batch)), ]) # extract the subset of data corresponding to that gene
      combinedData = sapply(geneData, mean)
      reduced[g, ] = combinedData
    }
    
  }
  dasl_reduced = t(reduced)

  
  # +++ 3) Process the data ----
  
  ## Filter out probes of low variance and mean.
  variances = apply(dasl_reduced, 2, var)
  means = apply(dasl_reduced, 2, mean)
  displayBasicStats(dasl_reduced)
  
  meanThres = summary(means)[2]
  varThres = summary(variances)[2]
  plot(means, variances, pch = 20)
  
  
  # +++ 4) Introduce GS codes and probes names and save ----
  ## Replace the sample names with anonymised names
  samplesGScodes = read.table(file="./raw/sampleGScodes.txt", h=T,sep = '\t')
  sampNames = rownames(dasl_reduced)
  anoNames = c()
  for(i in 1:length(sampNames)){
    anoNames = c(anoNames, paste0("GS",samplesGScodes[which(samplesGScodes[,1]==sampNames[i]),2]))
  }
  rownames(dasl_reduced) = anoNames
  
  ## Get the dasl unbiased dataset.
  load("./tmp/transcFiltered_unbiasedApproach.RData")
  dasl_complete = transcFiltered_unbiasedApproach
  
  # Save.
  save(dasl_complete, dasl_reduced, probesGenesMap2, file="./tmp/analysesDatasets_EMC.DASL.RData")
  
  
  
} # end of III - DASL