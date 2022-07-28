# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#                   Erasmus Data Preprocessing and Extraction
#
# Process the data from EMC to have it ready for subsequent analyses. That includes:
#   * For DASL data:
#     - QC and normalisation of the raw data (==> "./tmp/transc_preprocessed.RData")
#     - Extraction of probes for the unbiased analysis (==> "./tmp/transcFiltered_unbiasedApproach.RData")
#     - Extraction of probes for the cancer pathways/genes analysis (==> "./tmp/transcFiltered_cancerGenes.RData")
#     - Extraction of probes for the ALDH/Retinoic Acid pathway analysis (==> "./tmp/transcFiltered_aldhGenes.RData")
#   * For AUC data: ("./tmp/ic50_datasets.RData")
#     - QC of raw viability data
#     - Fitting log-logistic, exponential decay or linear models to the data
#     - Calculating normalized AUCs from the models
#   * Defining lists of samples that: (==> "./tmp/samplesLists.RData")
#     - contain all samples present in the DASL Data
#     - contain all GBM samples present in the DASL Data
#     - contain all Primary GBM samples present in the DASL Data
#     - contain all Recurrent GBM samples present in the DASL Data
#     - contain all samples present in the IC50 Dataset
#     - contain all Primary samples present in the DASL Data
#     - contain all Recurrent samples present in the DASL Data
#
# Started on the 12th of March, 2020. Copied from previous version
# 
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

setwd("/path/to/working_directory/")

library(lumi)
library(lumiHumanIDMapping)
library(arrayQualityMetrics)
library(sva)
library(org.Hs.eg.db)
library(biomaRt)
library(readxl)
library(ade4)
library(drc)
library(Bolstad2)
library(ggplot2)
library(reshape)


# oooooooooooooooooooooooooooooooooooooooooooooooooo
# I - Define samples lists ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

{
  diagnosis_raw = read.csv("./tmp/DASL_88samples_PA_diagnosis.csv", h=T, stringsAsFactors = F) # derived from "./raw/DASL_88 samples_PA_diagnosis_10072019_IN (2).xsls" ; "./raw/OS rates_GBM_EMC-220717_IN.xls" and "./raw/sampleGScodes.txt"
  diagnosis = diagnosis_raw[, -1]
  rownames(diagnosis) = diagnosis_raw[, 1]
  head(diagnosis)

  ### Lists with all samples ----
  samples_all88            = diagnosis_raw[, 1]
  samples_allGBM           = diagnosis_raw[ diagnosis_raw$Pathological_diagnosis=="GBM" , 1]
  samples_allPrimaryGBM    = diagnosis_raw[ diagnosis_raw$Pathological_diagnosis=="GBM" & diagnosis_raw$Recurrent=="No" , 1]
  samples_allRecurrentGBM  = diagnosis_raw[ diagnosis_raw$Pathological_diagnosis=="GBM" & diagnosis_raw$Recurrent=="Yes" , 1]
  sapply(list(samples_all88, samples_allGBM, samples_allPrimaryGBM, samples_allRecurrentGBM), length)
  
  ## Long-term and short-term samples
  samples_allSTsurvivors = diagnosis_raw[ diagnosis_raw$OS < 9  & !is.na(diagnosis_raw$OS), 1 ]
  samples_allLTsurvivors = diagnosis_raw[ diagnosis_raw$OS > 36 & !is.na(diagnosis_raw$OS), 1 ]
  sapply(list(samples_allSTsurvivors, samples_allLTsurvivors), length)
  
  
  ### Lists with samples that were used for the drug screening ----
  samples_drugExposedPrimary   = c("GS436",  "GS437", "GS452",  "GS454", "GS304", "GS149", "GS102",  "GS356",  "GS343p", "GS245",
                                   "GS261",  "GS289", "GS104p", "GS257", "GS284", "GS184", "GS79",   "GS203",  "GS295",  "GS274",	
                                   "GS357",  "GS281", "GS365",  "GS359", "GS330", "GS475", "GS280p", "GS413",  "GS461",  "GS463", 
                                   "GS254",  "GS423", "GS464")
  samples_drugExposedRecurrent = c("GS279c", "GS249", "GS224",  "GS186", "GS324", "GS446", "GS286",  "GS453",  "GS401",  "GS216c", "GS323p", "GS335")
  samples_drugExposedAll = c(samples_drugExposedPrimary, samples_drugExposedRecurrent)
  samples_drugExposedAll = samples_drugExposedAll[ sample.int(length(samples_drugExposedAll)) ] # randomise primary and recurrent within the list
  sapply(list(samples_drugExposedAll, samples_drugExposedPrimary, samples_drugExposedRecurrent), length)
  sum(samples_drugExposedAll %in% samples_all88) # = 45 if manually inputed samples all correspond to the ones in the diagnosis table
  
  ### Save all lists ----
  save(samples_all88, samples_allGBM, samples_allPrimaryGBM, samples_allRecurrentGBM, 
       samples_allSTsurvivors, samples_allLTsurvivors,
       samples_drugExposedPrimary, samples_drugExposedRecurrent, samples_drugExposedAll,
       diagnosis,
       file = "./tmp/samplesLists.RData")
}



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# II - Prepare the DASL Data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

# ===== II.1 Define gene symbols table =====
# Create a table that has
#   - the probes Id as rownames, 
#   - one column for the corresponding gene HGNC symbol, 
#   - one column for the names that will be attributed to the probe in the format geneSymbol.X where x=1,2... since a gene can have several corresponding probes
{
  probesGenes = read.table("./raw/probesIDandSYMBOL.txt", h=T) # mapping between probesets and corresponding gene
  probesGenes = (nuID2EntrezID(nuID=as.character(probesGenes[,2]), lib.mapping="lumiHumanIDMapping", returnAllInfo=TRUE)) # mapping between probesets and corresponding gene
  HGNCsymbol = mapIds(org.Hs.eg.db, probesGenes[,2], "SYMBOL", "ENTREZID")
  probesGenes = (cbind(probesGenes, as.character(HGNCsymbol)))
  head(probesGenes)
  typeof(probesGenes)
  counter = 1 # keep track of how many times a given gene has appeared
  probes.names = c()
  
  ## Define numbered gene names 
  probes.names = c( probes.names, paste(probesGenes[1,3], counter,sep='.') ) # symbol is used instead of HGNCsymbol because there are NAs in HGNCsymbol column
  for(probe in 2:nrow(probesGenes)){
    
    if(probesGenes[probe, 3] != probesGenes[(probe-1), 3]){ counter = 1 }
    probes.names = c( probes.names, paste(probesGenes[probe, 3], counter,sep='.') )
    counter = counter + 1
    
  }
  
  probesGenesMap = as.data.frame(cbind(probesGenes, probes.names))
  colnames(probesGenesMap) = c("accession", "entrezID", "probe.symbol", "probe.HGNCsymbol", "probe.name")
  head(probesGenesMap)

}


# ===== II.2 Processing =====

## Import the raw intensity files ----
batch1 <- lumiR.batch(c("./raw/BATCH1_A1222_DataExport_NoNormalization_BATCH1_FFonly.txt"),
                      columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")

batch2 <- lumiR.batch(c("./raw/BATCH2_A1395_SampleProbeProfile_NoNormalization_BATCH2_FFonly.txt"), 
                      columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")

batch3 <- lumiR.batch(c("./raw/BATCH3_A1733_Probe_Profile_No_Normalization_Xname_BATCH3_FFonly.txt"), 
                      columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")

batch4 <- lumiR.batch(c("./raw/BATCH4_A2259_Sample_Probe_Profile_nonorm.txt"),
                      columnNameGrepPattern=list(beadNum=NA), inputAnnotation=T, lib.mapping='lumiHumanIDMapping', dec=".")


## Overview of the raw data ----
# Get a sense of characteristic distributions
{
  # Expression data
  batchesExprs = list(exprs(batch1), exprs(batch2), exprs(batch3), exprs(batch4)) # expression data
  batchesMeans = lapply(batchesExprs, function(b) apply(b,1,mean) )
  batchesVar   = lapply(batchesExprs, function(b) apply(b,1,var) )
  par(mfrow=c(2,2))
  for(b in 1:length(batchesExprs)){  plot(batchesMeans[[b]], batchesVar[[b]], pch=20, xlab="Mean", ylab="Variance", main=paste0("Expression for batch ",b))  }
  par(mfrow=c(1,1))
  
  # Detection p-values
  batchesDpVal = list(detection(batch1), detection(batch2), detection(batch3), detection(batch4)) # detection p-values
  par(mfrow=c(2,2))
  for(b in 1:length(batchesExprs)){  
    pvalsBelow1 = t(batchesDpVal[[b]])[t(batchesDpVal[[b]])<=1]
    hist(pvalsBelow1, xlab = "p-value", main=paste0("Detection p-values for batch ",b), breaks=seq(0,1,0.01), xlim=c(0,1))
    abline(v=0.05, col="red", lty=2)
  }
  par(mfrow=c(1,1))
  lapply(batchesDpVal, function(b) sum(b>1) )
  lapply(batchesDpVal, function(b) (b[b>1]) )
  lapply(batchesDpVal, function(b) sum(b>0.05) )
}

# PCA
{
  batchesExprsComb = t(cbind(exprs(batch1), exprs(batch2), exprs(batch3), exprs(batch4))) # expression data
  batches = read.table("./raw/batchInfo.txt", h=T)
  samplesGScodes = read.table(file="./raw/sampleGScodes.txt", h=T,sep = '\t')
  sampNames = rownames(batchesExprsComb)
  anoNames = c()
  for(i in 1:length(sampNames)){
    anoNames = c(anoNames, paste0("GS",samplesGScodes[which(samplesGScodes[,1]==sampNames[i]),2]))
  }
  rownames(batchesExprsComb) = anoNames
  bInfo = batches$Batch
  sampNames = rownames(batches)
  anoNames = c()
  for(i in 1:length(sampNames)){
    anoNames = c(anoNames, paste0("GS",samplesGScodes[which(samplesGScodes[,1]==sampNames[i]),2]))
  }
  names(bInfo) = anoNames
  
  # Run PCA
  acp_raw = dudi.pca(batchesExprsComb, scannf = FALSE, nf = 10)
  pve <- 100*acp_raw$eig/sum(acp_raw$eig)
  cumsum(pve)
  par(mfrow=c(2,2))
  s.label(acp_raw$li, xax=1, yax=2)
  s.label(acp_raw$li, xax=3, yax=4)
  s.label(acp_raw$li, xax=5, yax=6)
  s.label(acp_raw$li, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for batches
  gcol = c("red","blue","green","brown")
  par(mfrow=c(2,2))
  s.class(dfxy = acp_raw$li, fac=as.factor(bInfo), col=gcol, xax=1, yax=2)
  s.class(dfxy = acp_raw$li, fac=as.factor(bInfo), col=gcol, xax=3, yax=4)
  s.class(dfxy = acp_raw$li, fac=as.factor(bInfo), col=gcol, xax=5, yax=6)
  s.class(dfxy = acp_raw$li, fac=as.factor(bInfo), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for Prim/Rec samples
  par(mfrow=c(2,2))
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=1, yax=2)
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=3, yax=4)
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=5, yax=6)
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for Pathological diagnosis
  gcol = rainbow(length(levels(as.factor(diagnosis$Pathological_diagnosis))))
  par(mfrow=c(2,2))
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=1, yax=2)
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=3, yax=4)
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=5, yax=6)
  s.class(dfxy = acp_raw$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
}


## Convert from intensities to expression values and quantile-normalize ----
qc1 = lumiExpresso(batch1, QC.evaluation = T)
qc2 = lumiExpresso(batch2, QC.evaluation = T)
qc3 = lumiExpresso(batch3, QC.evaluation = T)
qc4 = lumiExpresso(batch4, QC.evaluation = T)
expressions = cbind(exprs(qc1),exprs(qc2),exprs(qc3),exprs(qc4))

# Plot Var against mean
qcExpr = list(exprs(qc1),exprs(qc2),exprs(qc3),exprs(qc4))
qcMeans = lapply(qcExpr, function(b) apply(b,1,mean) )
qcVar   = lapply(qcExpr, function(b) apply(b,1,var) )
par(mfrow=c(2,2))
for(b in 1:length(qcExpr)){  plot(qcMeans[[b]], qcVar[[b]], pch=20, xlab="Mean", ylab="Variance", main=paste0("Expression for normalized batch ",b))  }
par(mfrow=c(1,1))


## Correct for batch effect ----
batches = read.table("./raw/batchInfo.txt", h=T)
batch = batches$Batch
modcombat = model.matrix(~1, data=batches)
DASL_Processed = ComBat(dat=expressions, batch=batch, mod=modcombat)

# See differences
{
  batchCorrected = t(DASL_Processed)
  sampNames = rownames(batchCorrected)
  anoNames = c()
  for(i in 1:length(sampNames)){
    anoNames = c(anoNames, paste0("GS",samplesGScodes[which(samplesGScodes[,1]==sampNames[i]),2]))
  }
  rownames(batchCorrected) = anoNames
  bInfo = batches$Batch
  sampNames = rownames(batches)
  anoNames = c()
  for(i in 1:length(sampNames)){
    anoNames = c(anoNames, paste0("GS",samplesGScodes[which(samplesGScodes[,1]==sampNames[i]),2]))
  }
  names(bInfo) = anoNames
  
  batchExprs_corrected = lapply(1:4, function(d) (batchCorrected[ batch==d, ]) )
  correctedMeans = lapply(batchExprs_corrected, function(b) apply(b,2,mean) )
  correctedVar   = lapply(batchExprs_corrected, function(b) apply(b,2,var) )
  par(mfrow=c(2,2))
  for(b in 1:length(batchExprs_corrected)){  plot(correctedMeans[[b]], correctedVar[[b]], 
                                                  pch=20, xlab="Mean", ylab="Variance", main=paste0("Expression for corrected batch ",b))  }
  par(mfrow=c(1,1))
  allBatches_mean = apply(batchCorrected,2,mean)
  allBatches_var  = apply(batchCorrected,2,var)
  plot(allBatches_mean, allBatches_var, pch=20, xlab="Mean", ylab="Variance", main="Expression for all data\nafter batch correction")
  
  acp_Corrected = dudi.pca((batchCorrected), scannf = FALSE, nf = 10)
  pve <- 100*acp_Corrected$eig/sum(acp_Corrected$eig)
  cumsum(pve)
  
  par(mfrow=c(2,2))
  s.label(acp_Corrected$li, xax=1, yax=2)
  s.label(acp_Corrected$li, xax=3, yax=4)
  s.label(acp_Corrected$li, xax=5, yax=6)
  s.label(acp_Corrected$li, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for batches
  gcol = c("red","blue","green","brown")
  par(mfrow=c(2,2))
  s.class(dfxy = acp_Corrected$li, fac=as.factor(bInfo), col=gcol, xax=1, yax=2)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(bInfo), col=gcol, xax=3, yax=4)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(bInfo), col=gcol, xax=5, yax=6)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(bInfo), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for Prim/Rec samples
  par(mfrow=c(2,2))
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=1, yax=2)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=3, yax=4)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=5, yax=6)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Recurrent), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for Pathological diagnosis
  gcol = rainbow(length(levels(as.factor(diagnosis$Pathological_diagnosis))))
  par(mfrow=c(2,2))
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=1, yax=2)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=3, yax=4)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=5, yax=6)
  s.class(dfxy = acp_Corrected$li, fac=as.factor(diagnosis$Pathological_diagnosis), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
}


## Filter out probes based on high detection p-value ----
{
  detectionPvalThres = 0.05 # threshold for considering a given expression value reliable
  detectionPvalFreq  = 0.75 # if a probe has p-value above the threshold beyond this frequency in the population, it will be removed
  
  # Identify probesets that should be removed due to high detection p-values in too many samples
  allPvals = cbind(detection(batch1), detection(batch2), detection(batch3), detection(batch4))
  probesToRemove_pvals = apply(allPvals, 1, FUN=function(p) (sum(p>detectionPvalThres)/ncol(allPvals))>=detectionPvalFreq )
  probesToRemove_pvals = names(probesToRemove_pvals[probesToRemove_pvals])
  length(probesToRemove_pvals)
  DASL_Processed_DetectPvalFiltered = DASL_Processed[ !(rownames(DASL_Processed) %in% probesToRemove_pvals), ]
}

## Introduce GS codes and probes names and save ----
{
  # Replace the sample names with anonymised names
  samplesGScodes = read.table(file="./raw/sampleGScodes.txt", h=T,sep = '\t')
  sampNames = colnames(DASL_Processed_DetectPvalFiltered)
  anoNames = c()
  for(i in 1:length(sampNames)){
    anoNames = c(anoNames, paste0("GS",samplesGScodes[which(samplesGScodes[,1]==sampNames[i]),2]))
  }
  colnames(DASL_Processed_DetectPvalFiltered) = anoNames
  
  # Replace probe ID with numbered gene names
  rownames(DASL_Processed_DetectPvalFiltered) = probesGenesMap$probe.name[ match(rownames(DASL_Processed_DetectPvalFiltered), 
                                                                                 rownames(probesGenesMap)) ]
}
# Save.
transc_preprocessed = t(DASL_Processed_DetectPvalFiltered)
dim(transc_preprocessed)
transc_preprocessed[1:5, 1:5]
save(transc_preprocessed, probesGenesMap, file="./tmp/transc_preprocessed.RData")


# ===== II.3 Select probes for unbiased approach =====

## Filter out probes based on redundancy, and low mean and variance ----
{
  # For genes that have several probes, select the probe with the highest variance
  probesGenes = probesGenesMap[ probesGenesMap$probe.name %in% colnames(transc_preprocessed), ] # account for the bad pvalue-filtered probes
  genesDuplicated = unique(probesGenes$probe.symbol[duplicated(probesGenes$probe.symbol)]) # genes that have several corresponding probes
  probesToKeep = as.character(probesGenes$probe.name[!(probesGenes$probe.symbol %in% genesDuplicated)]) # genes that don't have several corresponding probes
  selectedProbes = c()
  
  for(g in 1:length(genesDuplicated)){
    subset = transc_preprocessed[, probesGenes$probe.symbol==genesDuplicated[g]]
    variances = apply(subset,2,var)
    selectedProbes = c(selectedProbes, colnames(subset)[variances==max(variances)])
  }
  
  probesToKeep_genesDuplicates = c(probesToKeep, selectedProbes)
  
  transc_genesDuplicates = transc_preprocessed[, colnames(transc_preprocessed) %in% probesToKeep_genesDuplicates]
  dim(transc_genesDuplicates)
  
  
  # Filter out low expression signal and low variance probes
  variances = apply(transc_genesDuplicates, 2, var)
  means = apply(transc_genesDuplicates, 2, mean)
  
  meanThres = summary(means)[2]
  varThres = summary(variances)[2]
  plot(means, variances, pch = 20)
  abline(v=meanThres, h=varThres, lty=2, col="red")
  dim(transc_genesDuplicates[, variances>varThres])
  dim(transc_genesDuplicates[, means>meanThres])
  dim(transc_genesDuplicates[, ((variances>varThres) & (means>meanThres))])
  transcFiltered_unbiasedApproach = transc_genesDuplicates[, ((variances>varThres) & (means>meanThres)) ] # filter out probes with variance or mean below the 1st quartile of their respective distribution
  
  # Adapt probesGenesMap dataset
  probesGenesMap_unbiasedApproach = probesGenesMap[match(colnames(transcFiltered_unbiasedApproach), probesGenesMap$probe.name), ]
  probesGenesMap_unbiasedApproach$probe.name = probesGenesMap_unbiasedApproach$probe.symbol # the "gene.X" format is unnecessary here since we removed duplicate genes
  dim(probesGenesMap_unbiasedApproach)
  head(probesGenesMap_unbiasedApproach)
  
  # Change the gene names in the transc dataset
  colnames(transcFiltered_unbiasedApproach) = gsub("\\.\\d+$", "", colnames(transcFiltered_unbiasedApproach)) # the "gene.X" format is unnecessary here since we removed duplicate genes
  dim(transcFiltered_unbiasedApproach)
  transcFiltered_unbiasedApproach[1:5, 1:5]
  
  
}

## Save ----
save(transcFiltered_unbiasedApproach, meanThres, varThres,
     probesGenesMap_unbiasedApproach,
     file = "./tmp/transcFiltered_unbiasedApproach.RData")


# ===== II.4 Select probes for approach focusing on cancer genes =====

## Import lists of cancer genes ----
{
  agilentCGenes = read.table("./raw/cancerGenes_Agilent.txt", h=T)
  ampliSeqCGenes = read.table("./raw/cancerGenes_IlluminaAmpliSeq.txt", h=T)
  truSightCGenes = read.table("./raw/cancerGenes_IlluminaTrusight.txt", h=T)
  NSPatwhaysPanelCGenes = read.table("./raw/cancerGenes_NanoStringPathwaysPanel.txt", sep='\t', h=T)
  NSProgressionPanelCGenes = read.table("./raw/cancerGenes_NanoStringProgressionPanel.txt", sep='\t', h=T)
  qiagen1CGenes = read.table("./raw/cancerGenes_Qiagen1.txt", h=T)
  qiagen2CGenes = read.table("./raw/cancerGenes_Qiagen2.txt", sep=',', h=T)
  
  cancerGenesLists = list(NSPatwhaysPanelCGenes, NSProgressionPanelCGenes, # from nanoString nCounter Pancancer Panels
                          agilentCGenes, # from Agilent microarray 
                          ampliSeqCGenes, truSightCGenes, # from Illumina microarrays
                          qiagen1CGenes, qiagen2CGenes) # from Qiagen microarrays
  lapply(cancerGenesLists, dim)
  lapply(cancerGenesLists, colnames)
}

### Extract cancer genes from DASL data ----
## Probes for which the corresponding gene is in at least one of the cancer genes lists are identified and extracted.
## A matrix indicating in which genes list(s) each probe was found is also generated.
{
  daslInCGenesLists = matrix(data=NA, nrow=nrow(probesGenesMap), ncol=length(cancerGenesLists)) # mapping which probe gene was found in which list
  rownames(daslInCGenesLists) = probesGenesMap$probe.name
  colnames(daslInCGenesLists) = c("NSPatwhaysPanelCGenes", "NSProgressionPanelCGenes", 
                                  "Agilent1", "Illumina.AmpliSeq", "Illumina.TruSight", "Qiagen1", "Qiagen2")
  genesListInDasl = list(c(), c(), c(), c(), c(), c(), c())
  names(genesListInDasl) = c("NSPatwhaysPanelCGenes", "NSProgressionPanelCGenes", 
                             "Agilent", "AmpliSeq", "TruSight", "Qiagen1", "Qiagen2")
  
  probesNames = rownames(daslInCGenesLists)
  for(probe in 1:length(probesNames)){ # for each of the 30,000 probes
    # Test if the corresponding gene is in each of the cancer genes lists
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], 
                       c(as.character(NSPatwhaysPanelCGenes$Official.Symbol), as.character(NSPatwhaysPanelCGenes$Alias...Prev.Symbol)))
    ) > 0) |
    (length(grep( probesGenesMap$probe.symbol[probe], 
                  c(as.character(NSPatwhaysPanelCGenes$Official.Symbol), as.character(NSPatwhaysPanelCGenes$Alias...Prev.Symbol)))
    ) > 0)  ){ 
      daslInCGenesLists[probe, 1] = "Yes" 
      genesListInDasl[[1]] = c(genesListInDasl[[1]], probesGenesMap$probe.symbol[probe])
    }
    
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], 
                       c(as.character(NSProgressionPanelCGenes$Gene.Names), as.character(NSProgressionPanelCGenes$Alias...Prev.Symbol), as.character(NSProgressionPanelCGenes$HUGO.name.if.different)))
    ) > 0) |
    (length(grep( probesGenesMap$probe.symbol[probe], 
                  c(NSProgressionPanelCGenes$Gene.Names, NSProgressionPanelCGenes$Alias...Prev.Symbol, NSProgressionPanelCGenes$HUGO.name.if.different))
    ) > 0)){ 
      daslInCGenesLists[probe, 2] = "Yes" 
      genesListInDasl[[2]] = c(genesListInDasl[[2]], probesGenesMap$probe.symbol[probe])
    }
    
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], agilentCGenes[,1])) > 0) |
         (length(grep( probesGenesMap$probe.symbol[probe], agilentCGenes[,1])) > 0)){ 
      daslInCGenesLists[probe, 3] = "Yes" 
      genesListInDasl[[3]] = c(genesListInDasl[[3]], probesGenesMap$probe.symbol[probe])
    }
    
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], ampliSeqCGenes[,1])) > 0) |
         (length(grep( probesGenesMap$probe.symbol[probe], ampliSeqCGenes[,1])) > 0)){ 
      daslInCGenesLists[probe, 4] = "Yes" 
      genesListInDasl[[4]] = c(genesListInDasl[[4]], probesGenesMap$probe.symbol[probe])
    }
    
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], truSightCGenes[,1])) > 0) |
         (length(grep( probesGenesMap$probe.symbol[probe], truSightCGenes[,1])) > 0)){ 
      daslInCGenesLists[probe, 5] = "Yes" 
      genesListInDasl[[5]] = c(genesListInDasl[[5]], probesGenesMap$probe.symbol[probe])
    }
    
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], qiagen1CGenes[,1])) > 0) |
         (length(grep( probesGenesMap$probe.symbol[probe], qiagen1CGenes[,1])) > 0)){ 
      daslInCGenesLists[probe, 6] = "Yes" 
      genesListInDasl[[6]] = c(genesListInDasl[[6]], probesGenesMap$probe.symbol[probe])
    }
    
    if(  (length(grep( probesGenesMap$probe.HGNCsymbol[probe], c(as.character(qiagen2CGenes[,1]), as.character(qiagen2CGenes[,2])))) > 0) |
         (length(grep( probesGenesMap$probe.symbol[probe], c(as.character(qiagen2CGenes[,1]), as.character(qiagen2CGenes[,2])))))  ){ 
      daslInCGenesLists[probe, 7] = "Yes" 
      genesListInDasl[[7]] = c(genesListInDasl[[7]], probesGenesMap$probe.symbol[probe])
    }
  }
  
  apply(daslInCGenesLists,2, function(d) sum(!is.na(d))) # number of probes from DASL data that are in the cancer genes list
  sapply(genesListInDasl, function(l) length(unique(l))) # number of genes from the cancer genes lists that could be found in the DASL data
  sapply(cancerGenesLists, dim) # the previous 2 lines give much higher number because of all the duplicates in the genes and the fact that current, previous and HGNC symbols are screened to find a correspondance
  table(apply(daslInCGenesLists,1, function(p) sum(!is.na(p)))) # for a given probe/gene, count in how many cancer gene list it appears
  sum(table(apply(daslInCGenesLists,1, function(p) sum(!is.na(p))))[2:8]) # count the number of probes that appear in at least one cancer genes list
  
  
  ## Extract probes that correspond to a gene present in at least one of the cancer genes lists
  daslInCGenesLists = daslInCGenesLists[ apply(daslInCGenesLists,1, function(p) sum(!is.na(p))>0), ]
  dim(daslInCGenesLists)
  head(daslInCGenesLists)
  
  daslInCGenesLists_inTranscPreprocessed = rownames(daslInCGenesLists) %in% colnames(transc_preprocessed) # some of the identified daslInCGenesLists were previously removed from the dasl datasets due to high detection p-values
  transcFiltered_cancerGenes = transc_preprocessed[ , match(rownames(daslInCGenesLists)[daslInCGenesLists_inTranscPreprocessed], 
                                                            colnames(transc_preprocessed)) ]
  dim(transcFiltered_cancerGenes)
  transcFiltered_cancerGenes[1:5, 1:5]
  sum(is.na(transcFiltered_cancerGenes))
}

## Adapt probesGenesMap
probesGenesMap_cancerGenes = probesGenesMap[match(colnames(transcFiltered_cancerGenes), probesGenesMap$probe.name), ]

## Save the datasets
write.csv2(daslInCGenesLists, file = "./tmp/DASLprobes_in_CancerGenesLists.csv", row.names=T)
save(transcFiltered_cancerGenes, probesGenesMap_cancerGenes, file = "./tmp/transcFiltered_cancerGenes.RData")


# ===== II.5 Select probes for approach focusing on ALDH/Retinoic Acid Pathway =====
{
  GoI = c("ALDH1A1", "ALDH1A3", "SOX2", "OLIG2", "GFAP") # genes of interest
  GoIexpr = lapply(GoI, function(g) transc_preprocessed[, which(probesGenesMap$probe.HGNCsymbol==g | probesGenesMap$probe.symbol==g) ])
  names(GoIexpr) = GoI
  sapply(GoIexpr, dim)
  sapply(GoIexpr, length)
  
  # Define the corresponding probes into a dataframe
  transcFiltered_aldhGenes = cbind(GoIexpr$ALDH1A1, GoIexpr$ALDH1A3, GoIexpr$SOX2, GoIexpr$OLIG2, GoIexpr$GFAP)
  colnames(transcFiltered_aldhGenes) = c("ALDH1A1.1", "ALDH1A1.2", "ALDH1A3.1", "ALDH1A3.2", "SOX2.1", "SOX2.2", "OLIG2", "GFAP")
  dim(transcFiltered_aldhGenes)
  head(transcFiltered_aldhGenes)
}

## Save
save(transcFiltered_aldhGenes, file = "./tmp/transcFiltered_aldhGenes.RData")


# ===== II.6 Select probes corresponding to protein-coding genes =====
{
  ## Get data.
  load("./tmp/transc_preprocessed.RData")
  probesGenesMap_preprocessed = probesGenesMap[probesGenesMap$probe.name %in% colnames(transc_preprocessed),]
  dim(probesGenesMap_preprocessed)
  probesGenesLumi = as.data.frame(nuID2IlluminaID(nuID=rownames(probesGenesMap_preprocessed), lib.mapping="lumiHumanIDMapping", 
                                                  species = "Human", idType = "All"))
  probesGenesMap2 = cbind(probesGenesMap_preprocessed, probe.id = probesGenesLumi[match(rownames(probesGenesLumi), rownames(probesGenesMap_preprocessed)), "Probe_Id"])

  ## Identify which probes correspond to protein-coding transcripts.
  ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  listAttributes(ensembl)[grep("ncbi", listAttributes(ensembl)[,2], ignore.case = TRUE), ]
  genesMap <- getBM( filters = "illumina_humanht_12_v4", values = probesGenesMap2$probe.id, mart = ensembl,
                     attributes = c("hgnc_symbol", "entrezgene_id", "external_gene_name", 
                                    "transcript_biotype", "illumina_humanht_12_v4"))
                                    
  table(genesMap$transcript_biotype)
  head(genesMap)
  dim(genesMap)
  dim(probesGenesMap2)
  head(probesGenesMap2)
  
  genesMap2 = genesMap[genesMap$transcript_biotype == "protein_coding", ]
  dim(genesMap2)
  sum(duplicated(genesMap2$illumina_humanht_12_v4))
  head(genesMap2)
  genesMap2[duplicated(genesMap2$illumina_humanht_12_v4),]
  
  transc_proteinCoding = transc_preprocessed[, colnames(transc_preprocessed) %in% probesGenesMap2$probe.name[match(unique(genesMap2$illumina_humanht_12_v4), 
                                                                                                                   probesGenesMap2$probe.id)]]
  dim(transc_proteinCoding)
  
  ## Combine duplicates.
  head(probesGenesMap2)
  genesPositions = probesGenesMap2[match(colnames(transc_proteinCoding), probesGenesMap2$probe.name), c("probe.symbol", "probe.name")]
  transcFiltered_proteinCoding = as.data.frame(matrix(data = NA, nrow = nrow(transc_proteinCoding), ncol = length(unique(genesPositions$probe.symbol)),
                                                      dimnames = list(rownames(transc_proteinCoding), unique(genesPositions$probe.symbol))))
  for(g in 1:ncol(transcFiltered_proteinCoding)){
    geneData = as.data.frame(transc_proteinCoding[, colnames(transc_proteinCoding) %in% genesPositions$probe.name[genesPositions$probe.symbol == colnames(transcFiltered_proteinCoding)[g]]])
#    cat(g, "\n", geneData, "\n")
    transcFiltered_proteinCoding[, g] = apply(geneData, 1, mean)
  }
  transcFiltered_proteinCoding[1:5,1:5]
  save(transcFiltered_proteinCoding, genesMap2, probesGenesMap2, file="./tmp/transcFiltered_proteinCoding.RData")
  
}


# oooooooooooooooooooooooooooooooooooooooooooooooooo
# III - Prepare the AUC Values ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

# ===== III.1 - Import data =====
# "./tmp/viabilityIC50s_transformation2.xlsx" is derived from "./raw/20191021_rawdata2.xls" for the content of the data; and "./raw/20191015_Fu data2.xls" for the mapping between EMC and GS codes

## Raw data points
viability_raw = as.data.frame(read_excel("./tmp/viabilityIC50s_transformation2.xlsx", na="NA", sheet = 1, col_names = T))
rownames(viability_raw) = viability_raw[,1]
viability_raw = viability_raw[,2:ncol(viability_raw)]
viability_raw[1:10, 1:10]

## The mapping of which dose range was tested for each drug/sample combination
drugs_range = as.data.frame(read_excel("./tmp/viabilityIC50s_transformation2.xlsx", na="NA", sheet = 3, col_names = T))
rownames(drugs_range) = drugs_range[,1]
drugs_range = drugs_range[, 2:ncol(drugs_range)]
drugs_range

## The ranges of doses that were tested in the experiments
tested_ranges = as.data.frame(read_excel("./tmp/viabilityIC50s_transformation2.xlsx", na="NA", sheet = 2, col_names = T))
tested_ranges[1:5, 1:5]



# ===== III.2 - Transform into per drug datasets =====
# Will create a list, each element of which contains the data for a given drug
# The dataset will have columns for all possible concentrations, some of which will be mostly empty since 
# a given drug was not tested with all these doses;
# In the rows there will be the duplicates for each sample.

concentrations_chr = c("160 ?M", "80 ?M",	"16 ?M",	"8 ?M",	"1.6 ?M",	"0.8 ?M",	"0.16 ?M",	
                       "0.08 ?M",	"0.016 ?M",	"0.008 ?M",	"0.0016 ?M", "0.0008 ?M", "0.00008 ?M"	) # will be used as column names for the datasets
concentrations_num = c(160, 80, 16, 8, 1.6, 0.8, 0.16, 0.08, 0.016, 0.008, 0.0016, 0.0008, 0.00008)

drugs_position = gsub("_\\w+", "", rownames(viability_raw)) # to determine the positions of the duplicates in the raw data
samples_position = gsub("_\\w+", "", colnames(viability_raw)) # to determine the position of all tested doses for a given sample in the raw data
samples_list = unique(samples_position) # list the samples in the raw data, without duplicates

## Compute means(DMSO) for each sample
# will be used to compute %viabilities(sample, drug) = rawdata(sample, drug) / means(DMSO)
# where rawdata has the datapoints for a given duplicate
dmsoValues = matrix(data = NA, nrow = length(samples_list), ncol = 6,
                    dimnames = list(samples_list, c("con1", "con2", "con3", "con4", "con5", "con6")))
for(samp in samples_list){
  dmsoSamp = viability_raw[drugs_position=="DMSO", samples_position==samp]
  dmsoValues[samp,] = sapply(dmsoSamp, mean, na.rm=T)
}

## Prepare the rownames for the datasets
rowNames_viabsPerDrug = c()
for(samp in samples_list){  rowNames_viabsPerDrug = c(rowNames_viabsPerDrug, paste0(samp, "_dup1"), paste0(samp, "_dup2"))  }

## Create and fill the list
viabs_perDrug = list()
for(drug in unique(drugs_position)[2:length(unique(drugs_position))]){ # for each tested drug (position 1 is DMSO)
  
  data.drug = matrix(data = NA, nrow = length(rowNames_viabsPerDrug), ncol = length(concentrations_chr),
                     dimnames = list(rowNames_viabsPerDrug, concentrations_chr)) 

  samples_to_remove = c() # will be filled with names of samples that have no data for this drug, and removed from the dataset later
  for(samp in samples_list){ # for each sample
    vals <- viability_raw[drugs_position == drug, samples_position == samp] # extract the data for this sample and this drug
    
    if(drugs_range[drug , samp] == 0 | sum(is.na(vals)) == length(as.matrix(vals))){  samples_to_remove = c(samples_to_remove, samp)  } # if the sample was identified as to be removed in tested_ranges or if it has no data, add it to the to_remove list
    
    else{
      
      ## Determine which cells in the new dataset should be filled
      dScale = drugs_range[drug , samp] # identify which of the dose ranges was tested
      concentrations_idx = match(tested_ranges[dScale,2:7], colnames(data.drug)) # determine columns to be filled
      samples_idx = which(gsub("_\\w+","",rowNames_viabsPerDrug) == samp) # determine rows to be filled
      
      ## Compute and fill % viability values
      vals = sapply(1:6, function(val) vals[, val]/dmsoValues[samp, val]) # %viability(drug, sample, dose) = rawdata(drug, sample, dose) / means(DMSO, sample, dose)
      for(val in 1:6){  data.drug[samples_idx, concentrations_idx[val]] = vals[, val]  }

    }
  }
  
  ## Remove samples that have no data for this drug from the dataset
  data.drug = data.drug[!(gsub("_\\w+", "", rowNames_viabsPerDrug) %in% samples_to_remove), ]
  
  ## Add the new dataset to the list
  viabs_perDrug[[drug]] = data.drug
  
}



# ===== III.3 - Filter out obvious bad drugs and samples =====
# Here, drugs for which almost all samples seem to be unaffected (% viability doesn't vary much with dose increase)
# as well as drugs for which that effect is ambiguous are tested to determine whether it makes sense to include them
# in the subsequent analysis, or if the drug indeed has no effect and should be filtered out.
# These drugs are identified by looking at the curves before any fitting.
{
  
  ## Drugs that seem to not really affect the cells
  noEffectDrugsNames = c("Allopurinol", "Altretamine", "Anastrozole", "Capecitabine", "Cyclophosphamide", "Ifosfamide",
                         "Lenalidomide", "Letrozole", "Pentostatin", "Pomalidomide", "Thalidomide")
  noEffectDrugs = lapply(noEffectDrugsNames, function(drug) viabs_perDrug[[drug]])
  names(noEffectDrugs) = noEffectDrugsNames
  
  
  ## Drugs for which the effect is ambiguous
  ambigEffectDrugsNames = c("Cabazitaxel", "Cisplatin", "Dactinomycin", "Docetaxel", "Fulvestrant", "Methotrexate",
                            "Nelarabine", "Paclitaxel", "Pemetrexed", "Pralatrexate", "Procarbazine hydrochloride",
                            "Streptozocin", "Temozolomide", "Vinblastine sulfate", "Vincristine sulfate")
  ambigEffectDrugs = lapply(ambigEffectDrugsNames, function(drug) viabs_perDrug[[drug]])
  names(ambigEffectDrugs) = ambigEffectDrugsNames
  
  
  ## Function to fit a linear model to the data and extract the slope p-value
  getLMslopePvalue <- function(data, concentrations){
    # @params <data> is the dataset having the tested doses as columns and response as rows (with distinct lines per duplicate)
    # @params <concentrations> the vector of tested doses present in <data>
    # @return <LMslopePvalues> the vector of p-values
    
    sampNames = gsub("_\\w+","",rownames(data)) # extract only the sample name, without the "_dup1/2" part, to identify the samples
    LMslopePvalues = c()
    
    for(samp in unique(sampNames)){ # for each sample
      
      data.samp = data[gsub("_\\w+","",rownames(data)) == samp,] # extract the datapoints for this sample
      data.sample = c(data.samp[1,], data.samp[2,]) # put all datapoints in a vector rather than in a matrix with duplicates as rows
      
      samp.lm <- lm(data.sample ~ log(rep(concentrations,2)), na.action = na.omit)
      LMslopePvalues = c(LMslopePvalues, coefficients(summary(samp.lm))[2,4])
      
    }
    
    #print(unique(sampNames))
    names(LMslopePvalues) = unique(sampNames)
    
    return(LMslopePvalues)
    
  }
  
  ## Function to fit a linear model to the data and extract the slope coefficient
  getLMslopeCoef <- function(data, concentrations){
    # @params <data> is the dataset having the tested doses as columns and response as rows (with distinct lines per duplicate)
    # @params <concentrations> the vector of tested doses present in <data>
    # @return <LMslopeCoef> the vector of slope values
    
    sampNames = gsub("_\\w+","",rownames(data)) # extract only the sample name, without the "_dup1/2" part, to identify the samples
    LMslopeCoef = c()
    
    for(samp in unique(sampNames)){ # for each sample
      
      data.samp = data[gsub("_\\w+","",rownames(data)) == samp,] # extract the datapoints for this sample
      data.sample = c(data.samp[1,], data.samp[2,]) # put all datapoints in a vector rather than in a matrix with duplicates as rows
      
      samp.lm <- lm(data.sample ~ log(rep(concentrations,2)), na.action = na.omit)
      LMslopeCoef = c(LMslopeCoef, coefficients(samp.lm)[2])
      
    }
    
    #print(unique(sampNames))
    names(LMslopeCoef) = unique(sampNames)
    
    return(LMslopeCoef)
    
  }
  
  noEffectDrugs.pvalues = sapply(noEffectDrugs, function(drug) getLMslopePvalue(drug, concentrations_num))
  apply(noEffectDrugs.pvalues,2, function(drug) sum(drug<0.05, na.rm=T)/length(drug))
  boxplot(noEffectDrugs.pvalues)
  abline(h=0.05, lty=2, col="red")
  
  noEffectDrugs.coefs = sapply(noEffectDrugs, function(drug) getLMslopeCoef(drug, concentrations_num))
  apply(noEffectDrugs.coefs,2, function(drug) sum(drug<(0.0), na.rm=T)/length(drug))
  boxplot(noEffectDrugs.coefs)
  abline(h=0.0, lty=2, col="red")
  
  sapply(noEffectDrugsNames, 
         function(drug) sum(sapply(rownames(noEffectDrugs.coefs), 
                                   function(samp) (noEffectDrugs.coefs[samp, drug] < 0 & noEffectDrugs.pvalues[samp, drug] < 0.05  ))) 
  )
  
  
  ambigEffectDrugs.pvalues = sapply(ambigEffectDrugs, function(drug) getLMslopePvalue(drug, concentrations_num))
  sapply(ambigEffectDrugs.pvalues, function(drug) sum(drug<0.05, na.rm=T)/length(drug))
  boxplot(ambigEffectDrugs.pvalues)
  abline(h=0.05, lty=2, col="red")
  
  ambigEffectDrugs.coefs = sapply(ambigEffectDrugs, function(drug) getLMslopeCoef(drug, concentrations_num))
  sapply(ambigEffectDrugs.coefs, function(drug) sum(drug<(0.0), na.rm=T)/length(drug))
  boxplot(ambigEffectDrugs.coefs)
  abline(h=0.0, lty=2, col="red")
  
  sapply(ambigEffectDrugsNames, 
         function(drug) sum(sapply(names(ambigEffectDrugs.coefs[[drug]]), 
                                   function(samp) (ambigEffectDrugs.coefs[[drug]][samp] < 0 & ambigEffectDrugs.pvalues[[drug]][samp] < 0.05  ))) 
  )
  
  ### Remove drugs with low potential for analysis.
  # The drugs for which dose effect is amibiguous will be kept and samples
  # excluded on a case by case basis.
  noEffectDrugsNames = c(noEffectDrugsNames, "Cabazitaxel", "Docetaxel")
  length(noEffectDrugsNames)
  length(viabs_perDrug)
  viabs_perDrug_filtered = viabs_perDrug[which(!(names(viabs_perDrug) %in% noEffectDrugsNames))]
  length(viabs_perDrug_filtered)
  
} # end of III.3 - Filter out unusable drugs



# ===== III.4 - Fit response models =====
# Here, for each (drug, sample) combination:
#    * a sigmoid, an exponential decay and a linear model will be tested against the data
#    * if it succeeds new datapoints will be generated to produce a smooth curve of most
#         appropriate model, and the values for the coefficients of the model will be extracted
#    * if it fails for all models, the sample will be excluded from the analysis for the drug
#    * either way, the standard deviation from the raw datapoints will also be computed
{
  
## Define the new "dose range" of points for the smooth curve
newR = seq(min(concentrations_num), max(concentrations_num), length.out=500) # create 500 points evenly spaced between the lowest and highest tested doses
newRange = sort(c(newR, concentrations_num[2:(length(concentrations_num)-1)])) # also insert the initial tested doses
#newRange = log(newR) # the sigmoid curve is presented on a log scale


### Define functions used to fit the models.

## Function to fit reponse models to the data and extract relevant information
fitResponseModels <- function(data, concentrations, plotRange){
  # @params <data> is the dataset having the tested doses as columns and response as rows (with distinct lines per duplicate)
  # @params <concentrations> the vector of tested doses present in <data>
  # @params <plotRange> the vector of points that will be used to compute the new, smooth curve; these points should already be on the logged scale
  # @return <responseData> a list containing: 
  #                        - the <fittedValues> dataset for the smooth curves, where each row corresponds to a given sample, and each column to a given dose from <plotRange>
  #                        - the <metaData> dataset for models coefficients for each successfully fitted sample
  #                        - the <excludedSamples> vector listing all samples that were not properly modelled for this drug
  #                        - the <standardDeviations> dataset givind the SD from the duplicates for each of the originally tested doses
  
  responseModels = list() # will store each best successfully fitted response model
  modelsTypes = c() # will store which response model was selected for each sample
  sampNames = gsub("_\\w+","",rownames(data)) # extract only the sample name, without the "_dup1/2" part, to identify the samples
  
  for(samp in unique(sampNames)){ # for each sample
    
    cat(samp,", ")
    data.samp = data[gsub("_\\w+","",rownames(data)) == samp,] # extract the datapoints for this sample
    data.sample = c(data.samp[1,], data.samp[2,]) # put all datapoints in a vector rather than in a matrix with duplicates as rows
    respModels = list() # to store models that have been fitted and compare them
    
    
    ### Try to fit each of the 3 models to the data
    # Try to fit a sigmoid model to the data
    
    tryCatch(
      expr = {  respModels[["sigmoid"]] = drm(data.sample ~ rep(concentrations, 2), fct = LL.3())  }, 
      error = function(e){ 
        cat('\n')
        message(e) 
        cat('\n')
      } # if it fails, do nothing particular
    )
    
    # Try to fit an exponential decay model to the data
    tryCatch(
      expr = {  respModels[["exponential"]] = drm(data.sample[!is.na(data.sample)] ~ ( log(rep(concentrations, 2)) + rep(-log(min(concentrations)), 2*length(concentrations)) )[!is.na(data.sample)], fct = EXD.2())  },
      error = function(e){ 
        cat('\n')
        message(e) 
        cat('\n')
      } # if it fails, do nothing particular
    )
    
    # Try to fit a linear model to the data
    tryCatch(
      expr = {  respModels[["linear"]] = lm(data.sample ~ ( log(rep(concentrations, 2)) + rep(-log(min(concentrations)), 2*length(concentrations)) ), na.action = na.omit)  },
      error = function(e){ 
        cat('\n')
        message(e) 
        cat('\n')
      } # if it fails, do nothing particular
    )
    
    
    ### Chose the most appropriate model.
    if(length(respModels) > 0){ # if at least one model could be fitted
      
      ## Identify the model with the lowest residual standard error
      RSEs = c(  ifelse("sigmoid" %in% names(respModels), (summary(respModels[["sigmoid"]]))$rseMat[1], NA),
                 ifelse("exponential" %in% names(respModels), (summary(respModels[["exponential"]]))$rseMat[1], NA),
                 ifelse("linear" %in% names(respModels), (summary(respModels[["linear"]]))$sigma, NA)
      ) # extract the residual standard error from the models
      names(RSEs) = c("sigmoid", "exponential", "linear")
      chosenModel = c("sigmoid", "exponential", "linear")[which(RSEs == min(RSEs, na.rm = T))]
      #cat(RSEs, chosenModel, "\n")
      
      # When the exponential model gets too high at small dose, prefer the sigmoid model when available
      if(chosenModel == "exponential"){
        if(coef(respModels[["exponential"]])[1] > 3 & !is.na(RSEs[1])){
          if(coef(respModels[["sigmoid"]])[1] < 3){  chosenModel = "sigmoid"  }
        }
      }
      
      responseModels[[samp]] = respModels[[chosenModel]]
      modelsTypes = c(modelsTypes, chosenModel)
      
      
    } # end of if(length(respModels) > 0)
    
    
    ## Regardless of model fitting, compute standard deviation for the original data
    #stanDevs[samp, ] = apply(data.samp, 2, sd, na.rm=T) # compute the standard deviations from the raw datapoints
    
  } # end of for(samp in unique(sampNames))
  cat('\n')
  names(modelsTypes) = names(responseModels)
  
  ## Extract information from successful model
  fittedValues <- as.data.frame(t(sapply(names(responseModels), function(samp) computeSmoothCurve(responseModels[[samp]], modelsTypes[samp], plotRange, -log(min(concentrations))))))
  colnames(fittedValues) = as.character(log(plotRange))
  metaData <- as.data.frame(t(sapply(names(responseModels), function(samp) extractMetaData(responseModels[[samp]], modelsTypes[samp], -log(min(concentrations)))))) # model coefficients
  metaData[, 1:11] = sapply(1:11, function(i) as.numeric(as.character(metaData[,i]))) # collected coefficients are factors so they must be converted to numeric
  
  ## Remove poor models
  # Start with samples that could not be fit to any model
  nonFit = unique(sampNames)[!(unique(sampNames) %in% rownames(metaData))] 
  excludedSamples = rep("Unfittable", length(nonFit))
  names(excludedSamples) = nonFit
  
  # Then add samples for which the viability is monotonously increasing with dose
  for(samp in rownames(metaData)){ 
    if(metaData[samp, "model_type"] == "linear" & metaData[samp, "steepness_estimate"] > 0){  excludedSamples[samp] = "Model increasing with dose"  }
    else if(metaData[samp, "model_type"] %in% c("sigmoid","exponential") & metaData[samp, "steepness_estimate"] < 0){  excludedSamples[samp] = "Model increasing with dose"  }
  }
  
  # Then add samples for which the upper limit is above an arbitrary cut-off of 3
  for(samp in rownames(metaData)){ 
    if(metaData[samp, "upperLimit_estimate"] > 3.0 & !(samp %in% excludedSamples)){  excludedSamples[samp] = "Inaccurate model"  }
  }
  
  # Remove identified samples from the datasets
  metaData = metaData[!(rownames(metaData) %in% names(excludedSamples)), ]
  fittedValues = fittedValues[match(rownames(metaData), rownames(fittedValues)), ]
  
  
  ## Compile and terminate the function
  responseData = list(fittedValues, metaData, excludedSamples)#, stanDevs)
  names(responseData) = c("fittedValues", "metaData", "excludedSamples")#, "standardDeviations")
  
  return(responseData)
  
}

### Function to compute the smooth curve based on selected fitted model
computeSmoothCurve <- function(model, modelType, plotRange, offset){
  # @params <model> the selected fit model
  # @params <modelType> one of "sigmoid", "exponential", "linear", to determine the expression to use
  # @params <plotRange> the data points for which the values need to be calculated
  # @params <offset> the x offset for the exponential model
  # @returns <curveValues> the calculated values corresponding to <plotRange>
  
  if(modelType == "sigmoid"){
    curveValues = (coef(model)[2])/( 1 + exp((coef(model)[1])*(log(plotRange) - log((coef(model)[3])))))
  }else if(modelType == "exponential"){
    curveValues = (coef(model)[1]) * (exp((-(log(plotRange) + offset))/coef(model)[2]))
  }else if(modelType == "linear"){
    curveValues = (coef(model)[1]) + ((log(plotRange) + offset) * coef(model)[2])
    curveValues[curveValues < 0] = 0.0
  }
  
  return(curveValues)
  
}

### Function to extract the coefficients information from the fitted sigmoid model
extractMetaData <- function(model, modelType, offset){
  # !!! The exponential model was computed on the log(dose) scale with an offset. 
  #     The output coefficients correspond to that model, except for the EC50 estimate that 
  #     has been transformed back to the dose scale without the offset
  # @params <model> the response model
  # @params <modelType> one of "sigmoid", "exponential", "linear"
  # return <metaData> the dataset containing estimated value, std error and p-value for the coefficients <steepness> (or slope), 
  #                       the <upperLimit> (or intercept) and <ec50> for the model and corresponding curve, 
  #                       as well as residual standard error, corresponding degrees of freedom and type of the model
  
  modelSummary = summary(model)
  modelCoefs = modelSummary$coefficients
  
  if(modelType == "sigmoid"){
    metaData = c(modelCoefs[1,1], modelCoefs[1,2], modelCoefs[1,4], # b parameter is steepness
                 modelCoefs[2,1], modelCoefs[2,2], modelCoefs[2,4], # d parameter is upper limit
                 modelCoefs[3,1], modelCoefs[3,2], modelCoefs[3,4], # e parameter is EC50
                 modelSummary$rseMat[1], modelSummary$rseMat[2],
                 modelType)
  }else if(modelType == "exponential"){
    EC50 <- ED(model, respLev=c(50)) # EC50 has to be computed
    metaData = c(modelCoefs[2,1], modelCoefs[2,2], modelCoefs[2,4], # e parameter is steepness
                 modelCoefs[1,1], modelCoefs[1,2], modelCoefs[1,4], # d parameter is upper limit
                 exp(EC50[1]-offset), EC50[2], NA, 
                 modelSummary$rseMat[1], modelSummary$rseMat[2],
                 modelType)
  }else if(modelType == "linear"){
    metaData = c(modelCoefs[2,1], modelCoefs[2,2], modelCoefs[2,4], # slope of the linear regression
                 modelCoefs[1,1], modelCoefs[1,2], modelCoefs[1,4], # intercept of the linear regression
                 exp(((-0.5 * modelCoefs[1,1])/modelCoefs[2,1]) - offset), NA, NA, # EC50 has to be computed as exp[(((intercept/2) - intercept) / (slope)) - offset]
                 modelSummary$sigma, modelSummary$df[2],
                 modelType)
  }
  
  names(metaData) = c("steepness_estimate", "steepness_stdError", "steepness_pvalue",
                      "upperLimit_estimate", "upperLimit_stdError", "upperLimit_pvalue",
                      "EC50_estimate", "EC50_stdError", "EC50_pvalue",
                      "residual_stdErr", "on_DegreesOfFreedom",
                      "model_type")
  return(metaData)
}


### Actually run the functions
fittedData <- lapply(viabs_perDrug_filtered, function(drug) fitResponseModels(drug, concentrations_num, newRange))
names(fittedData) = names(viabs_perDrug_filtered)


# Count how many samples were not successfully modelled for each drug
sapply(fittedData, function(drug) length(drug$excludedSamples))

# Count how many of each model type were used
sapply(fittedData, function(drug) table(drug$metaData$model_type))

# Look at the distribution of Residual Standard Error
boxplot(sapply(fittedData, function(drug) (drug[["metaData"]][, "residual_stdErr"])))
abline(h=0.4, lty=2, col='red')

## Determine samples with a high Residual standard error
count_highResStdErr = (sapply(fittedData, function(drug) sum(drug[["metaData"]][, "residual_stdErr"]>0.4)))
highResStdErr = data.frame()
for(drug in names(count_highResStdErr[count_highResStdErr>0])){
  metaData_highResStdErr = fittedData[[drug]][["metaData"]][fittedData[[drug]][["metaData"]][, "residual_stdErr"]>0.4,]
  highResStdErr = rbind(highResStdErr, cbind("drug"=rep(drug, nrow(metaData_highResStdErr)),metaData_highResStdErr))
}
colnames(highResStdErr) = c("drug","steepVal","steepErr","steepPval","upVal","upErr","upPval","ec50Val","ec50Err","ec50Pval", "ResStdErr","DF","model")
highResStdErr

} # end of III.4 - Fit response models



# ===== III.5 - Print curves  =====

## Compute the standard deviations for each data point
rawData_standardDeviations = list()
for(drug in names(viabs_perDrug)){
  
  rawData = viabs_perDrug[[drug]]
  sampNames = gsub("_\\w+","",rownames(rawData)) # extract only the sample name, without the "_dup1/2" part, to identify the samples
  stanDevs = matrix(data = NA, nrow = length(unique(sampNames)), ncol = length(concentrations_num)) # will store the standard deviations from the raw data
  rownames(stanDevs) = unique(sampNames)
  
  for(samp in unique(sampNames)){ # for each sample
    data.samp = rawData[sampNames == samp,] # extract the datapoints for this sample
    stanDevs[samp, ] = apply(data.samp, 2, sd, na.rm=T) # compute sd
  }
  
  rawData_standardDeviations[[drug]] = stanDevs
  
}
# Determine the highest SD across all drugs and samples
maxSD = max( sapply(rawData_standardDeviations, function(drug) max(drug, na.rm=T)), na.rm = T )

## Plot the lines from the original data
rawData_curves = list()
for(drug in names(viabs_perDrug)){ # for each drug
  print(drug)
  
  # Extract relevant data for this drug
  rawData = viabs_perDrug[[drug]]
  SDs = rawData_standardDeviations[[drug]]
  colnames(SDs) = log(concentrations_num)
  samples = rownames(rawData)
  
  ## Compute the mean for a given sample
  meansData = matrix(data = NA, nrow = length(unique(gsub("_\\w+","",samples))), ncol = length(concentrations_num)) # create empty matrix
  rownames(meansData) = unique(gsub("_\\w+","",samples))
  for(samp in unique(gsub("_\\w+","",samples))){ # for each sample 
    
    data.sample = rawData[gsub("_\\w+","",samples) == samp,] # extract the datapoints for the duplicates
    meansData[samp, ] = apply(data.sample, 2, mean, na.rm=T) # compute the mean for the duplicate and put it in the matrix
  }
  
  ## Extract mean and sd into one dataset
  drug.data = as.data.frame(matrix(data = NA, nrow = sum(!is.na(SDs)), ncol = 4))
  colnames(drug.data) = c("sample", "x", "mean", "sd")
  idx = 1
  for(s in 1:nrow(SDs)){
    for(con in 1:ncol(SDs)){
      if(!is.na(SDs[s, con])){
        drug.data[idx, 1] = samples[s]
        drug.data[idx, 2] = log(concentrations_num[con])
        drug.data[idx, 3] = meansData[s, con]
        drug.data[idx, 4] = SDs[s, con]
        idx = idx + 1
      }
    }
  }
  rawData_curves[[drug]] = drug.data
  
  ## Plot the graph showing all curves for the excluded samples for this drug
  png(paste0("./results/AUCs_production/beforeFitting/", drug, "_curves.png"))
  print(
    ggplot( data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="lightgrey")) +
      labs(title = drug, x = "log(drug)", y = "% viability", colour = "sd") +
      geom_line(data = drug.data, aes(x = x, y = mean, group = sample)) +
      scale_colour_gradient2(limits=c(0, maxSD), low="blue", mid="white", high="red", midpoint = maxSD/2) +
      geom_point(data = as.data.frame(drug.data), mapping = aes(x=x, y=mean, color = sd))
  )
  dev.off()
  
}


## Plot the lines for samples that were successfully fitted
fittedData_curves = list()
for(drug in names(fittedData)){
  print(drug)
  
  ## Extract relevant datasets
  data = as.data.frame(fittedData[[drug]]["fittedValues"])
  colnames(data) = log(newRange)
  metaData = as.data.frame(fittedData[[drug]]["metaData"])
  SDs = rawData_standardDeviations[[drug]]
  SDs = SDs[!(rownames(SDs) %in% rownames(as.data.frame(fittedData[[drug]]["excludedSamples"]))),]
  colnames(SDs) = log(concentrations_num)
  
  ## Transform the standard deviations data
  stanDevs = matrix(data = NA, nrow = sum(!is.na(SDs)), ncol = 3)
  colnames(stanDevs) = c("x", "y", "sd")
  stanDevs_idx = 1
  for(s in 1:nrow(SDs)){
    for(con in 1:ncol(SDs)){
      if(!is.na(SDs[s, con])){
        stanDevs[stanDevs_idx, 1] = log(concentrations_num[con])
        stanDevs[stanDevs_idx, 2] = data[match(rownames(SDs)[s],rownames(data)), match(colnames(SDs)[con], colnames(data))]
        stanDevs[stanDevs_idx, 3] = SDs[s, con]
        stanDevs_idx = stanDevs_idx + 1
      }
    }
  }
  
  sampNames = rownames(data)
  modelType = metaData[match(sampNames, rownames(metaData)), 12]
  names(modelType) = sampNames
  data = cbind(sampNames, data, modelType)
  data_melted = melt(data, id.vars = c("sampNames","modelType"))
  fittedData_curves[[drug]] = data_melted
  
  ## Plot the graph showing all curves for the drug
  color_model <- c("black", "deepskyblue1","chartreuse3")
  names(color_model) <- c("sigmoid", "linear", "exponential")
  png(paste0("./results/AUCs_production/fittedCurves/", drug, "_modelCurves.png"))
  print(
    ggplot( data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="lightgrey")) +
      labs(title = drug, x = "log(drug)", y = "% viability", colour = "Model", fill="SD") +
      scale_color_manual(values=color_model) +
      geom_line(data = data_melted, aes(x = as.numeric(as.character(variable)), y = value, group = sampNames, colour = modelType)) +
      scale_fill_gradient2(limits=c(0, maxSD), low="blue", mid="white", high="red", midpoint = maxSD/2) +
      geom_point(data = as.data.frame(stanDevs), mapping = aes(x=x, y=y, fill = sd, shape = 21, stroke=0)) +
      scale_shape_identity()# + geom_vline(3.401197, color='red')
    
    
  )
  dev.off()
  
}


## Plot the lines for samples that are excluded from analysis
excludedData_curves = list()
for(drug in names(fittedData)){
  print(drug)
  excludedSamples = unlist(fittedData[[drug]]["excludedSamples"]) # get list of samples that were not fitted
  names(excludedSamples) = gsub("excludedSamples.", "", names(excludedSamples))
  
  if(length(excludedSamples) > 0){
    
    drug.rawdata = viabs_perDrug[[drug]]
    excludedData = drug.rawdata[(gsub("_\\w+","",rownames(drug.rawdata)) %in% names(excludedSamples)), ]
    excludedSDs = (rawData_standardDeviations[[drug]])[names(excludedSamples),]
    
    
    ## Compute the mean for a given sample
    meansData = matrix(data = NA, nrow = length(excludedSamples), ncol = length(concentrations_num))
    rownames(meansData) = names(excludedSamples)
    for(samp in names(excludedSamples)){
      data.sample = drug.rawdata[gsub("_\\w+","",rownames(drug.rawdata)) == samp,]
      meansData[samp, ] = apply(data.sample, 2, mean, na.rm=T)
    }
    
    if(length(excludedSamples) > 1){
      
      colnames(excludedSDs) = log(concentrations_num)
      
      ## Extract mean and sd into one dataset
      drug.data = as.data.frame(matrix(data = NA, nrow = sum(!is.na(excludedSDs)), ncol = 5))
      colnames(drug.data) = c("sample", "x", "mean", "sd", "reason")
      idx = 1
      for(s in 1:nrow(excludedSDs)){
        for(con in 1:ncol(excludedSDs)){
          if(!is.na(excludedSDs[s, con])){
            drug.data[idx, 1] = names(excludedSamples)[s]
            drug.data[idx, 2] = log(concentrations_num[con])
            drug.data[idx, 3] = meansData[s, con]
            drug.data[idx, 4] = excludedSDs[s, con]
            drug.data[idx, 5] = ifelse((excludedSamples)[s]=="Unfittable", "Unfittable",
                                       ifelse((excludedSamples)[s]=="Inaccurate model", "Inaccurate", "Increasing"))
            idx = idx + 1
          }
        }
      }
      
    }else if(length(excludedSamples) == 1){
      
      names(excludedSDs) = log(concentrations_num)
      
      ## Extract mean and sd into one dataset
      drug.data = as.data.frame(matrix(data = NA, nrow = sum(!is.na(excludedSDs)), ncol = 5))
      names(drug.data) = c("sample", "x", "mean", "sd", "reason")
      idx = 1
      for(con in 1:length(excludedSDs)){
        if(!is.na(excludedSDs[con])){
          drug.data[idx, 1] = names(excludedSamples)
          drug.data[idx, 2] = log(concentrations_num[con])
          drug.data[idx, 3] = meansData[con]
          drug.data[idx, 4] = excludedSDs[con]
          drug.data[idx, 5] = ifelse(excludedSamples=="Unfittable", "Unfittable",
                                     ifelse(excludedSamples=="Inaccurate model", "Inaccurate", "Increasing"))
          idx = idx + 1
          
        }
      }
      
    }
    excludedData_curves[[drug]] = drug.data
    
    ## Plot the graph showing all curves for the excluded samples for this drug
    color_reason <- c("black", "orange","forestgreen")
    names(color_reason) <- c("Increasing","Unfittable", "Inaccurate")
    png(paste0("./results/AUCs_production/excludedSamples/", drug, "_unfittedCurves.png"))
    print(
      ggplot( data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="lightgrey")) +
        labs(title = drug, x = "log(drug)", y = "% viability", colour = "Model:", fill = "SD") +
        scale_color_manual(values=color_reason) +
        geom_line(data = drug.data, aes(x = x, y = mean, group = sample, colour=reason)) +
        scale_fill_gradient2(limits=c(0, maxSD), low="blue", mid="white", high="red", midpoint = maxSD/2) +
        geom_point(data = as.data.frame(drug.data), mapping = aes(x=x, y=mean, fill = sd, shape = 21, stroke=0)) +
        scale_shape_identity()
    )
    dev.off()

  }# end of if(length(excludedSamples) > 0)
  
}


# ===== III.6 - Compute AUCs  =====

AUCs = matrix(NA, nrow = length(samples_list), ncol = length(viabs_perDrug_filtered))
dimnames(AUCs) = list(samples_list, names(viabs_perDrug_filtered))
cutoff_drugs = c("Afatinib",    "Aminolevulinic acid hydrochloride", "Amiodarone", "Azacitidine", "Bisutinib", "Cabozantinib",
                 "Celecoxib",   "Crizotininb", "Enzalutamide",   "Exemestone",  "Imatinib",    "Lapatinib",  "Lomustine", "Mitotane",
                 "Nilotinib",   "Pipobroman",  "Raloxifene",     "Regorafenib", "Sirolimus",   "Sorafenib",  "Sunitinib", "Tamoxifen citrate",
                 "Thioguanine", "Tretinoin",   "Uracil mustard", "Vandetanib",  "Vemurafenib", "Vismodegib", "Vorinostat")
cutoff_values = c(-5,   -5, -2.5, -5, -5,   -5, 
                  -2.5, -5, -2.5, -5, -2.5, -5, -2.5, -2.5, 
                  -5,   -5, -5,   -5, -5,   -2.5, -5, -5, 
                  -5, -2.5, -5,   -5, -5,   -2.5, -5)

for(d in 1:ncol(AUCs)){ # for each drug
  print(colnames(AUCs)[d])
  
  data = as.data.frame(fittedData[[d]]["fittedValues"])
  metaData = as.data.frame(fittedData[[d]]["metaData"])
  
  if(colnames(AUCs)[d] %in% cutoff_drugs){  
    computeRange = which(log(newRange) > cutoff_values[which(cutoff_drugs == colnames(AUCs)[d])])
  }else{ computeRange = 1:length(newRange)  }
  
  for(s in 1:nrow(data)){ # for each sample
      
      auc = (sintegral(log(newRange[computeRange]), data[s, computeRange], 500))$int # AUC for the sample
      auc0 = (sintegral(log(newRange[computeRange]), rep(metaData[s,4], length(computeRange)), 500))$int # max AUC for the corresponding model
      AUCs[match(rownames(data)[s], rownames(AUCs)), d] = auc / auc0
    
  }
  
}

## Save
write.csv2(AUCs, file="./tmp/AUCs_drugExposure.csv")
save(fittedData, rawData_standardDeviations, file="./tmp/responseModels_drugExposure.RData")
save(rawData_curves, fittedData_curves, excludedData_curves, file = "./tmp/responseModels_curves.RData")

