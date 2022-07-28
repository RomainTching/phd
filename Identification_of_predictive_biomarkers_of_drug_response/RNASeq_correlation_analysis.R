# oooooooooooooooooooooooooooooooooooooooooooooooooo
# I - Set up environment ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo


library(readxl)
library(ComplexHeatmap)
library(drc)
library(Bolstad2)
library(ggplot2)
library(reshape)
library(DESeq2)
library(topGO)

setwd("/path/to/working_directory")


### ===== ++ 1 - Load and prepare the data  =====

# This section is data specific, so you'll need to replace with whatever steps are needed to get:
#    - your normalized RNA-Seq/omics dataset, with samples as lines and genes/transcripts/features as columns
#    - the response(s) you want to correlate with gene expression
#    - a "gene universe" matrix of two columns: one with the gene IDs as presented in the omics dataset, and the other with the gene symbols corresponding to these IDs. If the gene names are already used as column names of the omics dataset, just create a matrix with twice the same column.
# Sample names should match in both the omics and the responses datasets.
# I leave what I did as an example, but understand that IT DOESN'T NECESSARILY MATCH WHAT YOU NEED TO DO.
{
  ### ++++ a. RNA-Seq ----
  {
    ## Load data.
    # Response data and MGMT status.
    response_rna = as.data.frame(read_excel("./raw/080920_list with GLIOTRAIN samples for Romain_IN.xlsx", na="NA", sheet = 1, col_names = T))
    response_OmaCyta = as.data.frame(read_excel("./raw/doseresponse_omacetaxine_cytarabine.xlsx", na="NA", sheet = 1, col_names = T))
    response_rna = as.data.frame(cbind(response_rna, response_OmaCyta[, -1]))
    response_rna$GT.code <- gsub("01-WT", "02-CP", response_rna$GT.code)
    response_rna$GT.code <- gsub("-\\d+-\\w+-\\d+$", "", response_rna$GT.code)
    colnames(response_rna) <- c("GS.number", "GT.code", "TMZ_IC50", "TMZ_AUC", "TMZ_viability_100uM", "MGMT_status", "Oma_IC50", "Cyta_IC50", "Oma_AUC", "Cyta_AUC")
    rownames(response_rna) <- response_rna$GS.number
    
    # RNA-Seq data reduction and normalization. Leads to a processed file which is saved and can be loaded.
    {
      # # RNA-Seq data.
      # rnaseq_annots = read.table("../../../../GBMdatabase/GTdata/workingFolder/tmp/conversion_rnaSeqIDs.txt", h=T, sep='\t')
      # EMC_annots = read.csv("../../../../GBMAnalysis/EMCcollab/TMZ_Lomustine/Data/tmp/RNA_seq_annotation_Ensembl_BioMart_v89_protein_coding.csv", h=T)
      # rnaseq_raw = read.csv("../../../../GBMdatabase/GTdata/workingFolder/tmp/rnaSeq_raw_resequencedUpdated_200110.csv", header = T, row.names = 1)
      # rnaseq_samplesMap = read.csv("../../../../GBMdatabase/GTdata/workingFolder/tmp/samplesLabels_RNASeq.csv", h=T)
      # rnaseq_samplesBatch = read.csv("../../../../GBMdatabase/GTdata/workingFolder/tmp/samplesBatches_RNASeq.csv", h=T)
      # 
      # 
      # ## Prepare the raw count dataset.
      # colnames(rnaseq_raw) = rnaseq_samplesMap$Gliotrain_ID[match(colnames(rnaseq_raw), rnaseq_samplesMap$DILA_ID_RNA)]
      # rna_data_CP = rnaseq_raw[, grep("CP", colnames(rnaseq_raw))]
      # colnames(rna_data_CP) = gsub("-\\d+-\\w+-\\d+-\\w+$", "", colnames(rna_data_CP))
      # rna_data_CP = rna_data_CP[apply(rna_data_CP, 1, sd) != 0, ]
      # 
      # # Remove faulty cell lines.
      # faultyLines = c("GT-03-25", "GT-03-27", "GT-03-33", "GT-03-34", "GT-03-56", "GT-03-39", "GT-03-09")
      # rna_data_CP = rna_data_CP[ , !(colnames(rna_data_CP) %in% faultyLines)]
      # 
      # 
      # ## Reduce dataset.
      # rnaseq_annots$geneID[match(EMC_annots$Ensembl_ID, rnaseq_annots$id)] = EMC_annots$Gene_name # replace gene names in rnaseq_annots with those from EMC
      # 
      # # Make sure to only keep genes present in the data.
      # rnaseq_annots = rnaseq_annots[rnaseq_annots$id %in% rownames(rna_data_CP),]
      # EMC_annots = EMC_annots[EMC_annots$Ensembl_ID %in% rownames(rna_data_CP), ]
      # rnaseq_annots = rnaseq_annots[match(rownames(rna_data_CP), rnaseq_annots$id),] # put them in the same order as in the dataset
      # 
      # geneNames = unique(EMC_annots$Gene_name) # and determine how many unique coding-protein genes there are.
      # rna_data_reduced = as.data.frame(matrix(data = NA, nrow = length(geneNames), ncol = ncol(rna_data_CP),
      #                                         dimnames = list(geneNames, colnames(rna_data_CP)))) # create the empty dataframe that will receive the "reduced" data
      # 
      # for(g in 1:nrow(rna_data_reduced)){ # for each protein coding gene,
      #   if(g%%2000 == 0){  print(g)  }
      # 
      #   geneData = rna_data_CP[match(EMC_annots$Ensembl_ID[EMC_annots$Gene_name == rownames(rna_data_reduced)[g]], rownames(rna_data_CP)), ] # extract the subset of data corresponding to that gene
      # 
      #   # And insert that data in the reduced dataset, adding it per sample if there are more than one transcripts corresponding to that gene.
      #   if(nrow(geneData) == 1){
      #     rna_data_reduced[g, ] = geneData
      #   }else{
      #     combinedData = sapply(geneData, sum)
      #     cat(g, " combine\n")
      #     rna_data_reduced[g, ] = combinedData
      #   }
      # 
      # }
      # rna_data_reduced = rna_data_reduced[apply(rna_data_reduced, 1, sd) != 0, ]
      # 
      # 
      # ## Normalize with vst.
      # coldata_rna <- data.frame(GTid = colnames(rna_data_reduced))
      # ddsMat <- DESeqDataSetFromMatrix(countData = rna_data_reduced, colData = coldata_rna, design = ~1)
      # dim(ddsMat)
      # ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 10, ] # filter out transcripts that have less than 10 count total
      # dim(ddsMat)
      # vstMat <- DESeq2::vst(ddsMat, blind = FALSE)
      # rnaseq_normCP_combined <- t(assay(vstMat))
      # save(rnaseq_normCP_combined, EMC_annots, rnaseq_annots, file = "../../../GTdata_analysis/Data/tmp/rnaseq_normCP_duplicateCombined_protCoding.RData")
      # write.table(rnaseq_normCP_combined, file = "../../../GTdata_analysis/Data/tmp/rnaseq_normCP_duplicateCombined_protCoding.tsv", sep = '\t', col.names = NA)
    }
    load(file = "../../../GTdata_analysis/Data/tmp/rnaseq_normCP_duplicateCombined_protCoding.RData")
    
    ## Extract only data for TMZ analysis.
    rna_data <- rnaseq_normCP_combined #[match(response_rna$GT.code, rownames(rnaseq_normCP_combined)), ] # keep only non-faulty cell lines
    rownames(rna_data) <- response_rna$GS.number[match(rownames(rna_data), response_rna$GT.code)] # replace GT IDs with GS numbers
    rna_data_all <- rna_data[, sapply(1:ncol(rna_data), function(g) var(rna_data[,g])>0)] # remove genes that are constant across the cell lines
    rna_data_M <- rna_data_all[match(response_rna$GS.number[response_rna$MGMT_status=="Methylated"], rownames(rna_data_all)),]
    rna_data_M <- rna_data_M[, sapply(1:ncol(rna_data_M), function(g) var(rna_data_M[,g])>0)] # remove genes that are constant across these cell lines
    rna_data_UM <- rna_data_all[match(response_rna$GS.number[response_rna$MGMT_status=="Unmethylated"], rownames(rna_data_all)),]
    rna_data_UM <- rna_data_UM[, sapply(1:ncol(rna_data_UM), function(g) var(rna_data_UM[,g])>0)] # remove genes that are constant across these cell lines
    table(response_rna$MGMT_status)
    lapply(list(rna_data_all, rna_data_UM, rna_data_M), dim)
    colnames(rna_data_all)[!(colnames(rna_data_all) %in% colnames(rna_data_UM))]
    colnames(rna_data_all)[!(colnames(rna_data_all) %in% colnames(rna_data_M))]
    
    ## Mapping between transcripts ENSEMBL IDs and gene symbols.
    geneUniverse_rna <- data.frame(probe.name = EMC_annots$Gene_name, probe.symbol = EMC_annots$Gene_name) # to fit the runTopGO function
    geneUniverse_rna_all = geneUniverse_rna[match(colnames(rna_data_all), geneUniverse_rna$probe.name),]
    geneUniverse_rna_M = geneUniverse_rna[match(colnames(rna_data_M), geneUniverse_rna$probe.name),]
    geneUniverse_rna_UM = geneUniverse_rna[match(colnames(rna_data_UM), geneUniverse_rna$probe.name),]
    lapply(list(geneUniverse_rna_all, geneUniverse_rna_M, geneUniverse_rna_UM), dim)
    
  } # end of RNA-Seq-related data preparation

  
  ### ++++ b. DASL ----
  {
    ## Response data and MGMT status.
    response_dasl = as.data.frame(read.table("./raw/DASL_samples_response_data.txt", header = TRUE, sep = '\t'))
    rownames(response_dasl) = response_dasl$EMC_ID
    response_dasl <- response_dasl[, -c(1,2)]
    table(response_dasl$MGMT_status_culture)
    for(i in 1:nrow(response_dasl)){  
      if(!is.na(response_dasl$MGMT_status_culture[i])){
        if(response_dasl$MGMT_status_culture[i] == "Methylated (+UM)"){ 
          response_dasl$MGMT_status_culture[i] = "Methylated" 
        }  
      }
    }

    
    ## DASL data.
    #load("../../Ioannis_drugs/Data/tmp/transcFiltered_proteinCoding.RData")
    load("../../../GTdata_analysis/Data/tmp/analysesDatasets_EMC.DASL.RData")
    transcFiltered_proteinCoding = dasl_reduced
    
    
    ## Extract only data for TMZ analysis.
    sum(rownames(response_dasl) %in% rownames(transcFiltered_proteinCoding))
    rownames(response_dasl)[!(rownames(response_dasl) %in% rownames(transcFiltered_proteinCoding))]
    rownames(response_dasl)[!(rownames(response_dasl) %in% rownames(transcFiltered_proteinCoding))] <- c("GS104p", "GS279c", "GS280p", "GS323p", "GS343p", "GS79", "GS216c")
    rownames(response_dasl)[!(rownames(response_dasl) %in% rownames(transcFiltered_proteinCoding))]
    response_dasl_noNA = response_dasl[!is.na(response_dasl$MGMT_status_culture),]
    response_dasl_M = response_dasl_noNA[response_dasl_noNA$MGMT_status_culture == "Methylated", ]
    response_dasl_UM = response_dasl_noNA[response_dasl_noNA$MGMT_status_culture == "Unmethylated", ]
    
    dasl_all = as.data.frame(transcFiltered_proteinCoding[match(rownames(response_dasl), rownames(transcFiltered_proteinCoding)), ])
    dasl_M = dasl_all[match(rownames(response_dasl_M), rownames(dasl_all)), ]
    dasl_UM = dasl_all[match(rownames(response_dasl_UM), rownames(dasl_all)), ]
    
    
    ## Check out distribution of variances.
    var_all = sapply(dasl_all, var)
    var_M = sapply(dasl_M, var)
    var_UM = sapply(dasl_UM, var)
    vars = list(var_all, var_M, var_UM)
    par(mfrow = c(3, 1))
    sapply(vars, function(v) hist(v, breaks = 100))
    par(mfrow = c(1, 1))
    sapply(vars, min)
    
    
    ## Turn into matrices
    dasl_all = as.matrix(dasl_all)
    dasl_M = as.matrix(dasl_M)
    dasl_UM = as.matrix(dasl_UM)
    
    ## Mapping between transcripts ENSEMBL IDs and gene symbols.
    geneUniverse_dasl <- data.frame(probe.name = colnames(dasl_all), probe.symbol = colnames(dasl_all))
    
  } # end of DASL-related data preparation
  
}# end of data preprocessing.



### ===== ++ 2 - Define functions  =====

computeCorrelation <- function(dat, resp){
  ## Function testing each feature in the dataset for correlation with the response.
  # The rownames of the dataset and names of the response vector need to match!
  # @param <dat> the tested dataset, with samples as lines and features (genes, proteins...) as columns
  # @param <resp> the response vector tested for correlation
  # @return <corGenes> a dataframe containing the name of the gene, 
  #                                     correlation coefficient and method (Pearson/Spearman), 
  #                                     and p-value and q-value associated
  #                     sorted by increasing q-value.
  
  ## Create empty <corGenes>
  corGenes <- data.frame(gene = colnames(dat), method = rep(NA, ncol(dat)), cor.value = rep(NA, ncol(dat)),
                         p.value = rep(NA, ncol(dat)), q.value = rep(NA, ncol(dat)))
  
  
  ## Determine if Pearson or Spearman correlation should be used. 
  # Determine if <resp> vector is normally distributed.
  resp.test <- shapiro.test(resp)
  resp.isNormal <- (resp.test$p.value > 0.1)
  
  # Determine method accordingly.
  if(resp.isNormal){
    # for(g in 1:ncol(dat)){
    #   print(g)
    #   shapiro.pvalue <- (shapiro.test(dat[,g]))$p.value
    #   
    #   #print(shapiro.test(dat[,g]))
    #   corGenes$method[g] <- ifelse(shapiro.pvalue > 0.1, "pearson", "spearman")
    # }
    shapiro.pvalues <- sapply(1:ncol(dat), function(g) (shapiro.test(dat[,g]))$p.value)
    corGenes$method <- ifelse(shapiro.pvalues > 0.1, "pearson", "spearman")
  }else{  corGenes$method <- rep("spearman", nrow(corGenes))  } # if response is not normally distributed, only Spearman correlation can be computed
  
  ## Perform correlation test and FDR correction.
  for(g in 1:ncol(dat)){
    gCor <- cor.test(dat[, g], resp, method = corGenes$method[g])
    corGenes$cor.value[g] <- gCor$estimate
    corGenes$p.value[g] <- gCor$p.value
  }
  corGenes$q.value <- p.adjust(corGenes$p.value, "fdr")
  
  ## Sort result and end function.
  corGenes <- corGenes[order(corGenes$q.value, corGenes$p.value),]
  return(corGenes)
  
} # end of computeCorrelation()

runTopGO <- function(geneUniverse, investigatedGenes, outPath, title="TopGo_analysis", pvalues = NULL, ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20){
  {
    ## Performs an enrichment analysis using TopGO
    #
    # @params <geneUniverse> vector containing the complete list of genes symbols out of which the genes of interest from <investigatedGenes> were identified
    # @params <investigatedGenes> vector containing the list of genes of interest
    # @params <pvalues> vector containing the p-values for the investigated genes; relevant for methods relying on score, but can be left as NULL for methods relying on count
    # @params <ontology> the ontology to use for the analysis; can be one of "BP", "MF" or "CC"
    # @params <algorithm> name of the algorithm to explore the ontology; can be one of "classic", "elim", "weight", "weight01", "lea", "parentchild"; see package vignette
    # @params <algorithm> name of the statistical test used to compare GO terms; can be one of "fisher", "ks", "t", "globaltest", "sum"; see package vignette
    # @params <topNods> the number of top candidate GO terms to include in output visualization
  } # header of runTopGO()
  
  cat("=====  Performing enrichment analysis  =====\n")
  
  ## Transform the information from geneUniverse, investigatedGenes and pvalues into a gene list indicating the genes of interest.
  if(is.null(pvalues)){
    geneList = factor(as.integer(geneUniverse$probe.name %in% investigatedGenes))
    names(geneList) = geneUniverse$probe.symbol
    selectionFunction <- function(geneList){  return(geneList==1)  }
  }else{
    geneList = rep(1, nrow(geneUniverse))
    names(geneList) = geneUniverse$probe.symbol
    geneList[match(investigatedGenes, geneUniverse$probe.name)] = pvalues
    selectionFunction <- function(geneList){  return(geneList < 0.05)  }
  }
  
  ## Creating the topGOdata object.
  GOdata = new("topGOdata", ontology=ontology, allGenes=geneList, geneSel=selectionFunction,
               annot=annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  print(GOdata)
  
  
  ## Run enrichment tests
  resultenrichment <- runTest(GOdata, algorithm=algorithm, statistic=statistic)
  
  ## Visualize results
  summaryTableGo <- GenTable(GOdata, Pvalues=resultenrichment, orderBy="Pvalues", topNodes=topNods, numChar=100)
  write.csv(summaryTableGo, file = paste0(outPath, title, "_enrichmentAnalysisTopNodes.csv"))
  
  geneData(resultenrichment)
  
  showSigOfNodes(GOdata, score(resultenrichment), firstSigNodes = topNods, useInfo = "all")
  dev.copy2pdf(file = paste0(outPath, title, "_enrichmentAnalysisGraph.pdf"));
  dev.off()
  
  return(GOdata)
  
  
} # end of runTopGO()

runAnalysis <- function(dataset, resp, outPath, title, geneUniverse, filterStat = NULL, filterValue = 0.05){
  ## Function to run the several steps of the analysis: correlation, selection of significant genes,
  #                                                     enrichment analysis, heatmap and output.
  # @param <dataset> the expression dataset with samples as rows and genes as columns
  # @param <resp> the response vector (AUC, IC50...). names(resp) should match rownames(dataset)
  # @param <outPath> path to the output folder
  # @param <title> title of the run/analysis
  # @param <geneUniverse> a data frame containing a mapping between gene IDs in <dataset> and gene symbols
  # @param <filterStat> the statistic (column) used to determine genes of interest; can be either "cor.value", "correlation", "p.value", "p-value", "q.value" or "adjusted p-value"
  # @param <filterValue> the value used as a cut-off in the <filterStat> to determine genes of interest;
  #                                 - if <filterStat> is "cor.value" or "correlation", genes having ABOVE <filterValue> will be considered relevant
  #                                 - if <filterStat> is any other statistic (p- or q-value), genes that have the value BELOW it will be kept
  # @return <correlations_data> the dataframe containing correlation estimates and p-values
  
  
  ## Compute correlation and p-values between response and each gene.
  cat("===== Computing correlations =====\n")
  correlations_data <- computeCorrelation(dataset, resp)
  
  ## Define "significant genes".
  cat("===== Define significant genes =====\n")
  # Look at correlation values.
  cat("Here is the beginning of correlation analyses results, sorted by adjusted p-value significance:\n")
  print(head(correlations_data, 20))
  
  if(is.null(filterStat)){ # if the criterion to filter relevant genes is not predefined, ask the user.
    # Determine if interesting genes should be selected based on correlation coefficient, p-value or adjusted p-value.
    colSel <- TRUE
    while(colSel){
      columnSelection <- readline(prompt = "Which column should be used to define interesting genes? (cor.value = correlation, p.value = p-value, q.value = adjusted p-value): ")
      if(columnSelection %in% c("cor.value", "correlation", "p.value", "p-value", "q.value", "adjusted p-value")){  colSel <- FALSE  }
    }
    # Define cut-off for the value.
    valSel <- TRUE
    while(valSel){
      question <- ifelse(columnSelection %in% c("cor.value", "correlation"),
                         "Minimum absolute value for correlation to be of interest: ",
                         paste0(columnSelection, " must be lower than: "))
      thresholdSelection <- readline(prompt = question)
      thresholdSelection <- as.double(thresholdSelection)
      if(!is.na(thresholdSelection)){  valSel <- FALSE  }
    }
  }else{ # otherwise apply filtering parameters as provided
    columnSelection <- filterStat
    thresholdSelection <- filterValue
  }
  # Filter genes of interest.
  if(columnSelection %in% c("cor.value", "correlation")){
    correlations_signif <- correlations_data[ abs(correlations_data$cor.value) >= thresholdSelection, ]
  }else if(columnSelection %in% c("p.value", "p-value")){
    correlations_signif <- correlations_data[ correlations_data$p.value < thresholdSelection, ]
  }else{
    correlations_signif <- correlations_data[ correlations_data$q.value < thresholdSelection, ]
  }
  cat("*** ", nrow(correlations_signif), "genes selected ***\n")
  
  ## Produce heatmap.
  cat("===== Produce heatmap =====\n")
  dataset_signif <- dataset[, correlations_signif$gene]
  png(filename = paste0(outPath, title, "_heatmap.png"))
  heatmap_results <- heatmap(dataset_signif, labCol = FALSE)
  dev.off()
  
  ## Enrichment analysis.
  runTopGO(geneUniverse = geneUniverse, investigatedGenes = correlations_signif$gene,
           outPath = outPath, title = title)
  
  ## Save.
  cat("===== Wrap up =====\n")
  save(correlations_data, correlations_signif, dataset_signif, heatmap_results, 
       file = paste0(outPath, title, "_results.RData"))
  write.table(correlations_signif, file = paste0(outPath, title, "_significant_correlations.tsv"), sep = '\t', row.names = FALSE)
  
  cat("End of analysis ", title, "\n")
  return(correlations_data)
  
} # end of runAnalysis()



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# II - Run Analyses ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  results.dir <- "./results/2021.06.22_pval.05/"
  results.stat <- "p.value"
  results.val <- 0.05
  
  ### ===== ++ 1 - Run analysis for RNA-Seq data =====
  
  ### Correlation with AUC.
  resp <- response_rna$TMZ_AUC
  names(resp) <- rownames(response_rna)
  resp <- resp[match(rownames(rna_data_all), names(resp))]
  cor_TMZAUC_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
                                outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_TMZAUC_all",
                                filterStat = results.stat, filterValue = results.val)
  cor_TMZAUC_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
                               geneUniverse = geneUniverse_rna_UM, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_TMZAUC_UM",
                               filterStat = results.stat, filterValue = results.val)
  cor_TMZAUC_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
                              geneUniverse = geneUniverse_rna_M, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_TMZAUC_M",
                              filterStat = results.stat, filterValue = results.val)
  ## All samples.
  # cor_AUC_all <- computeCorrelation(rna_data_all, resp)
  # print(cor_AUC_all)
  # cor_AUC_allSignif <- cor_AUC_all[cor_AUC_all$p.value < 0.05, ]
  # write.table(cor_AUC_allSignif, file = "./results/2020.10.01/rnaseq_AUC_all_correlations")
  # geneUniverse <- rnaseq_annots
  # colnames(geneUniverse) <- c("probe.name", "probe.symbol") # to fit the runTopGO function
  # runTopGO(geneUniverse = geneUniverse, investigatedGenes = cor_AUC_allSignif$gene,
  #          outPath = "./results/2020.10.01/", title = "rnaseq_AUC_all")
  # data_AUC_allSignif <- rna_data_all[, cor_AUC_allSignif$gene]
  # png(filename = "./results/2020.10.01/rnaseq_AUC_all_heatmap.png")
  # heatmap(data_AUC_allSignif, labCol = FALSE)
  # dev.off()
  
  
  ### Correlation with IC50.
  # resp <- response_rna$TMZ_IC50
  # names(resp) <- rownames(response_rna)
  # cor_IC50_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
  #                             outPath = "./results/2021.01.25/rnaseq/", title = "rnaseq_IC50_all")
  # cor_IC50_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
  #                            geneUniverse = geneUniverse_rna_UM, outPath = "./results/2021.01.25/rnaseq/", title = "rnaseq_IC50_UM")
  # cor_IC50_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
  #                           geneUniverse = geneUniverse_rna_M, outPath = "./results/2021.01.25/rnaseq/", title = "rnaseq_IC50_M")
  
  
  ### Correlation with viability at 100uM.
  resp <- response_rna$TMZ_viability_100uM
  names(resp) <- rownames(response_rna)
  resp <- resp[match(rownames(rna_data_all), names(resp))]
  cor_100uM_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
                               outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_TMZ100uM_all",
                               filterStat = results.stat, filterValue = results.val)
  cor_100uM_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
                              geneUniverse = geneUniverse_rna_UM, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_TMZ100uM_UM",
                              filterStat = results.stat, filterValue = results.val)
  cor_100uM_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
                             geneUniverse = geneUniverse_rna_M, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_TMZ100uM_M",
                             filterStat = results.stat, filterValue = results.val)
  
  
  ### Correlation with Omacetaxine AUC.
  resp <- response_rna$Oma_AUC
  names(resp) <- rownames(response_rna)
  resp <- resp[match(rownames(rna_data_all), names(resp))]
  cor_OmaAUC_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
                                outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_OmaAUC_all",
                                filterStat = results.stat, filterValue = results.val)
  cor_OmaAUC_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
                               geneUniverse = geneUniverse_rna_UM, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_OmaAUC_UM",
                               filterStat = results.stat, filterValue = results.val)
  cor_OmaAUC_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
                              geneUniverse = geneUniverse_rna_M, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_OmaAUC_M",
                              filterStat = results.stat, filterValue = results.val)
  
  ### Correlation with Omacetaxine IC50.
  resp <- response_rna$Oma_IC50
  names(resp) <- rownames(response_rna)
  resp <- resp[match(rownames(rna_data_all), names(resp))]
  cor_OmaIC50_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
                                outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_OmaIC50_all",
                                filterStat = results.stat, filterValue = results.val)
  cor_OmaIC50_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
                               geneUniverse = geneUniverse_rna_UM, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_OmaIC50_UM",
                               filterStat = results.stat, filterValue = results.val)
  cor_OmaIC50_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
                              geneUniverse = geneUniverse_rna_M, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_OmaIC50_M",
                              filterStat = results.stat, filterValue = results.val)
  
  
  ### Correlation with Cytarabine AUC.
  resp <- response_rna$Cyta_AUC
  names(resp) <- rownames(response_rna)
  resp <- resp[match(rownames(rna_data_all), names(resp))]
  cor_CytaAUC_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
                                outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_CytaAUC_all",
                                filterStat = results.stat, filterValue = results.val)
  cor_CytaAUC_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
                               geneUniverse = geneUniverse_rna_UM, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_CytaAUC_UM",
                               filterStat = results.stat, filterValue = results.val)
  cor_CytaAUC_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
                              geneUniverse = geneUniverse_rna_M, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_CytaAUC_M",
                              filterStat = results.stat, filterValue = results.val)
  
  ### Correlation with Cytarabine IC50.
  resp <- response_rna$Cyta_IC50
  names(resp) <- rownames(response_rna)
  resp <- resp[match(rownames(rna_data_all), names(resp))]
  cor_CytaIC50_all <- runAnalysis(dataset = rna_data_all, resp = resp, geneUniverse = geneUniverse_rna_all,
                                 outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_CytaIC50_all",
                                 filterStat = results.stat, filterValue = results.val)
  cor_CytaIC50_UM <- runAnalysis(dataset = rna_data_UM, resp = resp[match(rownames(rna_data_UM), names(resp))], 
                                geneUniverse = geneUniverse_rna_UM, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_CytaIC50_UM",
                                filterStat = results.stat, filterValue = results.val)
  cor_CytaIC50_M <- runAnalysis(dataset = rna_data_M, resp = resp[match(rownames(rna_data_M), names(resp))],
                               geneUniverse = geneUniverse_rna_M, outPath = paste0(results.dir, "rnaseq/"), title = "rnaseq_CytaIC50_M",
                               filterStat = results.stat, filterValue = results.val)
  
  
  
  ### ===== ++ 2 - Run analysis for DASL data =====
  
  ### Correlation with AUC.
  resp <- response_dasl$AUC
  names(resp) <- rownames(response_dasl)
  resp <- resp[match(rownames(dasl_all), names(resp))]
  resp_M <- response_dasl_M$AUC
  names(resp_M) <- rownames(response_dasl_M)
  resp_UM <- response_dasl_UM$AUC
  names(resp_UM) <- rownames(response_dasl_UM)
  
  cor_AUC_all <- runAnalysis(dataset = dasl_all, resp = resp, geneUniverse = geneUniverse_dasl,
                             outPath = paste0(results.dir, "dasl/"), title = "dasl_AUC_all",
                             filterStat = results.stat, filterValue = results.val)
  cor_AUC_UM <- runAnalysis(dataset = dasl_UM, resp = resp_UM, geneUniverse = geneUniverse_dasl, 
                            outPath = paste0(results.dir, "dasl/"), title = "dasl_AUC_UM",
                            filterStat = results.stat, filterValue = results.val)
  cor_AUC_M <- runAnalysis(dataset = dasl_M, resp = resp_M, geneUniverse = geneUniverse_dasl, 
                           outPath = paste0(results.dir, "dasl/"), title = "dasl_AUC_M",
                           filterStat = results.stat, filterValue = results.val)
  
  
  ### Correlation with IC50.
  # resp <- response_dasl$IC50_uM
  # names(resp) <- rownames(response_dasl)
  # resp_M <- response_dasl_M$IC50_uM
  # names(resp_M) <- rownames(response_dasl_M)
  # resp_UM <- response_dasl_UM$IC50_uM
  # names(resp_UM) <- rownames(response_dasl_UM)
  # 
  # cor_IC50_all <- runAnalysis(dataset = dasl_all, resp = resp, geneUniverse = geneUniverse_dasl,
  #                                      outPath = "./results/2021.01.25/dasl/", title = "dals_IC50_all")
  # cor_IC50_UM <- runAnalysis(dataset = dasl_UM, resp = resp_UM, geneUniverse = geneUniverse_dasl, 
  #                                     outPath = "./results/2021.01.25/dasl/", title = "dals_IC50_UM")
  # cor_IC50_M <- runAnalysis(dataset = dasl_M, resp = resp_M, geneUniverse = geneUniverse_dasl, 
  #                                    outPath = "./results/2021.01.25/dasl/", title = "dals_IC50_M")
  
  
  ### Correlation with viability at 80uM.
  resp <- response_dasl$percent_cell_viability_at_80uM
  names(resp) <- rownames(response_dasl)
  resp <- resp[match(rownames(dasl_all), names(resp))]
  resp_M <- response_dasl_M$percent_cell_viability_at_80uM
  names(resp_M) <- rownames(response_dasl_M)
  resp_UM <- response_dasl_UM$percent_cell_viability_at_80uM
  names(resp_UM) <- rownames(response_dasl_UM)
  
  cor_80uM_all <- runAnalysis(dataset = dasl_all, resp = resp, geneUniverse = geneUniverse_dasl,
                              outPath = paste0(results.dir, "dasl/"), title = "dasl_80uM_all",
                              filterStat = results.stat, filterValue = results.val)
  cor_80uM_UM <- runAnalysis(dataset = dasl_UM, resp = resp_UM, geneUniverse = geneUniverse_dasl, 
                             outPath = paste0(results.dir, "dasl/"), title = "dasl_80uM_UM",
                             filterStat = results.stat, filterValue = results.val)
  cor_80uM_M <- runAnalysis(dataset = dasl_M, resp = resp_M, geneUniverse = geneUniverse_dasl, 
                            outPath = paste0(results.dir, "dasl/"), title = "dasl_80uM_M",
                            filterStat = results.stat, filterValue = results.val)
  
  
  
} # end of II - Run Analyses



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# III - Compile results ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  ## Define all names that will be used to read files and extract relevant data.
  run_titles_rna <- c("rnaseq_TMZ100uM", "rnaseq_TMZAUC", "rnaseq_CytaAUC", 
                      "rnaseq_CytaIC50", "rnaseq_OmaAUC", "rnaseq_OmaIC50")
  run_titles_dasl <- c("dasl_80uM", "dasl_AUC")
  run_cohorts <- c("all", "M", "UM")
  
  
  ## Compile significant genes together.
  significant_genes <- list()
  for(run in run_titles_rna){
    titles = paste(run, run_cohorts, sep = '_')
    for(titleRun in titles){
      correlations = as.data.frame(read.table(file = paste0(results.dir, "rnaseq/", titleRun, "_significant_correlations.tsv"), header = TRUE, sep = '\t'))
      significant_genes[[titleRun]] = correlations[, "gene"]
    }
  }
  for(run in run_titles_dasl){
    titles = paste(run, run_cohorts, sep = '_')
    for(titleRun in titles){
      correlations = as.data.frame(read.table(file = paste0(results.dir, "dasl/", titleRun, "_significant_correlations.tsv"), header = TRUE, sep = '\t'))
      significant_genes[[titleRun]] = correlations[, "gene"]
    }
  }
  maxGeneNumber = max(sapply(significant_genes, length))
  results_significantGenes = matrix(data = NA, ncol = length(significant_genes), nrow = maxGeneNumber)
  for(r in 1:length(significant_genes)){
    results_significantGenes[ 1:length(significant_genes[[r]]), r] = significant_genes[[r]]
  }
  colnames(results_significantGenes) = names(significant_genes)
  results_significantGenes[1:5, 1:5]
  
  
  ## Compile enrichment analysis results form significant genes together.
  enrichment_genes <- list()
  for(run in run_titles_rna){
    titles = paste(run, run_cohorts, sep = '_')
    for(titleRun in titles){
      GOresults = as.data.frame(read.csv(file = paste0(results.dir, "rnaseq/", titleRun, "_enrichmentAnalysisTopNodes.csv"), header = TRUE))
      enrichment_genes[[titleRun]] = GOresults[, "Term"]
    }
  }
  for(run in run_titles_dasl){
    titles = paste(run, run_cohorts, sep = '_')
    for(titleRun in titles){
      GOresults = as.data.frame(read.csv(file = paste0(results.dir, "dasl/", titleRun, "_enrichmentAnalysisTopNodes.csv"), header = TRUE))
      enrichment_genes[[titleRun]] = GOresults[, "Term"]
    }
  }
  maxGeneNumber = max(sapply(enrichment_genes, length))
  results_enrichmentGenes = matrix(data = NA, ncol = length(enrichment_genes), nrow = maxGeneNumber)
  for(r in 1:length(enrichment_genes)){
    results_enrichmentGenes[ 1:length(enrichment_genes[[r]]), r] = enrichment_genes[[r]]
  }
  colnames(results_enrichmentGenes) = names(enrichment_genes)
  results_enrichmentGenes[1:5, 1:5]
  
  
  ## Define "signatures".
  # As genes that are found significant in runs on the same dataset, cohort and drug (but response data is different).
  runCouples = c("rnaseq_TMZ", "dasl", "rnaseq_Cyta", "rnaseq_Oma")
  signature_genes = list()
  for(run in runCouples){
    drugSubset = results_significantGenes[, grep(run, colnames(results_significantGenes))]
    for(cohort in run_cohorts){
      cohortSubset = drugSubset[, grep(paste0('_',cohort), colnames(drugSubset))]
      signature = intersect(na.omit(cohortSubset[, 1]), na.omit(cohortSubset[, 2]))
      signature_genes[[paste(run, cohort, sep = '_')]] = signature
    }
  }
  maxGeneNumber = max(sapply(signature_genes, length))
  results_signatureGenes = matrix(data = NA, ncol = length(signature_genes), nrow = maxGeneNumber)
  for(r in 1:length(signature_genes)){
    if(length(signature_genes[[r]]) > 0 ){
      results_signatureGenes[ 1:length(signature_genes[[r]]), r] = signature_genes[[r]]
    }
  }
  colnames(results_signatureGenes) = names(signature_genes)
  results_signatureGenes[1:5, 1:5]
  
  
  ## Run enrichment analyses on signature genes.
  # Run enrichment analyses.
  gene_universes = list(geneUniverse_rna_all, geneUniverse_rna_M, geneUniverse_rna_UM, # rnaseq_TMZ
                        geneUniverse_dasl, geneUniverse_dasl, geneUniverse_dasl, # dasl
                        geneUniverse_rna_all, geneUniverse_rna_M, geneUniverse_rna_UM, # rnaseq_Oma
                        geneUniverse_rna_all, geneUniverse_rna_M, geneUniverse_rna_UM) # rnaseq_Cyta
  for(i in 1:ncol(results_signatureGenes)){
    if(length(na.omit(results_signatureGenes[, i])) > 0){
      runTopGO(geneUniverse = gene_universes[[i]], investigatedGenes = na.omit(results_signatureGenes[, i]),
               outPath = paste0(results.dir,"signatures_enrichment/"), title = colnames(results_signatureGenes)[i])
    }
  }
  
  # Compile at the same place.
  enrichment_signatures <- list()
  # for(run in runCouples){
  #   titles = paste(run, run_cohorts, sep = '_')
  #   for(titleRun in titles){
  #     if()
  for(i in 1:ncol(results_signatureGenes)){
    if(length(na.omit(results_signatureGenes[, i])) > 1){
      GOresults = as.data.frame(read.csv(file = paste0(results.dir, "signatures_enrichment/", colnames(results_signatureGenes)[i], "_enrichmentAnalysisTopNodes.csv"), header = TRUE))
      enrichment_signatures[[colnames(results_signatureGenes)[i]]] = GOresults[, "Term"]
    }
  }
  maxGeneNumber = max(sapply(enrichment_signatures, length))
  results_enrichmentSignatures = matrix(data = NA, ncol = length(enrichment_signatures), nrow = maxGeneNumber)
  for(r in 1:length(enrichment_signatures)){
    results_enrichmentSignatures[ 1:length(enrichment_signatures[[r]]), r] = enrichment_signatures[[r]]
  }
  colnames(results_enrichmentSignatures) = names(enrichment_signatures)
  results_enrichmentSignatures[1:5, 1:5]
  
  
  ## Save tables.
  write.table(results_significantGenes, file = paste0(results.dir, "results_significantGenes.tsv"), sep = '\t', row.names = FALSE, na = "")
  write.table(results_enrichmentGenes, file = paste0(results.dir, "results_enrichmentGenes.tsv"), sep = '\t', row.names = FALSE, na = "")
  write.table(results_signatureGenes, file = paste0(results.dir, "results_signatureGenes.tsv"), sep = '\t', row.names = FALSE, na = "")
  write.table(results_enrichmentSignatures, file = paste0(results.dir, "results_enrichmentSignatures.tsv"), sep = '\t', row.names = FALSE, na = "")
  
} # end of III - compare results



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# IV - Validate signatures ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

{
  results_significantGenes = read.table(file = paste0(results.dir, "results_significantGenes.tsv"), sep = '\t', header = T)
  results_enrichmentGenes = read.table(file = paste0(results.dir, "results_enrichmentGenes.tsv"), sep = '\t', header = T)
  results_signatureGenes = read.table(file = paste0(results.dir, "results_signatureGenes.tsv"), sep = '\t', header = T)
  results_enrichmentSignatures = read.table(file = paste0(results.dir, "results_enrichmentSignatures.tsv"), sep = '\t', header = T)
  signatures = as.data.frame(results_signatureGenes)
  names(signatures)
  
  survival_rna = read.table("./tmp/TMZstudy_survival.csv", sep = ',', h=T, row.names = 1)
  survivals_rna = survival_rna[, 2]
  names(survivals_rna) = rownames(survival_rna)
  
  load("../../Ioannis_drugs/Data/tmp/samplesLists.RData")
  survivals_dasl = diagnosis$OS[match(rownames(dasl_all), rownames(diagnosis))]
  names(survivals_dasl) = rownames(diagnosis)[match(rownames(dasl_all), rownames(diagnosis))]

  
  ### ++ 1 - Validate by computing combination score ----
  {
    validateSignature <- function(signature, dataset, survival, bootstrap = 500){
      ## Function to test identified gene signature in the dataset.
      # For each sample, abs(expression value) for all genes in the signature are multiplied. 
      # Distribution of these signatures is compared to a bootstrap of fake "signatures" calculated 
      # from randomly selected genes in the dataset.
      #
      # @param <signature> vector of genes composing the signature
      # @param <dataset> data used for validation
      # @param <bootstrap> number of randomly generated "signatures" used for validation
      # @param <title> the title used for exporting files
      # @return a list containing 
      #                - the genes used for both the tested signature and the randomly generated ones
      #                - the corresponding score 
      #                - 
      
      
      ## Make sure only genes existing in <dataset> are kept in <signature>
      if(any(!(signature %in% colnames(dataset)))){
        cat("!! ", signature[which(!(signature %in% colnames(dataset)))], 
            " were removed from tested signature since they are not available in the validation dataset !!\n")
        signature = signature[which((signature %in% colnames(dataset)))]
      }
      cat("Signature length = ", length(signature), '\n')
      
      
      ## Generate random signatures of the same length as tested signature.
      signaturesList = sapply(1:bootstrap, function(s) sample(colnames(dataset), length(signature)))
      signaturesList = cbind(signature, signaturesList)
      colnames(signaturesList) = c("Signature", paste0("Random_", 1:bootstrap))
      
      
      ## Compute scores.
      signaturesScores = matrix(data = NA, nrow = nrow(dataset), ncol = ncol(signaturesList),
                                dimnames = list(rownames(dataset), colnames(signaturesList)))
      for(s in 1:ncol(signaturesList)){
        
        dat = dataset[, match(signaturesList[,s], colnames(dataset))] # get data subset
        signaturesScores[, s] = apply(dat, # calculate all scores for this signature by
                                      1, # for each sample
                                      function(samp) mean(abs(samp))) # calculate the average of the absolute value for all genes in subset
        
      }
      
      
      ## Determine correlation with <survival>.
      survivalCorrelations = matrix(data = NA, nrow = ncol(signaturesList), ncol = 2,
                                    dimnames = list(colnames(signaturesList), c("cor.estimate","p.value")))
      for(s in 1:nrow(survivalCorrelations)){
        test = cor.test(signaturesScores[,s], survival[match(rownames(signaturesScores), names(survival))], method = "spearman")
        survivalCorrelations[s,1] = test$estimate
        survivalCorrelations[s,2] = test$p.value
      }
      
      return(list(signaturesList, signaturesScores, survivalCorrelations))
      
    }
    
    bootstrapNB = 500
    
    ## ++++ a. RNA-Seq signatures for TMZ ----
    {
      # With all identified genes.
      validation_RNASeq_TMZ_all_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_all),
                                                           survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                           dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_all_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_TMZ_all_RNASeq[[3]][1:2,]
      
      
      validation_RNASeq_TMZ_UM_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_UM),
                                                          survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
                                                          dataset = rna_data_UM, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_UM_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_TMZ_UM_RNASeq[[3]][1:2,]
      
      validation_RNASeq_TMZ_M_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_M),
                                                         survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                         dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_M_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_TMZ_M_RNASeq[[3]][1:2,]
      
      # With first few genes.
      validation_RNASeq_TMZ_all_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_all[1:20]),
                                                           survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                           dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_all_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_TMZ_all_RNASeq[[3]][1:2,]
      
      validation_RNASeq_TMZ_UM_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_UM[1:5]),
                                                          survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
                                                          dataset = rna_data_UM, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_UM_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_TMZ_UM_RNASeq[[3]][1:2,]
      
      validation_RNASeq_TMZ_M_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_M[1:60]),
                                                         survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                         dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_M_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_TMZ_M_RNASeq[[3]][1:2,]
    }
    
    ## ++++ b. RNA-Seq signatures for Omacetaxine ----
    {
      # With all identified genes.
      validation_RNASeq_Oma_all_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Oma_all),
                                                           survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                           dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Oma_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Oma_all_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Oma_all_RNASeq[[3]][1:2,]
      
      # validation_RNASeq_Oma_UM_RNASeq = validateSignature(signature = na.omit(signatures$RNASeq_Oma_UM),
      #                                                     survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
      #                                                     dataset = rna_data_UM, bootstrap = bootstrapNB)
      # hist(validation_RNASeq_Oma_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      # abline(v=validation_RNASeq_Oma_UM_RNASeq[[2]][,1], col = "red")
      
      validation_RNASeq_Oma_M_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Oma_M),
                                                         survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                         dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Oma_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Oma_M_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Oma_M_RNASeq[[3]][1:2,]
      
      # With first few genes.
      validation_RNASeq_Oma_all_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Oma_all[1:12]),
                                                           survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                           dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Oma_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Oma_all_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Oma_all_RNASeq[[3]][1:2,]
      
      # validation_RNASeq_Oma_UM_RNASeq = validateSignature(signature = na.omit(signatures$RNASeq_Oma_UM[1:10]),
      #                                                     survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
      #                                                     dataset = rna_data_UM, bootstrap = 100)
      # hist(validation_RNASeq_Oma_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      # abline(v=validation_RNASeq_Oma_UM_RNASeq[[2]][,1], col = "red")
      
      validation_RNASeq_Oma_M_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Oma_M[1:10]),
                                                         survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                         dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Oma_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Oma_M_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Oma_M_RNASeq[[3]][1:2,]
    }
    
    ## ++++ c. RNA-Seq signatures for Cytarabine ----
    {
      # With all identified genes.
      validation_RNASeq_Cyta_all_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_all),
                                                            survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                            dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_all_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Cyta_all_RNASeq[[3]][1:2,]
      
      validation_RNASeq_Cyta_UM_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_UM),
                                                           survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
                                                           dataset = rna_data_UM, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_UM_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Cyta_UM_RNASeq[[3]][1:2,]
      
      validation_RNASeq_Cyta_M_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_M),
                                                          survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                          dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_M_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Cyta_M_RNASeq[[3]][1:2,]
      
      # With first few genes.
      validation_RNASeq_Cyta_all_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_all[1:15]),
                                                            survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                            dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_all_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Cyta_all_RNASeq[[3]][1:2,]
      
      validation_RNASeq_Cyta_UM_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_UM[1:30]),
                                                           survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
                                                           dataset = rna_data_UM, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_UM_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Cyta_UM_RNASeq[[3]][1:2,]
      
      validation_RNASeq_Cyta_M_RNASeq = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_M[1:12]),
                                                          survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                          dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_M_RNASeq[[2]][,1], col = "red")
      validation_RNASeq_Cyta_M_RNASeq[[3]][1:2,]
    }
    
    ## ++++ d. DASL signatures ----
    {
      # With all identified genes.
      validation_DASL_TMZ_all_DASL = validateSignature(signature = na.omit(signatures$dals_all),
                                                       survival = survivals_dasl[match(rownames(dasl_all), names(survivals_dasl))],
                                                       dataset = dasl_all, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_all_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_all_DASL[[2]][,1], col = "red")
      validation_DASL_TMZ_all_DASL[[3]][1:2,]
      
      validation_DASL_TMZ_UM_DASL = validateSignature(signature = na.omit(signatures$dals_UM),
                                                      survival = survivals_dasl[match(rownames(dasl_UM), names(survivals_dasl))],
                                                      dataset = dasl_UM, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_UM_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_UM_DASL[[2]][,1], col = "red")
      validation_DASL_TMZ_UM_DASL[[3]][1:2,]
      
      validation_DASL_TMZ_M_DASL = validateSignature(signature = na.omit(signatures$dals_M),
                                                     survival = survivals_dasl[match(rownames(dasl_M), names(survivals_dasl))],
                                                     dataset = dasl_M, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_M_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_M_DASL[[2]][,1], col = "red")
      validation_DASL_TMZ_M_DASL[[3]][1:2,]
      
      
      # With first few genes.
      validation_DASL_TMZ_all_DASL = validateSignature(signature = na.omit(signatures$dals_all[1:27]), # can't find good subset
                                                       survival = survivals_dasl[match(rownames(dasl_all), names(survivals_dasl))],
                                                       dataset = dasl_all, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_all_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_all_DASL[[2]][,1], col = "red")
      validation_DASL_TMZ_all_DASL[[3]][1:2,]
      
      validation_DASL_TMZ_UM_DASL = validateSignature(signature = na.omit(signatures$dals_UM[1:36]), # all genes in signature
                                                      survival = survivals_dasl[match(rownames(dasl_UM), names(survivals_dasl))],
                                                      dataset = dasl_UM, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_UM_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_UM_DASL[[2]][,1], col = "red")
      validation_DASL_TMZ_UM_DASL[[3]][1:2,]
      
      validation_DAS_TMZL_M_DASL = validateSignature(signature = na.omit(signatures$dals_M[1:25]),# problem with bootstrap?
                                                     survival = survivals_dasl[match(rownames(dasl_M), names(survivals_dasl))],
                                                     dataset = dasl_M, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_M_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_M_DASL[[2]][,1], col = "red")
      validation_DASL_TMZ_M_DASL[[3]][1:2,]
    }
    
    ## ++++ e. Cross-datasets ----
    {
      # RNAseq TMZ on DASL.
      validation_RNASeq_TMZ_all_DASL = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_all[1:1000]), # all genes in signature, +
                                                           survival = survivals_dasl[match(rownames(dasl_all), names(survivals_dasl))],
                                                           dataset = dasl_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_all_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_all_DASL[[2]][,1], col = "red")
      validation_RNASeq_TMZ_all_DASL[[3]][1:2,]
      
      validation_RNASeq_TMZ_UM_DASL = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_UM[1:124]), # good, ++
                                                        survival = survivals_dasl[match(rownames(dasl_UM), names(survivals_dasl))],
                                                        dataset = dasl_UM, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_UM_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_UM_DASL[[2]][,1], col = "red")
      validation_RNASeq_TMZ_UM_DASL[[3]][1:2,]
      
      validation_RNASeq_TMZ_M_DASL = validateSignature(signature = na.omit(signatures$rnaseq_TMZ_M[1:38]), # can't find good one
                                                         survival = survivals_dasl[match(rownames(dasl_M), names(survivals_dasl))],
                                                         dataset = dasl_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_TMZ_M_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_TMZ_M_DASL[[2]][,1], col = "red")
      validation_RNASeq_TMZ_M_DASL[[3]][1:2,]
      
      # RNAseq Cytarabine on DASL.
      validation_RNASeq_Cyta_all_DASL = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_all[1:30]), # can't find good one
                                                         survival = survivals_dasl[match(rownames(dasl_all), names(survivals_dasl))],
                                                         dataset = dasl_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_all_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_all_DASL[[2]][,1], col = "red")
      validation_RNASeq_Cyta_all_DASL[[3]][1:2,]
      
      validation_RNASeq_Cyta_UM_DASL = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_UM[1:3]), # fine, +
                                                        survival = survivals_dasl[match(rownames(dasl_UM), names(survivals_dasl))],
                                                        dataset = dasl_UM, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_UM_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_UM_DASL[[2]][,1], col = "red")
      validation_RNASeq_Cyta_UM_DASL[[3]][1:2,]
      
      validation_RNASeq_Cyta_M_DASL = validateSignature(signature = na.omit(signatures$rnaseq_Cyta_M[1:42]), # fine, +
                                                       survival = survivals_dasl[match(rownames(dasl_M), names(survivals_dasl))],
                                                       dataset = dasl_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Cyta_M_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Cyta_M_DASL[[2]][,1], col = "red")
      validation_RNASeq_Cyta_M_DASL[[3]][1:2,]
      
      # RNAseq Omacetaxine on DASL.
      validation_RNASeq_Oma_all_DASL = validateSignature(signature = na.omit(signatures$rnaseq_Oma_all[1:18]), # fine, -
                                                          survival = survivals_dasl[match(rownames(dasl_all), names(survivals_dasl))],
                                                          dataset = dasl_all, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Oma_all_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Oma_all_DASL[[2]][,1], col = "red")
      validation_RNASeq_Oma_all_DASL[[3]][1:2,]
      
      validation_RNASeq_Oma_M_DASL = validateSignature(signature = na.omit(signatures$rnaseq_Oma_M[1:10]), # can't find a good one
                                                        survival = survivals_dasl[match(rownames(dasl_M), names(survivals_dasl))],
                                                        dataset = dasl_M, bootstrap = bootstrapNB)
      hist(validation_RNASeq_Oma_M_DASL[[2]], breaks = 100, main = "Scores for RNA-Seq_M gene signature", xlab = "Score")
      abline(v=validation_RNASeq_Oma_M_DASL[[2]][,1], col = "red")
      validation_RNASeq_Oma_M_DASL[[3]][1:2,]
      
      # DASL TMZ on RNASeq.
      validation_DASL_TMZ_all_RNASeq = validateSignature(signature = na.omit(signatures$dals_all[1:150]), # fine, -- spread
                                                         survival = survivals_rna[match(rownames(rna_data_all), names(survivals_rna))],
                                                         dataset = rna_data_all, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_all_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_all_RNASeq[[2]][,1], col = "red")
      validation_DASL_TMZ_all_RNASeq[[3]][1:2,]
      
      validation_DASL_TMZ_UM_RNASeq = validateSignature(signature = na.omit(signatures$dals_UM), # all genes in signature
                                                         survival = survivals_rna[match(rownames(rna_data_UM), names(survivals_rna))],
                                                         dataset = rna_data_UM, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_UM_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_UM_RNASeq[[2]][,1], col = "red")
      validation_DASL_TMZ_UM_RNASeq[[3]][1:2,]
      
      validation_DASL_TMZ_M_RNASeq = validateSignature(signature = na.omit(signatures$dals_M[1:4]), # can't find a good one.
                                                        survival = survivals_rna[match(rownames(rna_data_M), names(survivals_rna))],
                                                        dataset = rna_data_M, bootstrap = bootstrapNB)
      hist(validation_DASL_TMZ_M_RNASeq[[2]], breaks = 100, main = "Scores for RNA-Seq_all gene signature", xlab = "Score")
      abline(v=validation_DASL_TMZ_M_RNASeq[[2]][,1], col = "red")
      validation_DASL_TMZ_M_RNASeq[[3]][1:2,]

      
    }
    
    
  } # end of IV.1 - Validate by computing combination score

  
  
  ### ++ 2 - Validate with external list ----
  tmz_external_signature = c("ACSL3", "APOD", "BIRC5", "CDT1", "CHD7", "CHI3L1", "CYR61", "DDX56", "EIF2B5", "EPAS1", 
                             "FEN1", "FTSJ2", "GSTP1", "HEPN1", "HIST1H2AC", "ILF3", "LGALS3", "LST1", "MCM10", "MEX3B", 
                             "MIF", "MSANTD3.TMEFF1", "MT1X", "MTA2", "MTRF1", "MVP", "NCAM1", "NFIL3", "NOL4L", "OSMR", 
                             "PKN1", "PTBP1", "RCC2", "SH3BGRL3", "TPI1", "UCHL5", "VEZF1", "WT1")
  sum(signatures %in% tmz_external_signature)
  
} # end of IV - Validate signatures
