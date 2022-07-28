# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#                                  GLIOTRAIN Data Analysis: DEA
#
# Run DEA on RNA-Seq and GISTIC CNV data.
# Results are exported, submitted to enrichment analysis 
#    and tested against other datasets: RNA-Seq/CNV, logR values, EMC DASL data, TCGA datasets
#
# Started on the 11th of May, 2020
# 
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

### Environment set-up ----
setwd("/path/to/working_directory/")

library(biomaRt)
library(org.Hs.eg.db)
library(glmpca)
library(DESeq2)
library(ggplot2)
library(topGO)
library(MASS)


# oooooooooooooooooooooooooooooooooooooooooooooooooo
# I - Preparations ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

### I.1 Load and prepare data ----
{
  
  rna_batches = read.table("./raw/subset1_GTDATAreformatted_rnaseq_batches.tsv", sep = '\t', h = T)
  wgs_batches = read.table("./raw/subset1_GTDATAreformatted_wgs_batches.tsv", sep = '\t', h = T)
  
  
  ### Clinical data
  clinic_data = read.table("./raw/subset1_GTDATAreformatted_clinical_data.tsv", sep = '\t', h = T)
  clinic_data$Survivor_Type = ifelse(clinic_data$Overall_Survival_Overall_Survival_months <= 12, 
                                     "ST", 
                                     ifelse(clinic_data$Overall_Survival_Overall_Survival_months > 36,
                                            "LT", "IT"))
  
  ## Add batch information to clinic_data
  for(p in 1:nrow(clinic_data)){
    if(clinic_data$Subject_ID[p] %in% rna_batches$SubjectID){  clinic_data$rna_batch[p] = rna_batches$RNASeq_batches[match(clinic_data$Subject_ID[p], rna_batches$SubjectID)]  }
    if(clinic_data$Subject_ID[p] %in% wgs_batches$SubjectID){  clinic_data$wgs_batch[p] = wgs_batches$WGS_batches[match(clinic_data$Subject_ID[p], wgs_batches$SubjectID)]  }
  }
  
  # ## Separate IT from ST/LT
  # clinic_IT = clinic_data[clinic_data$Survivor_Type == "IT", ]
  # clinic_STLT = clinic_data[clinic_data$Survivor_Type != "IT", ]
  
  
  
  ### RNA-Seq data
  {
    ## Load data.
    rna_data_raw = read.table("./raw/subset1_GTDATAreformatted_RNASeq_rawCounts.tsv", sep = '\t', h = T)
    rna2genes = rna_data_raw[, c(1,2)]
    load("./tmp/rnaseq_norm_BatchCenterCP.RData")
    rnaseq_norm[1:5,1:5]
    
    ## Make adjustments.
    rownames(rna_data_raw) = rna_data_raw$Transcript_ID
    rna_data_raw = rna_data_raw[, -c(1,2)]
    colnames(rna_data_raw) = gsub("\\.\\d+\\.\\w+\\.\\d+\\.\\w+$", "", colnames(rna_data_raw))
    colnames(rna_data_raw) = gsub("\\.", "\\-", colnames(rna_data_raw))
    rna_data_vst = t(rnaseq_norm[-grep("CP", rownames(rnaseq_norm)), ])
    colnames(rna_data_vst) = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", colnames(rna_data_vst))
    
    remove(rnaseq_norm)
  }
  
  ### CNV data
  {
    # wgs_data_raw = read.table("./raw/subset1_GTDATAreformatted_CopyNumber_logRvalues.tsv", sep = '\t', h = T)
    # rownames(wgs_data_raw) = wgs_data_raw$Region
    # wgs_data_raw = wgs_data_raw[, -1]
    # colnames(wgs_data_raw) = gsub("\\.\\d+\\.\\w+\\.\\d+\\.\\w+$", "", colnames(wgs_data_raw))
    # colnames(wgs_data_raw) = gsub("\\.", "\\-", colnames(wgs_data_raw))
    
    wgs_data_cnv = read.table("./raw/cnv_data_focalEvents.tsv", sep = '\t', h = T)
    wgs2genes = read.table("./raw/logR2genes_map.tsv", sep = '\t', h = T)
    rownames(wgs_data_cnv) = wgs_data_cnv$region_name
    wgs_data_cnv = wgs_data_cnv[, -1]
    colnames(wgs_data_cnv) = gsub("\\.\\d+\\.\\w+\\.\\d+\\.\\w+\\.\\w+$", "", colnames(wgs_data_cnv))
    colnames(wgs_data_cnv) = gsub("\\.", "\\-", colnames(wgs_data_cnv))
    
    # ## Separate IT from ST/LT
    # wgs_cnv_IT = wgs_data_cnv[, colnames(wgs_data_cnv) %in% rownames(clinic_IT)]
    # wgs_cnv_STLT = wgs_data_cnv[, colnames(wgs_data_cnv) %in% rownames(clinic_STLT)]
  }
  
  
  ### RPPA data
  
  
  
  
  ### DASL data
  {
    load("../../EMCcollab/Ioannis_drugs/Data/tmp/samplesLists.RData")
    load("../../EMCcollab/Ioannis_drugs/Data/tmp/transc_preprocessed.RData")
    mgmt_dasl = read.table("../../EMCcollab/Ioannis_drugs/Data/raw/MGMT_status_GS_cell_cultures_IN.txt", 
                           h = TRUE, stringsAsFactors = FALSE)
    mgmt_tmz = read.csv("../../EMCcollab/TMZ_Lomustine/Data/tmp/TMZstudy_MGMTstatus.csv", 
                        h = TRUE, stringsAsFactors = FALSE)
    
    mgmt_dasl <- mgmt_dasl[na.omit(match(samples_allPrimaryGBM, mgmt_dasl$GS_number)),]
    clini_dasl <- diagnosis[match(samples_allPrimaryGBM, rownames(diagnosis)), ]
    clini_dasl$mgmt <- rep(NA, nrow(clini_dasl))
    clini_dasl$mgmt[match(mgmt_dasl$GS_number, rownames(clini_dasl))] <- mgmt_dasl$MGMT_status
    sum(!is.na(clini_dasl$mgmt))
    
    mgmt_tmz <- mgmt_tmz[na.omit(match(samples_allPrimaryGBM, mgmt_tmz$GS.number)),]
    which(mgmt_tmz$GS.number %in% rownames(clini_dasl[is.na(clini_dasl$mgmt),]))
    mgmt_tmz[c(7,10,14),]
    clini_dasl$mgmt[match(c("GS209", "GS370", "GS102"), rownames(clini_dasl))] = c("Unmethylated", "Methylated", "Methylated")
    clini_dasl["GS104c",] <- clini_dasl["GS104p",]
    clini_dasl["GS175", "OS"] <- 999 # from TMZstudy_survivalData file, still alive
    clini_dasl["GS175", "mgmt"] <- "Methylated" # from TMZstudy_survivalData file
    clini_dasl <- clini_dasl[!is.na(clini_dasl$OS), ]
    clini_dasl$Survivor_Type = ifelse(clini_dasl$OS <= 12, "ST", ifelse(clini_dasl$OS > 36, "LT", "IT"))
    
    ## Microarray data
    colnames(transc_preprocessed)[1:5]
    dasl_data = t(transc_preprocessed)
    emc_genes <- gsub("\\.\\d$", "", colnames(dasl_data))
    emc_genes[grep('\\.', emc_genes)]
    remove(transc_preprocessed)
    
    # ## Separate IT from ST/LT
    # clini_dasl_IT <- clini_dasl[clini_dasl$Survivor_Type == "IT", c("OS", "mgmt", "Survivor_Type")]
    # clini_dasl_STLT <- clini_dasl[clini_dasl$Survivor_Type != "IT", c("OS", "mgmt", "Survivor_Type")]
    # dasl_IT <- transc_preprocessed[ rownames(transc_preprocessed) %in% rownames(clini_dasl_IT), ]
    # dasl_STLT <- transc_preprocessed[ rownames(transc_preprocessed) %in% rownames(clini_dasl_STLT), ]
    
  }
  
  
  
  ### TCGA data
  {
    ## Clinical data.
    clini_tcga <- read.table("../../TCGA_GBM/2020.09.08/data_clinical_TCGA.txt", header = T)
    clini_tcga$Survivor_Type <- ifelse(clini_tcga$OS_MONTHS <= 12, "ST", ifelse(clini_tcga$OS_MONTHS > 36, "LT", "IT"))
    
    ## CNV data.
    cnv_tcga <- read.csv("../../TCGA_GBM/2020.09.08/data_CNA.txt", header = T, sep = '\t')
    colnames(cnv_tcga) <- gsub("\\.", "-", colnames(cnv_tcga))
    
    ## RNA-Seq data.
    rna_zscores_tcga <- read.table("../../TCGA_GBM/2020.09.08/data_RNA_Seq_v2_mRNA_median_Zscores.txt", header = T, sep = '\t')
    colnames(rna_zscores_tcga) <- gsub("\\.", "-", colnames(rna_zscores_tcga))
    
    # ## Separate IT from ST/LT
    # clini_tcga_IT <- clini_tcga[clini_tcga$Survivor_Type == "IT", ]
    # clini_tcga_STLT <- clini_tcga[clini_tcga$Survivor_Type != "IT", ]
    # cnv_tcga_IT <- cnv_tcga[, c(1,2,3, which(colnames(cnv_tcga) %in% clini_tcga_IT$SAMPLE_ID))]
    # cnv_tcga_STLT <- cnv_tcga[, c(1,2,3, which(colnames(cnv_tcga) %in% clini_tcga_STLT$SAMPLE_ID))]
    
  }
  
  
  
  ### Fetch ENSEMBL db to determine map position, ENSEMBL, ENTREZ and HUGO gene IDs.
  ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "89")
  genesMap <- getBM(mart = ensembl,
                    attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene", "external_gene_name", "chromosome_name", "start_position", "end_position"))
  genesMap <- genesMap[genesMap$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
                                                       "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
  head(genesMap)
  
  
  
  ### Finish.
  {
    ## Check out distributions.
    sapply(list(clinic_data, clini_dasl, clini_tcga), dim)
    sapply(list(clinic_data$Survivor_Type, clini_dasl$Survivor_Type, clini_tcga$Survivor_Type), table)
    sapply(list(rna_data_raw, rna_data_vst, wgs_data_cnv, dasl_data, rna_zscores_tcga, cnv_tcga), dim)
    sapply(list(rna_data_raw, rna_data_vst, wgs_data_cnv, dasl_data, rna_zscores_tcga, cnv_tcga), function(d) d[1:3,1:3])
    sort(table(wgs2genes$chr))
    
    save(clinic_data, clini_dasl, clini_tcga,
         rna_data_raw, rna_data_vst, wgs_data_cnv, dasl_data, rna_zscores_tcga, cnv_tcga,
         wgs2genes, rna_batches, rna_batches, genesMap,
         file = "./tmp/DEA_all_datasets.RData")
  }
  

  
}# end I.1 Load and prepare data


### I.2 Validation Functions ----
{
  
  ### The functions here will be used to see if a similar trend as what is observed in the results of a DEA 
  ###    can be obseved in other datasets. There is are two distinct functions per dataset: one is to cross-
  ###    validate using ST/LT comparison, the other is to cross-validate using IT group with a continuous 
  ###    response variable.
  ### Since they all have similar input and output, but work slightly differently based on the type of data
  ###    and response variable, the following "header" can be considered to be shared across all of them.
  
  
  ## Function to search the $cross-validating$ dataset for the significantly DE genes identified in a DE analysis.
  # Tries to find the DE genes in the $cross-validating$ dataset. 
  #   If it is found, runs a test to see if the $cross-validating$ data also reflects a difference between groups. 
  #     If it does, checks if the trend (monotonous direction) matches that of the provided DE genes.
  #
  # WARNING: assumes ordinal levels are properly defined.
  #
  # @param <DEres> a dataset resulting from a DE analysis. Must contain a "geneName" column, 
  #                  and a "difference" numeric column giving amplitude and directionality of the groups 
  #                  (e.g. foldchange, correlation...)
  # @param <genesMap> a biomart data to identify the locus on which the gene is, based on ENSEMBL ID (for RNA-Seq and CNV data)
  # @param <validationData> the $cross-validating$ dataset
  # @param <CNVregions> a 3-columns ("chromosome", "start", "end") dataframe with rownames matching that
  #                     of <validationData> to provide loci information.
  # @param <responseVariable> the response variable for the test of the <validationData>. 
  #                              Must be a vector named in the same manner as <validationData>
  # @return <validationOutput> a dataframe indicating for each gene in <DEres> the p-value of the test
  #                             on <validationData> being DE as well (-1 if gene not in the data, -2 if gene is not 
  #                             recognized at all, -3 if no variation across samples), 
  #                             and if that DE is in the same direction as in the DEA dataset
  
  
  validate_with_GT.CNV.STLT <- function(DEres, genesMap, validationData, CNVregions, responseVariable){
    cat("#### COMPARING RESULTS WITH GT CNV ST vs LT ####\n")
    
    # Extract ENSEMBL ID (remove the version '.xx' extension in rownames).
    #ensID = gsub("\\.\\d+", "", rownames(DEres))
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)),
                                   matching.region = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      
      if(is.na(posInMart)){  
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding region in cnv.
        cnvReg = which(CNVregions$chromosome == paste0("chr", genesMap$chromosome_name[posInMart]) & 
                         CNVregions$start <= genesMap$start_position[posInMart] & 
                         CNVregions$end >= genesMap$end_position[posInMart])
        
        if(isEmpty(cnvReg)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          cat(rownames(DEres)[g], " was found in the CNV data at ", cnvReg, "\n")
          if(length(cnvReg) > 1){ cnvReg = min(cnvReg)  } # TODO: find a better way to handle this
          validationOutput[g, "matching.region"] = rownames(CNVregions)[cnvReg]
          
          # Create a dataset that combines...
          responseVariable <- responseVariable[match(colnames(validationData), names(responseVariable))]
          dat <- data.frame(cnv = as.numeric(validationData[cnvReg, ]), # ... the cnv data for the corresponding region
                            response = responseVariable, # ... the response variable
                            response_num = ifelse(responseVariable == "ST", 0, 1)) # ... a numeric conversion of the response variable
          
          # Test if CNV is dependent on <responseVariable>.
          #print(as.numeric(dat$cnv))
          #print(dat)
          kTest <- kruskal.test(dat$cnv, dat$response)
          validationOutput[g, "test.pvalue"] = kTest$p.value
          
          if(kTest$p.value < 0.05){
            
            # Test if direction of CNV is the same as DE.
            spearCor = cor.test(dat$cnv, dat$response_num, method = "spearman")
            direction = ifelse( (spearCor$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
            validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(spearCor$estimate, 3), 
                                                    ", pval=", round(spearCor$p.value, 4))
            
          }
          
        }# end of "gene was found in the CNV data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  }# end of validate_with_GT.CNV.STLT()
  
  validate_with_GT.CNV.IT <- function(DEres, genesMap, validationData, CNVregions, responseVariable){
    cat("#### COMPARING RESULTS WITH GT CNV IT ####\n")
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)),
                                   matching.region = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      
      if(is.na(posInMart)){  
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding region in cnv.
        cnvReg = which(CNVregions$chromosome == paste0("chr", genesMap$chromosome_name[posInMart]) & 
                         CNVregions$start <= genesMap$start_position[posInMart] & 
                         CNVregions$end >= genesMap$end_position[posInMart])
        
        if(isEmpty(cnvReg)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          cat(rownames(DEres)[g], " was found in the CNV data at ", cnvReg, "\n")
          if(length(cnvReg) > 1){ cnvReg = min(cnvReg)  } # TODO: find a better way to handle this
          validationOutput[g, "matching.region"] = rownames(CNVregions)[cnvReg]
          
          # Create a dataset that combines...
          dat <- data.frame(cnv = as.numeric(validationData[cnvReg, ]), # ... the cnv data for the corresponding region
                            response = responseVariable[match(colnames(validationData), names(responseVariable))])  # ... the response variable
          
          # Test if CNV is dependent on <responseVariable>.
          corTest <- cor.test(dat$cnv, dat$response, method = "spearman")
          validationOutput[g, "test.pvalue"] = corTest$p.value
          
          if(corTest$p.value < 0.05){
            
            # Test if direction of CNV is the same as DE.
            direction = ifelse( (corTest$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
            validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(corTest$estimate, 3), 
                                                    ", pval=", round(corTest$p.value, 4))
            
          }
          
        }# end of "gene was found in the CNV data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  }# end of validate_with_GT.CNV.IT()
  
  validate_with_GT.RNASeq.STLT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH GT RNASeq ST vs LT ####\n")
    
    # Extract ENSEMBL ID (remove the version '.xx' extension in rownames).
    #ensID = gsub("\\.\\d+", "", rownames(DEres))
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      
      if(is.na(posInMart)){    
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding transcript in RNA-Seq data
        ensID = gsub("\\.\\d+$", "", rownames(validationData))
        transcriptIdx = which(ensID == genesMap$ensembl_gene_id[posInMart])
        
        if(isEmpty(transcriptIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          cat(rownames(DEres)[g], " was found in the RNA-Seq data at ", transcriptIdx, "\n")
          if(length(transcriptIdx) > 1){  transcriptIdx = min(transcriptIdx)  } # TODO: find a better way to handle this
          
          if(sd(as.numeric(validationData[transcriptIdx, -c(1,2)])) == 0){ # if gene has no variation across patients
            validationOutput[g, "test.pvalue"] = -3
          }else{
            # Create a dataset that combines...
            responseVariable = responseVariable[match(colnames(validationData), names(responseVariable))]
            dat <- data.frame(rnaseq = as.numeric(validationData[transcriptIdx, ]), # ... the RNA-Seq data for the corresponding region
                              response = responseVariable[match(colnames(validationData), names(responseVariable))],  # ... the response variable
                              response_num = ifelse(responseVariable == "ST", 0, 1)) # ... a numeric conversion of the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            tTest <- t.test(dat$rnaseq[dat$response == "ST"], dat$rnaseq[dat$response == "LT"])
            validationOutput[g, "test.pvalue"] = tTest$p.value
            
            if(tTest$p.value < 0.05){
              
              # Test if direction of CNV is the same as DE.
              spearCor = cor.test(dat$rnaseq, dat$response_num, method = "spearman")
              direction = ifelse( (spearCor$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(spearCor$estimate, 3), 
                                                      ", pval=", round(spearCor$p.value, 4))
              
            }
          }
        }# end of "gene was found in the RNA-Seq data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  }# end of validate_with_GT.RNASeq.STLT()
  
  validate_with_GT.RNASeq.IT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH GT RNA-Seq IT ####\n")
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){   
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding transcript in RNA-Seq data
        ensID = gsub("\\.\\d+$", "", rownames(validationData))
        transcriptIdx = which(ensID == genesMap$ensembl_gene_id[posInMart])
        
        if(isEmpty(transcriptIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          cat(rownames(DEres)[g], " was found in the RNA-Seq data at ", transcriptIdx, "\n")
          if(length(transcriptIdx) > 1){ transcriptIdx = min(transcriptIdx)  } # TODO: find a better way to handle this
          
          if(sd(as.numeric(validationData[transcriptIdx, -c(1,2)])) == 0){ # if gene has no variation across patients
            validationOutput[g, "test.pvalue"] = -3
          }else{
            # Create a dataset that combines...
            dat <- data.frame(rnaseq = as.numeric(validationData[transcriptIdx, ]), # ... the RNA-Seq data for the corresponding region
                              response = responseVariable[match(colnames(validationData), names(responseVariable))])  # ... the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            corTest <- cor.test(dat$rnaseq, dat$response, method = "pearson")
            validationOutput[g, "test.pvalue"] = corTest$p.value
            
            if(corTest$p.value < 0.05){
              
              # Test if direction of CNV is the same as DE.
              direction = ifelse( (corTest$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(corTest$estimate, 3), 
                                                      ", pval=", round(corTest$p.value, 4))
            }
          }
        }# end of "gene was found in the RNA-Seq data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  }# end of validate_with_GT.RNASeq.IT()
  
  validate_with_GT.RPPA.STLT <- function(DEres, genesMap, validationData, CNVregions, responseVariable){  } # end of validate_with_GT.RPPA.STLT()
  
  validate_with_GT.RPPA.IT <- function(DEres, genesMap, validationData, CNVregions, responseVariable){  } # end of validate_with_GT.RPPA.IT()
  
  validate_with_DASL.STLT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH DASL ST vs LT ####\n")
    
    
    validationGenes <- gsub("\\.\\d+$", "", rownames(validationData))
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){  
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding transcript in DASL data
        geneIdx = which(validationGenes %in% c(genesMap$hgnc_symbol[posInMart], genesMap$external_gene_name[posInMart]))
        
        if(isEmpty(geneIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          cat(rownames(DEres)[g], " was found in the DASL data at ", geneIdx, "\n")
          
          for(probe in 1:length(geneIdx)){
            
            # Create a dataset that combines...
            responseVariable = responseVariable[match(colnames(validationData), names(responseVariable))]
            dat <- data.frame(dasl = as.numeric(validationData[geneIdx[probe] , ]), # ... the RNA-Seq data for the corresponding region
                              response = responseVariable,  # ... the response variable
                              response_num = ifelse(responseVariable == "ST", 0, 1)) # ... a numeric conversion of the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            tTest <- t.test(dat$dasl[dat$response == "ST"], dat$dasl[dat$response == "LT"])
            validationOutput[g, "test.pvalue"] = ifelse(probe == 1, tTest$p.value, 
                                                        paste0(validationOutput[g, "test.pvalue"], " ; ", tTest$p.value))
            
            if(tTest$p.value < 0.05){
              
              # Test if direction of validation data is the same as DE.
              spearCor = cor.test(dat$dasl, dat$response_num, method = "spearman")
              direction = ifelse( (spearCor$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              directionToPrint = paste0(direction, " direction: rho=", round(spearCor$estimate, 3), ", pval=", round(spearCor$p.value, 4))
              validationOutput[g, "comment"] = ifelse(probe == 1, directionToPrint, 
                                                      paste0(validationOutput[g, "comment"], " ; ", directionToPrint))
              
            }
            
          }
          
        }# end of "gene was found in the DASL data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  } # end of validate_with_DASL.STLT()
  
  validate_with_DASL.IT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH DASL IT ####\n")
    
    
    validationGenes <- gsub("\\.\\d+$", "", rownames(validationData))
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){   
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding transcript in DASL data
        geneIdx = which(validationGenes %in% c(genesMap$hgnc_symbol[posInMart], genesMap$external_gene_name[posInMart]))
        
        if(isEmpty(geneIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          cat(rownames(DEres)[g], " was found in the DASL data at ", geneIdx, "\n")
          
          for(probe in 1:length(geneIdx)){
            
            # Create a dataset that combines...
            dat <- data.frame(dasl = as.numeric(validationData[ geneIdx[probe], ]), # ... the RNA-Seq data for the corresponding region
                              response = responseVariable[match(colnames(validationData), names(responseVariable))])  # ... the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            corTest <- cor.test(dat$dasl, dat$response, method = "pearson")
            validationOutput[g, "test.pvalue"] = ifelse(probe == 1, corTest$p.value, 
                                                        paste0(validationOutput[g, "test.pvalue"], " ; ", corTest$p.value))
            
            if(corTest$p.value < 0.05){
              
              # Test if direction of validation data is the same as DE.
              direction = ifelse( (corTest$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              directionToPrint = paste0(direction, " direction: rho=", round(corTest$estimate, 3), ", pval=", round(corTest$p.value, 4))
              validationOutput[g, "comment"] = ifelse(probe == 1, directionToPrint, 
                                                      paste0(validationOutput[g, "comment"], " ; ", directionToPrint))
              
            }
            
          }
          
        }# end of "gene was found in the DASL data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  } # end of validate_with_DASL.IT()
  
  validate_with_TCGA.RNASeq.STLT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH TCGA RNASeq ST vs LT ####\n")
    
    # Extract ENSEMBL ID (remove the version '.xx' extension in rownames).
    #ensID = gsub("\\.\\d+", "", rownames(DEres))
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){   
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding transcript in RNA-Seq data
        transcriptIdx = which(validationData$Hugo_Symbol %in% c(genesMap$hgnc_symbol[posInMart], genesMap$external_gene_name[posInMart]) |
                                validationData$Entrez_Gene_Id == genesMap$entrezgene[posInMart])
        
        if(isEmpty(transcriptIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          #cat(rownames(DEres)[g], " was found in the RNA-Seq data at ", transcriptIdx, "\n")
          if(length(transcriptIdx) > 1){  transcriptIdx = min(transcriptIdx)  } # TODO: find a better way to handle this
          
          if(sd(as.numeric(validationData[transcriptIdx, -c(1,2)])) == 0){ # if gene has no variation across patients
            validationOutput[g, "test.pvalue"] = -3
          }else{
            # Create a dataset that combines...
            responseVariable = responseVariable[match(colnames(validationData)[-c(1,2)], names(responseVariable))]
            dat <- data.frame(rnaseq = as.numeric(validationData[transcriptIdx, -c(1,2)]), # ... the RNA-Seq data for the corresponding region
                              response = responseVariable,  # ... the response variable
                              response_num = ifelse(responseVariable == "ST", 0, 1)) # ... a numeric conversion of the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            tTest <- t.test(dat$rnaseq[dat$response == "ST"], dat$rnaseq[dat$response == "LT"])
            validationOutput[g, "test.pvalue"] = tTest$p.value
            
            if(tTest$p.value < 0.05){
              
              # Test if direction of CNV is the same as DE.
              spearCor = cor.test(dat$rnaseq, dat$response_num, method = "spearman")
              direction = ifelse( (spearCor$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(spearCor$estimate, 3), 
                                                      ", pval=", round(spearCor$p.value, 4))
              
            }
          }
        }# end of "gene was found in the RNA-Seq data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  }# end of validate_with_TCGA.RNASeq.STLT()
  
  validate_with_TCGA.RNASeq.IT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH TCGA RNA-Seq IT ####\n")
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){  
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding transcript in RNA-Seq data
        transcriptIdx = which(validationData$Hugo_Symbol %in% c(genesMap$hgnc_symbol[posInMart], genesMap$external_gene_name[posInMart]) |
                                validationData$Entrez_Gene_Id == genesMap$entrezgene[posInMart])
        
        if(isEmpty(transcriptIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          #cat(rownames(DEres)[g], " was found in the RNA-Seq data at ", transcriptIdx, "\n")
          if(length(transcriptIdx) > 1){ transcriptIdx = min(transcriptIdx)  } # TODO: find a better way to handle this
          
          if(sd(as.numeric(validationData[transcriptIdx, -c(1,2)])) == 0){ # if gene has no variation across patients
            validationOutput[g, "test.pvalue"] = -3
          }else{
            # Create a dataset that combines...
            dat <- data.frame(rnaseq = as.numeric(validationData[transcriptIdx, -c(1,2)]), # ... the RNA-Seq data for the corresponding region
                              response = responseVariable[match(colnames(validationData)[-c(1,2)], names(responseVariable))])  # ... the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            corTest <- cor.test(dat$rnaseq, dat$response, method = "pearson")
            validationOutput[g, "test.pvalue"] = corTest$p.value
            
            if(corTest$p.value < 0.05){
              
              # Test if direction of CNV is the same as DE.
              direction = ifelse( (corTest$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(corTest$estimate, 3), 
                                                      ", pval=", round(corTest$p.value, 4))
            }
          }
          
        }# end of "gene was found in the RNA-Seq data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  }# end of validate_with_TCGA.RNASeq.IT()
  
  validate_with_TCGA.CNV.STLT <- function(DEres, genesMap, validationData, responseVariable){
    cat("#### COMPARING RESULTS WITH TCGA CNV ST vs LT ####\n")
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){  
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding gene in TCGA CNV data
        geneIdx = which(validationData$Hugo_Symbol %in% c(genesMap$hgnc_symbol[posInMart], genesMap$external_gene_name[posInMart]) |
                          validationData$Entrez_Gene_Id == genesMap$entrezgene[posInMart])
        
        if(isEmpty(geneIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          #cat(rownames(DEres)[g], " was found in the TCGA CNV data at ", geneIdx, "\n")
          
          for(probe in 1:length(geneIdx)){
            
            # Create a dataset that combines...
            responseVariable = responseVariable[match(colnames(validationData)[-c(1,2,3)], names(responseVariable))]
            dat <- data.frame(cnv = as.numeric(validationData[geneIdx[probe], -c(1,2,3)]), # ... the CNV data for the corresponding gene
                              response = responseVariable,  # ... the response variable
                              response_num = ifelse(responseVariable == "ST", 0, 1)) # ... a numeric conversion of the response variable
            
            # Test if validation data is dependent on <responseVariable>.
            kTest <- kruskal.test(dat$cnv, dat$response)
            validationOutput[g, "test.pvalue"] = kTest$p.value
            
            if(kTest$p.value < 0.05){
              
              # Test if direction of CNV is the same as DE.
              spearCor = cor.test(dat$cnv, dat$response_num, method = "spearman")
              direction = ifelse( (spearCor$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              directionToPrint = paste0(direction, " direction: rho=", round(spearCor$estimate, 3), ", pval=", round(spearCor$p.value, 4))
              validationOutput[g, "comment"] = ifelse(probe == 1, directionToPrint, 
                                                      paste0(validationOutput[g, "comment"], " ; ", directionToPrint))
              
            }
            
          }
          
        }# end of "gene was found in the TCGA CNV data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  } # end of validate_with_TCGA.CNV.STLT()
  
  validate_with_TCGA.CNV.IT <- function(DEres, genesMap, validationData, CNVregions, responseVariable){
    cat("#### COMPARING RESULTS WITH TCGA CNV IT ####\n")
    
    # Create an empty table which will be filled with the results.
    validationOutput <- data.frame(test.pvalue = rep(NA, nrow(DEres)),
                                   comment = rep(NA, nrow(DEres)))
    rownames(validationOutput) <- rownames(DEres)
    
    # Run the test for each gene in <DEres>.
    cat("** Testing for each DE gene **\n")
    for(g in 1:nrow(DEres)){ 
      
      gName = DEres$geneName[g]
      posInMart = ifelse(gName %in% genesMap$external_gene_name, min(which(genesMap$external_gene_name == gName)),
                         ifelse(gName %in% genesMap$ensembl_gene_id, min(which(genesMap$ensembl_gene_id == gName)),
                                ifelse(gName %in% genesMap$hgnc_symbol, min(which(genesMap$hgnc_symbol == gName)), NA))) # find entry in biomart
      if(is.na(posInMart)){   
        cat(rownames(DEres)[g], " not found in BioMart\n")  
        validationOutput[g, "test.pvalue"] = -2
      }else{
        
        # Find corresponding gene in TCGA CNV data
        geneIdx = which(validationData$Hugo_Symbol %in% c(genesMap$hgnc_symbol[posInMart], genesMap$external_gene_name[posInMart]) |
                          validationData$Entrez_Gene_Id == genesMap$entrezgene[posInMart])
        
        if(isEmpty(geneIdx)){  
          #cat(rownames(DEres)[g], " not represented in the CNV data\n")  
          validationOutput[g, "test.pvalue"] = -1
        }else{
          #cat(rownames(DEres)[g], " was found in the TCGA CNV data at ", geneIdx, "\n")
          
          for(probe in 1:length(geneIdx)){
            
            # Create a dataset that combines...
            dat <- data.frame(cnv = as.numeric(validationData[ geneIdx[probe], -c(1,2,3)]), # ... the CNV data for the corresponding gene
                              response = responseVariable[match(colnames(validationData)[-c(1,2,3)], names(responseVariable))])  # ... the response variable
            
            
            # Test if CNV is dependent on <responseVariable>.
            corTest <- cor.test(dat$cnv, dat$response, method = "spearman")
            validationOutput[g, "test.pvalue"] = corTest$p.value
            
            if(corTest$p.value < 0.05){
              
              # Test if direction of CNV is the same as DE.
              direction = ifelse( (corTest$estimate * DEres[g, "difference"]) < 0, "Opposite", "Same")
              validationOutput[g, "comment"] = paste0(direction, " direction: rho=", round(corTest$estimate, 3), 
                                                      ", pval=", round(corTest$p.value, 4))
              
            }
            
          }
          
        }# end of "gene was found in the TCGA CNV data\n"
      }# end of "gene was found in biomart\n"
      
      
    }# end of for loop
    
    return(validationOutput)
    
    
  } # end of validate_with_TCGA.CNV.IT()
  
} # end of I.2 Functions

### I.3 Other Functions ----
{
  
  runTopGO <- function(geneUniverse, investigatedGenes, outPath, title="TopGo_analysis", pvalues = NULL, ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 10){
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
  
  
}


# oooooooooooooooooooooooooooooooooooooooooooooooooo
# II - DE Analysis for RNA-Seq Data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

{
# # ===== II.1 Run analysis 
# 
# # ST vs LT 
# {
#   ## Prepare objects
#   coldata <- clinic_data[clinic_data$Subject_ID %in% colnames(rna_data), c("Subject_ID", "Survivor_Type", "rna_batch", "Centers")]
#   coldata_de = coldata[coldata$Survivor_Type != "IT", ]
#   coldata_de$rna_batch <- as.factor(coldata_de$rna_batch)
#   rna_data_de = rna_data[, match(coldata_de$Subject_ID, colnames(rna_data))]
#   
#   ddsMat_STvsLT <- DESeqDataSetFromMatrix(countData = rna_data_de, colData = coldata_de, design = ~ rna_batch + Centers + Survivor_Type)
#   dim(ddsMat_STvsLT)
#   ddsMat_STvsLT <- ddsMat_STvsLT[rowSums(counts(ddsMat_STvsLT)) > 10, ] # filter out transcripts that have less than 10 count total
#   dim(ddsMat_STvsLT)
#   
#   ## Run DE analysis
#   dds_STvsLT <- DESeq(ddsMat_STvsLT)
#   res <- results(dds_STvsLT, contrast = c("Survivor_Type", "ST", "LT"), alpha = 0.05)
#   sum(res$padj < 0.05, na.rm = TRUE)
#   hist(res$padj, breaks = 100)
#   summary(res)
#   
#   resultsNames(dds_STvsLT)
#   resLFC <- lfcShrink(dds_STvsLT, coef="Survivor_Type_ST_vs_LT", type="apeglm") # first, install the Bioconductor package 'apeglm'
#   resLFC
#   summary(resLFC, alpha = 0.05)
#   
#   resSignif_STvsLT <- res[which(res$padj < 0.05), ]
#   resSignif_STvsLT <- resSignif_STvsLT[order(resSignif_STvsLT$padj),]
#   resSignif_STvsLT$geneName <- rna2genes[match(rownames(resSignif_STvsLT), rna2genes$Transcript_ID), 2]
#   head(resSignif_STvsLT)
#   plotMA(resSignif_STvsLT, ylim=c(-2,2))
# } # deprecated
# 
# # ST vs LT accounting for MGMT 
# {
#   ## Prepare objects
#   coldata <- clinic_data[clinic_data$Subject_ID %in% colnames(rna_data), c("Subject_ID", "Survivor_Type", "MGMT_Methylation_Status", "rna_batch", "Centers")]
#   coldata_de = coldata[coldata$Survivor_Type != "IT", ]
#   coldata_de$rna_batch <- as.factor(coldata_de$rna_batch)
#   rna_data_de = rna_data[, match(coldata_de$Subject_ID, colnames(rna_data))]
#   
#   ddsMat_MGMT <- DESeqDataSetFromMatrix(countData = rna_data_de, colData = coldata_de, design = ~ rna_batch + Centers + MGMT_Methylation_Status + Survivor_Type + MGMT_Methylation_Status:Survivor_Type)
#   dim(ddsMat_MGMT)
#   ddsMat_MGMT <- ddsMat_MGMT[rowSums(counts(ddsMat_MGMT)) > 10, ] # filter out transcripts that have less than 10 count total
#   dim(ddsMat_MGMT)
#   
#   ## Run DE analysis
#   dds_MGMT <- DESeq(ddsMat_MGMT)
#   summary(dds_MGMT)
#   resSurv <- results(dds_MGMT, contrast = c("Survivor_Type", "ST", "LT"), alpha = 0.05)
#   resMGMT <- results(dds_MGMT, contrast = c("MGMT_Methylation_Status", "unmethylated", "methylated"), alpha = 0.05)
#   resultsNames(dds_MGMT)
#   resMGMTSurv <-results(dds_MGMT, name = "MGMT_Methylation_Statusunmethylated.Survivor_TypeST", alpha = 0.05)
#   
#   hist(resSurv$padj, breaks = 100)
#   sum(resSurv$padj < 0.05, na.rm = TRUE)
#   summary(resSurv)
#   resSignif_MGMT_Surv <- resSurv[which(resSurv$padj < 0.05), ]
#   resSignif_MGMT_Surv <- resSignif_MGMT_Surv[order(resSignif_MGMT_Surv$padj),]
#   
#   hist(resMGMT$padj, breaks = 100)
#   sum(resMGMT$padj < 0.05, na.rm = TRUE)
#   summary(resMGMT)
#   resSignif_MGMT_MGMT <- resMGMT[which(resMGMT$padj < 0.05), ]
#   resSignif_MGMT_MGMT <- resSignif_MGMT_MGMT[order(resSignif_MGMT_MGMT$padj),]
#   
#   hist(resMGMTSurv$padj, breaks = 100)
#   sum(resMGMTSurv$padj < 0.05, na.rm = TRUE)
#   summary(resMGMTSurv)
#   resSignif_MGMT_MGMTSurv <- resMGMTSurv[which(resMGMTSurv$padj < 0.05), ]
#   resSignif_MGMT_MGMTSurv <- resSignif_MGMT_MGMTSurv[order(resSignif_MGMT_MGMTSurv$padj),]
#   resSignif_MGMT_Surv$geneName <- rna2genes[match(rownames(resSignif_MGMT_Surv), rna2genes$Transcript_ID), 2]
#   resSignif_MGMT_MGMT$geneName <- rna2genes[match(rownames(resSignif_MGMT_MGMT), rna2genes$Transcript_ID), 2]
#   resSignif_MGMT_MGMTSurv$geneName <- rna2genes[match(rownames(resSignif_MGMT_MGMTSurv), rna2genes$Transcript_ID), 2]
#   
#   
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_Surv)[1], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_Surv)[2], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_Surv)[3], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMT)[1], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMT)[2], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMT)[3], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMTSurv)[1], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMTSurv)[2], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMTSurv)[3], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
#              normalized = TRUE, transform = TRUE, pch=20)
#   
#   plotMA(resSignif_MGMT_Surv, ylim=c(-2,2))
#   plotMA(resSignif_MGMT_MGMT, ylim=c(-2,2))
#   plotMA(resSignif_MGMT_MGMTSurv, ylim=c(-2,2))
#   
# }
# 
# # All Numeric OS accounting for MGMT 
# {
#   ## Prepare objects
#   coldata <- clinic_data[clinic_data$Subject_ID %in% colnames(rna_data), c("Subject_ID", "Overall_Survival_Overall_Survival_months", "MGMT_Methylation_Status", "rna_batch", "Centers")]
#   coldata_de = coldata
#   coldata_de$rna_batch <- as.factor(coldata_de$rna_batch)
#   rna_data_de = rna_data[, match(coldata_de$Subject_ID, colnames(rna_data))]
#   
#   ddsMat_OS <- DESeqDataSetFromMatrix(countData = rna_data_de, colData = coldata_de, design = ~ rna_batch + Centers + MGMT_Methylation_Status + Overall_Survival_Overall_Survival_months + MGMT_Methylation_Status:Overall_Survival_Overall_Survival_months)
#   dim(ddsMat_OS)
#   ddsMat_OS <- ddsMat_OS[rowSums(counts(ddsMat_OS)) > 10, ] # filter out transcripts that have less than 10 count total
#   dim(ddsMat_OS)
#   
#   ## Run DE analysis
#   dds_OS <- DESeq(ddsMat_OS)
#   summary(dds_OS)
#   resultsNames(dds_OS)
#   resSurv <- results(dds_OS, name = "Overall_Survival_Overall_Survival_months", alpha = 0.05)
#   resMGMT <- results(dds_OS, contrast = c("MGMT_Methylation_Status", "unmethylated", "methylated"), alpha = 0.05)
#   resMGMTSurv <- results(dds_OS, name = "MGMT_Methylation_Statusunmethylated.Overall_Survival_Overall_Survival_months", alpha = 0.05)
#   
#   hist(resSurv$padj, breaks = 100)
#   sum(resSurv$padj < 0.05, na.rm = TRUE)
#   summary(resSurv)
#   resSignif_OS_Surv <- resSurv[which(resSurv$padj < 0.05), ]
#   resSignif_OS_Surv <- resSignif_OS_Surv[order(resSignif_OS_Surv$padj),]
#   
#   hist(resMGMT$padj, breaks = 100)
#   sum(resMGMT$padj < 0.05, na.rm = TRUE)
#   summary(resMGMT)
#   resSignif_OS_MGMT <- resMGMT[which(resMGMT$padj < 0.05), ]
#   resSignif_OS_MGMT <- resSignif_OS_MGMT[order(resSignif_OS_MGMT$padj),]
#   
#   hist(resMGMTSurv$padj, breaks = 100)
#   sum(resMGMTSurv$padj < 0.05, na.rm = TRUE)
#   summary(resMGMTSurv)
#   resSignif_OS_MGMTSurv <- resMGMTSurv[which(resMGMTSurv$padj < 0.05), ]
#   resSignif_OS_MGMTSurv <- resSignif_OS_MGMTSurv[order(resSignif_OS_MGMTSurv$padj),]
#   
#   plotMA(resSignif_OS_Surv, ylim=c(-2,2))
#   abline(h = 0.2, lty = 2)
#   plotMA(resSignif_OS_MGMT, ylim=c(-2,2))
#   plotMA(resSignif_OS_MGMTSurv, ylim=c(-2,2))
#   abline(h = -0.9, lty = 2)
#   
#   sum(abs(resSignif_OS_Surv$log2FoldChange) > 0.2)
#   sum(abs(resSignif_OS_MGMTSurv$log2FoldChange) > 0.9)
#   
#   resSignif_OS_Surv <- resSignif_OS_Surv[abs(resSignif_OS_Surv$log2FoldChange) > 0.2, ]
#   resSignif_OS_MGMTSurv <- resSignif_OS_MGMTSurv[abs(resSignif_OS_MGMTSurv$log2FoldChange) > 0.9, ]
#   resSignif_OS_Surv$geneName <- rna2genes[match(rownames(resSignif_OS_Surv), rna2genes$Transcript_ID), 2]
#   resSignif_OS_MGMT$geneName <- rna2genes[match(rownames(resSignif_OS_MGMT), rna2genes$Transcript_ID), 2]
#   resSignif_OS_MGMTSurv$geneName <- rna2genes[match(rownames(resSignif_OS_MGMTSurv), rna2genes$Transcript_ID), 2]
#   
# }
# 
# 
# # ===== II.2 Explore results 
# 
# # Look at overlap 
# {
#   STvsLT_genes <- rownames(resSignif_STvsLT)
#   MGMT_genes <- unique(c(rownames(resSignif_MGMT_Surv), rownames(resSignif_MGMT_MGMT), rownames(resSignif_MGMT_MGMTSurv)))
#   OS_genes <- unique(c(rownames(resSignif_OS_Surv), rownames(resSignif_OS_MGMT), rownames(resSignif_OS_MGMTSurv)))
#   
#   intersect_categs = intersect(STvsLT_genes, MGMT_genes)
#   intersect_surviv = intersect(STvsLT_genes, OS_genes)
#   intersect_MGMTst = intersect(OS_genes, MGMT_genes)
#   
#   sapply(list(STvsLT_genes, MGMT_genes, OS_genes), length)
#   sapply(list(intersect_categs, intersect_surviv, intersect_MGMTst), length)
#   
#   length(unique(c(intersect_categs, intersect_surviv, intersect_MGMTst)))
# }
# 
# 
# par(mfrow = c(3,3))
# for(g in 1:9){  plotCounts(dds, gene = rownames(resSignif_STvsLT)[g], intgroup = "Survivor_Type")  }
# for(g in 10:18){  plotCounts(dds, gene = rownames(resSignif_STvsLT)[g], intgroup = "Survivor_Type")  }
# par(mfrow = c(1,1))
# 
# # Perform Enrichment Analysis 
# {
#   ### Run Analyses.
#   geneUniverse = data.frame(probe.name = rna2genes$Gene_Symbol, probe.symbol = rna2genes$Gene_Symbol)
#   
#   ## STvsLT.
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_STvsLT$geneName, pvalues = resSignif_STvsLT$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_STvsLT", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
#   
#   ## MGMT.
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_Surv$geneName, pvalues = resSignif_MGMT_Surv$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_MGMT_Surv", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_MGMT$geneName, pvalues = resSignif_MGMT_MGMT$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_MGMT_MGMT", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_MGMTSurv$geneName, pvalues = resSignif_MGMT_MGMTSurv$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_MGMT_MGMTSurv", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
#   
#   ## OS.
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_OS_Surv$geneName, pvalues = resSignif_OS_Surv$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_OS_Surv", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_OS_MGMT$geneName, pvalues = resSignif_OS_MGMT$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_OS_MGMT", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
#   runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_OS_MGMTSurv$geneName, pvalues = resSignif_OS_MGMTSurv$padj,
#            outPath = "./results/2020.09.08/", title="TopGo_analysis_DESeq2_OS_MGMTSurv", 
#            ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
# }
# 
# # ===== II.3 Compare to other datasets 
# 
# ### Compare to GT CNV. 
# {
# 
#   
#   # Break down the rownames of cnv data to get loci information.
#   cnvRegions <- as.data.frame(t(data.frame(sapply(rownames(wgs_data_cnv), function(reg) strsplit(reg, '_')))))
#   colnames(cnvRegions) <- c("chromosome", "start", "end")
#   cnvRegions$start <- as.numeric(as.character(cnvRegions$start))
#   cnvRegions$end <- as.numeric(as.character(cnvRegions$end))
#   
#   
#   
#   
#   # Define responses variables to use.
#   cliniDE_STLT = clinic_data[clinic_data$Survivor_Type %in% c("ST", "LT"), c("Subject_ID", "Survivor_Type", "MGMT_Methylation_Status")]
#   SVLT_Type = factor(cliniDE_STLT$Survivor_Type, levels = c("ST", "LT"), ordered = TRUE)
#   SVLT_MGMT = factor(cliniDE_STLT$MGMT_Methylation_Status, levels = c("unmethylated", "methylated"), ordered = TRUE)
#   SVLT_SVMT = factor(paste0(cliniDE_STLT$Survivor_Type, cliniDE_STLT$MGMT_Methylation_Status),
#                      levels = c("STunmethylated", "STmethylated", "LTunmethylated", "LTmethylated"), ordered = TRUE)
#   names(SVLT_Type) <- names(SVLT_MGMT) <- names(SVLT_SVMT) <- cliniDE_STLT$Subject_ID
#   OS_all = clinic_data$Overall_Survival_Overall_Survival_months
#   OS_all_MGMT = clinic_data$MGMT_Methylation_Status
#   names(OS_all) <- names(OS_all_MGMT) <- clinic_data$Subject_ID
#   OS_all = OS_all[match(colnames(wgs_data_cnv), names(OS_all))]
#   OS_all_MGMT = OS_all_MGMT[match(colnames(wgs_data_cnv), names(OS_all_MGMT))]
#   
#   validate_GTCNV_STvsLT = validateWithCNV(DEres = resSignif_STvsLT, genesMap = genesMap, 
#                                           CNVdata = wgs_data_cnv[, match(names(SVLT_Type), colnames(wgs_data_cnv))],
#                                           CNVregions = cnvRegions, responseVariable = SVLT_Type)
#   
#   validate_GTCNV_MGMT_Surv = validateWithCNV(DEres = resSignif_MGMT_Surv, genesMap = genesMap, 
#                                              CNVdata = wgs_data_cnv[, match(names(SVLT_Type), colnames(wgs_data_cnv))],
#                                              CNVregions = cnvRegions, responseVariable = SVLT_Type)
#   validate_GTCNV_MGMT_MGMT = validateWithCNV(DEres = resSignif_MGMT_MGMT, genesMap = genesMap, 
#                                              CNVdata = wgs_data_cnv[, match(names(SVLT_MGMT), colnames(wgs_data_cnv))],
#                                              CNVregions = cnvRegions, responseVariable = SVLT_MGMT)
#   validate_GTCNV_MGMT_MGMTSurv = validateWithCNV(DEres = resSignif_MGMT_MGMTSurv, genesMap = genesMap, 
#                                                  CNVdata = wgs_data_cnv[, match(names(SVLT_SVMT), colnames(wgs_data_cnv))],
#                                                  CNVregions = cnvRegions, responseVariable = SVLT_SVMT)
#   
#   validate_GTCNV_OS_Surv = validateWithCNV(DEres = resSignif_OS_Surv, genesMap = genesMap, 
#                                            CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#                                            CNVregions = cnvRegions, responseVariable = OS_all)
#   validate_GTCNV_OS_MGMTSurv = validateWithCNV(DEres = resSignif_OS_MGMTSurv, genesMap = genesMap, 
#                                                CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#                                                CNVregions = cnvRegions, responseVariable = OS_all)
#   validate_GTCNV_OS_MGMT = validateWithCNV(DEres = resSignif_OS_MGMT, genesMap = genesMap, 
#                                            CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#                                            CNVregions = cnvRegions, responseVariable = OS_all_MGMT)
#   
#   validate_GTCNV_sets = list(validate_GTCNV_STvsLT,
#                              validate_GTCNV_MGMT_Surv, validate_GTCNV_MGMT_MGMT, validate_GTCNV_MGMT_MGMTSurv,
#                              validate_GTCNV_OS_Surv, validate_GTCNV_OS_MGMTSurv, validate_GTCNV_OS_MGMT)
#   names(validate_GTCNV_sets) = c("validate_GTCNV_STvsLT",
#                                  "validate_GTCNV_MGMT_Surv", "validate_GTCNV_MGMT_MGMT", "validate_GTCNV_MGMT_MGMTSurv",
#                                  "validate_GTCNV_OS_Surv", "validate_GTCNV_OS_MGMTSurv", "validate_GTCNV_OS_MGMT")
#   lapply(validate_GTCNV_sets, function(d) (d[d$kruskal.pvalue != -1,]))
#   lapply(validate_GTCNV_sets, function(d) paste(sum(d$kruskal.pvalue != -1), "/", nrow(d)))
#   
#   
#   resSignif_MGMT_Surv["ENSG00000204666.3",]
#   gloup = wgs_data_cnv["chr19_48225001_54624999", match(names(SVLT_MGMT), colnames(wgs_data_cnv))]
#   cbind(cliniDE_STLT, as.numeric(gloup))[order(cliniDE_STLT$Survivor_Type),]
#   
# }# end of ### Compare to GT CNV.
# 
# ### Compare to EMC DASL. 
# 
# ## Prepare data.
# 
# 
# 
# ## Compare survival-related genes.
# {
#   match(resSignif_MGMT_Surv$geneName, emc_genes)
#   cat(sum(is.na(match(resSignif_MGMT_Surv$geneName, emc_genes))), nrow(resSignif_MGMT_Surv))
#   
#   dasl_surv <- transc_preprocessed[match(rownames(clini_dasl), rownames(transc_preprocessed)), 
#                                    na.omit(match(resSignif_MGMT_Surv$geneName, emc_genes))]
#   daslValidate_surv <- data.frame(ttest_pvalue = rep(NA, ncol(dasl_surv)), 
#                                   comment = rep(NA, ncol(dasl_surv)), 
#                                   gene = colnames(dasl_surv))
#   rownames(daslValidate_surv) <- rownames(resSignif_MGMT_Surv)[which(!is.na(match(resSignif_MGMT_Surv$geneName, emc_genes)))]
#   daslValidate_surv$lChange_rnaSeq <- resSignif_MGMT_Surv$log2FoldChange[match(rownames(daslValidate_surv), rownames(resSignif_MGMT_Surv))]
#   for(g in 1:ncol(dasl_surv)){
#     dat = data.frame(survivorType = clini_dasl$Survivor_Type, expr = dasl_surv[, g])
#     ttest = t.test(dat[dat$survivorType=="ST", "expr"], dat[dat$survivorType=="LT", "expr"])
#     daslValidate_surv$ttest_pvalue[g] = ttest$p.value
#     if(ttest$p.value < 0.05){
#       diff = mean(dat[dat$survivorType=="LT", "expr"]) - mean(dat[dat$survivorType=="ST", "expr"])
#       if((diff * resSignif_MGMT_Surv[rownames(daslValidate_surv)[g], "log2FoldChange"]) < 0 ){
#         daslValidate_surv$comment[g] = paste("Opposite direction:", diff)
#       }
#       
#     }
#   }
# }
# 
# ## Compare MGMT-related genes.
# {
#   match(resSignif_MGMT_MGMT$geneName, emc_genes)
#   cat(sum(is.na(match(resSignif_MGMT_MGMT$geneName, emc_genes))), nrow(resSignif_MGMT_MGMT))
#   
#   dasl_mgmt <- transc_preprocessed[match(rownames(clini_dasl), rownames(transc_preprocessed)), 
#                                    na.omit(match(resSignif_MGMT_MGMT$geneName, emc_genes))]
#   daslValidate_mgmt <- data.frame(ttest_pvalue = rep(NA, ncol(dasl_mgmt)), 
#                                   comment = rep(NA, ncol(dasl_mgmt)), 
#                                   gene = colnames(dasl_mgmt))
#   rownames(daslValidate_mgmt) <- rownames(resSignif_MGMT_MGMT)[which(!is.na(match(resSignif_MGMT_MGMT$geneName, emc_genes)))]
#   daslValidate_mgmt$lChange_rnaSeq <- resSignif_MGMT_MGMT$log2FoldChange[match(rownames(daslValidate_mgmt), rownames(resSignif_MGMT_MGMT))]
#   for(g in 1:ncol(dasl_mgmt)){
#     dat = data.frame(survivorType = clini_dasl$Survivor_Type, expr = dasl_mgmt[, g])
#     ttest = t.test(dat[dat$survivorType=="ST", "expr"], dat[dat$survivorType=="LT", "expr"])
#     daslValidate_mgmt$ttest_pvalue[g] = ttest$p.value
#     if(ttest$p.value < 0.05){
#       diff = mean(dat[dat$survivorType=="LT", "expr"]) - mean(dat[dat$survivorType=="ST", "expr"])
#       if((diff * resSignif_MGMT_Surv[rownames(daslValidate_mgmt)[g], "log2FoldChange"]) < 0 ){
#         daslValidate_mgmt$comment[g] = paste("Opposite direction:", diff)
#       }
#       
#     }
#   }
# }
# 
# ### Compare to TCGA CNV. 
# {
#   ## Load data
#   
#   
#   sum(is.na(match(resSignif_MGMT_MGMT$geneName, cnv_tcga$Hugo_Symbol)))
#   resSignif_MGMT_MGMT$geneName[is.na(match(resSignif_MGMT_MGMT$geneName, cnv_tcga$Hugo_Symbol))]
#   sum(is.na(match(resSignif_MGMT_Surv$geneName, cnv_tcga$Hugo_Symbol)))
#   resSignif_MGMT_Surv$geneName[is.na(match(resSignif_MGMT_Surv$geneName, cnv_tcga$Hugo_Symbol))]
#   
#   
#   
#   ## Function to search TCGA CNV data for the significantly DE genes from GT RNA-Seq DEA.
#   # Tries to find the loci of the DE genes in the GISTIC focal events dataset. If it is found,
#   # runs a test to see if the CNV data also reflects a difference between groups. If it does,
#   # checks if the direction of copy number matches that of the RNA-Seq DE.
#   #
#   # WARNING: assumes ordinal levels are properly defined.
#   #
#   # @param <DEres> a dataset resulting from a DESeq2 DE analysis
#   # @param <genesMap> a biomart data to identify the locus on which the gene is, based on ENSEMBL ID
#   # @param <CNVdata> the GISTIC focal events dataset
#   # @param <CNVregions> a 3-columns ("chromosome", "start", "end") dataframe with rownames matching that
#   #                     of <CNVdata> to provide loci information.
#   # @param <responseVariable> the response variable for the test of the CNV DE. Must be a vector named in the same manner as <CNVdata>
#   # @return <CNVvalidationData> a dataframe indicating for each gene in <DEres> the p-value of the test
#   #                             on CNV data being DE as well (-1 if gene not in CNV data, -2 if gene is not recognized
#   #                             at all), and if that DE is in the same direction as in the RNA-Seq data
#   validateWithTCGACNV <- function(DEres, genesMap, CNVdata, responseVariable){
#     
#     # Extract ENSEMBL ID (remove the version '.xx' extension in rownames).
#     ensID = gsub("\\.\\d+", "", rownames(DEres))
#     
#     # Create an empty table which will be filled with the results.
#     CNVvalidationData <- data.frame(kruskal.pvalue = rep(NA, nrow(DEres)),
#                                     comment = rep(NA, nrow(DEres)))
#     rownames(CNVvalidationData) <- rownames(DEres)
#     
#     # Run the test for each gene in <DEres>.
#     for(g in 1:nrow(DEres)){ 
#       
#       ## Find gene in TCGA CNV data.
#       posInCNV <- match(DEres$geneName[g], CNVdata$Hugo_Symbol)
#       if(is.na(posInCNV)){ # if it could be found by name, check Entrez ID
#         posInCNV <- match(genesMap$entrezgene[match(ensID[g], genesMap$ensembl_gene_id)], CNVdata$Entrez_Gene_Id)
#       }
#       cat("Begin. g = ", g, " ; posInCNV = ", posInCNV, "\n")
#       
#       
#       if(!is.na(posInCNV)){
#         
#         # Create a dataset that combines...
#         dat <- data.frame(cnv = as.factor(CNVdata[posInCNV, 3:ncol(CNVdata)]), # the cnv data for the corresponding region
#                           response = responseVariable[match(colnames(CNVdata)[3:ncol(CNVdata)], names(responseVariable))]) # and the response
#         
#         # Test if CNV is dependent on <responseVariable>.
#         #print(as.numeric(dat$cnv))
#         #print(dat$response)
#         kTest <- kruskal.test(as.numeric(dat$cnv), dat$response)
#         CNVvalidationData[g, "kruskal.pvalue"] = kTest$p.value
#         
#         if(kTest$p.value < 0.06){
#           
#           # Test if direction of CNV is the same as DE.
#           spearCor = cor.test(as.numeric(dat$cnv), as.numeric(dat$response), method = "spearman")
#           if(spearCor$estimate * DEres[g, "log2FoldChange"] < 0){ # if the direction of differential expression does not match
#             CNVvalidationData[g, "comment"] = paste0("Opposite direction: rho=", round(spearCor$estimate, 3), 
#                                                      ", pval=", round(spearCor$p.value, 4))
#           }else{  
#             CNVvalidationData[g, "comment"] = paste0("Same direction: rho=", round(spearCor$estimate, 3), 
#                                                      ", pval=", round(spearCor$p.value, 4)) 
#           }
#         }
#       }
#     }# end of for loop
#     
#     return(CNVvalidationData)
#     
#   }# end of validateWithTCGACNV()
#   
#   # Define responses variables to use.
#   cliniDE_STLT <- clini_TCGA[clini_TCGA$SAMPLE_ID %in% colnames(cnv_tcga),]
#   cliniDE_STLT <- cliniDE_STLT[(cliniDE_STLT$OS_MONTHS>36 | cliniDE_STLT$OS_MONTHS<12),]
#   cliniDE_STLT$Survivor_Type <- ifelse(cliniDE_STLT$OS_MONTHS < 12, "ST", "LT")
#   cliniDE_STLT$Survivor_Type <- factor(cliniDE_STLT$Survivor_Type, levels = c("ST", "LT"), ordered = TRUE)
#   table(cliniDE_STLT$Survivor_Type, useNA = "ifany")
#   table(cliniDE_STLT$MGMT_STATUS, useNA = "ifany")
#   table(cliniDE_STLT$Survivor_Type, cliniDE_STLT$MGMT_STATUS, useNA = "ifany")
#   
#   SVLT_Type = factor(cliniDE_STLT$Survivor_Type, levels = c("ST", "LT"), ordered = TRUE)
#   SVLT_MGMT = factor(cliniDE_STLT$MGMT_STATUS, levels = c("UNMETHYLATED", "METHYLATED"), ordered = TRUE)
#   SVLT_SVMT = factor(paste0(cliniDE_STLT$Survivor_Type, cliniDE_STLT$MGMT_STATUS),
#                      levels = c("STUNMETHYLATED", "STMETHYLATED", "LTUNMETHYLATED", "LTMETHYLATED"), ordered = TRUE)
#   names(SVLT_Type) <- names(SVLT_MGMT) <- names(SVLT_SVMT) <- cliniDE_STLT$SAMPLE_ID
#   SVLT_MGMT <- na.omit(SVLT_MGMT)
#   SVLT_SVMT <- na.omit(SVLT_SVMT)
#   
#   
#   validate_TCGACNV_STvsLT = validateWithTCGACNV(DEres = resSignif_STvsLT, genesMap = genesMap, 
#                                                 CNVdata = cnv_tcga[, c(1,2,match(names(SVLT_Type), colnames(cnv_tcga)))],
#                                                 responseVariable = SVLT_Type)
#   
#   validate_TCGACNV_MGMT_Surv = validateWithTCGACNV(DEres = resSignif_MGMT_Surv, genesMap = genesMap, 
#                                                    CNVdata = cnv_tcga[, c(1,2,match(names(SVLT_Type), colnames(cnv_tcga)))],
#                                                    responseVariable = SVLT_Type)
#   validate_TCGACNV_MGMT_MGMT = validateWithTCGACNV(DEres = resSignif_MGMT_MGMT, genesMap = genesMap, 
#                                                    CNVdata = cnv_tcga[, c(1,2,match(names(SVLT_MGMT), colnames(cnv_tcga)))],
#                                                    responseVariable = SVLT_MGMT)
#   validate_TCGACNV_MGMT_MGMTSurv = validateWithTCGACNV(DEres = resSignif_MGMT_MGMTSurv, genesMap = genesMap, 
#                                                        CNVdata = cnv_tcga[, c(1,2,match(names(SVLT_SVMT), colnames(cnv_tcga)))],
#                                                        responseVariable = SVLT_SVMT)
#   
#   # validate_GTCNV_OS_Surv = validateWithCNV(DEres = resSignif_OS_Surv, genesMap = genesMap, 
#   #                                          CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#   #                                          CNVregions = cnvRegions, responseVariable = OS_all)
#   # validate_GTCNV_OS_MGMTSurv = validateWithCNV(DEres = resSignif_OS_MGMTSurv, genesMap = genesMap, 
#   #                                              CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#   #                                              CNVregions = cnvRegions, responseVariable = OS_all)
#   # validate_GTCNV_OS_MGMT = validateWithCNV(DEres = resSignif_OS_MGMT, genesMap = genesMap, 
#   #                                          CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#   #                                          CNVregions = cnvRegions, responseVariable = OS_all_MGMT)
#   
#   validate_GTCNV_sets = list(validate_GTCNV_STvsLT,
#                              validate_GTCNV_MGMT_Surv, validate_GTCNV_MGMT_MGMT, validate_GTCNV_MGMT_MGMTSurv,
#                              validate_GTCNV_OS_Surv, validate_GTCNV_OS_MGMTSurv, validate_GTCNV_OS_MGMT)
#   names(validate_GTCNV_sets) = c("validate_GTCNV_STvsLT",
#                                  "validate_GTCNV_MGMT_Surv", "validate_GTCNV_MGMT_MGMT", "validate_GTCNV_MGMT_MGMTSurv",
#                                  "validate_GTCNV_OS_Surv", "validate_GTCNV_OS_MGMTSurv", "validate_GTCNV_OS_MGMT")
#   lapply(validate_GTCNV_sets, function(d) (d[d$kruskal.pvalue != -1,]))
#   lapply(validate_GTCNV_sets, function(d) paste(sum(d$kruskal.pvalue != -1), "/", nrow(d)))
#   
#   
#   resSignif_MGMT_Surv["ENSG00000204666.3",]
#   gloup = wgs_data_cnv["chr19_48225001_54624999", match(names(SVLT_MGMT), colnames(wgs_data_cnv))]
#   cbind(cliniDE_STLT, as.numeric(gloup))[order(cliniDE_STLT$Survivor_Type),]
#   
# }# end of ### Compare to TCGA CNV.
# 
# 
# ### Compare to TCGA RNA-Seq 
# {
#   ## Load data
#   
#   
#   # Fetch ENSEMBL db to determine genomic source of the transcripts.
#   ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "89")
#   genesMap <- getBM(mart = ensembl,
#                     attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene", "external_gene_name", "chromosome_name", "start_position", "end_position"))
#   head(genesMap)
#   
#   validateWithTCGARNA <- function(DEres, genesMap, RNAdata, responseVariable){
#     
#     # Extract ENSEMBL ID (remove the version '.xx' extension in rownames).
#     ensID = gsub("\\.\\d+", "", rownames(DEres))
#     
#     # Create an empty table which will be filled with the results.
#     RNAvalidationData <- data.frame(MannWhitneyU.pvalue = rep(NA, nrow(DEres)),
#                                     comment = rep(NA, nrow(DEres)))
#     rownames(RNAvalidationData) <- rownames(DEres)
#     
#     # Run the test for each gene in <DEres>.
#     for(g in 1:nrow(DEres)){ 
#       
#       ## Find gene in normalized TCGA RNA-Seq data.
#       posInData <- match(DEres$geneName[g], RNAdata$Hugo_Symbol)
#       if(is.na(posInData)){ # if it could be found by name, check Entrez ID
#         posInData <- match(genesMap$entrezgene[match(ensID[g], genesMap$ensembl_gene_id)], RNAdata$Entrez_Gene_Id)
#       }
#       cat("Begin. g = ", g, " ; posInCNV = ", posInData, "\n")
#       
#       
#       if(!is.na(posInData)){
#         
#         # Create a dataset that combines...
#         dat <- data.frame(rna = as.numeric(RNAdata[posInData, 3:ncol(RNAdata) ]), # the cnv data for the corresponding region
#                           response = responseVariable[match(colnames(RNAdata)[3:ncol(RNAdata)], names(responseVariable))]) # and the response
#         
#         # Test if CNV is dependent on <responseVariable>.
#         #print(as.numeric(dat$cnv))
#         #print(dat$response)
#         # if(g>292){
#         #   print(dat$rna)
#         #   print(dat$response)
#         # }
#         if(sum(is.na(dat$rna)) == nrow(dat)){
#           RNAvalidationData[g, "comment"] = "No data"
#           wTestPval <- 1
#         }else{
#           wTest <- wilcox.test(dat$rna ~ dat$response)
#           wTestPval <- wTest$p.value
#           RNAvalidationData[g, "MannWhitneyU.pvalue"] = wTestPval
#         }
#         if(wTestPval < 0.05){
#           
#           # Test if direction of CNV is the same as DE.
#           spearCor = cor.test(as.numeric(dat$rna), as.numeric(dat$response), method = "spearman")
#           if(spearCor$estimate * DEres[g, "log2FoldChange"] < 0){ # if the direction of differential expression does not match
#             RNAvalidationData[g, "comment"] = paste0("Opposite direction: rho=", round(spearCor$estimate, 3), 
#                                                      ", pval=", round(spearCor$p.value, 4))
#           }else{  
#             RNAvalidationData[g, "comment"] = paste0("Same direction: rho=", round(spearCor$estimate, 3), 
#                                                      ", pval=", round(spearCor$p.value, 4)) 
#           }
#         }
#       }
#     }# end of for loop
#     
#     return(RNAvalidationData)
#     
#   }# end of validateWithTCGACNV()
#   
#   # Define responses variables to use.
#   cliniDE_STLT <- clini_TCGA[clini_TCGA$SAMPLE_ID %in% colnames(rna_zscores_tcga),]
#   cliniDE_STLT <- cliniDE_STLT[(cliniDE_STLT$OS_MONTHS>36 | cliniDE_STLT$OS_MONTHS<12),]
#   cliniDE_STLT$Survivor_Type <- ifelse(cliniDE_STLT$OS_MONTHS < 12, "ST", "LT")
#   cliniDE_STLT$Survivor_Type <- factor(cliniDE_STLT$Survivor_Type, levels = c("ST", "LT"), ordered = TRUE)
#   table(cliniDE_STLT$Survivor_Type, useNA = "ifany")
#   table(cliniDE_STLT$MGMT_STATUS, useNA = "ifany")
#   table(cliniDE_STLT$Survivor_Type, cliniDE_STLT$MGMT_STATUS, useNA = "ifany")
#   
#   SVLT_Type = factor(cliniDE_STLT$Survivor_Type, levels = c("ST", "LT"), ordered = TRUE)
#   SVLT_MGMT = factor(cliniDE_STLT$MGMT_STATUS, levels = c("UNMETHYLATED", "METHYLATED"), ordered = TRUE)
#   SVLT_SVMT = factor(paste0(cliniDE_STLT$Survivor_Type, cliniDE_STLT$MGMT_STATUS),
#                      levels = c("STUNMETHYLATED", "STMETHYLATED", "LTUNMETHYLATED", "LTMETHYLATED"), ordered = TRUE)
#   names(SVLT_Type) <- names(SVLT_MGMT) <- names(SVLT_SVMT) <- cliniDE_STLT$SAMPLE_ID
#   SVLT_MGMT <- na.omit(SVLT_MGMT)
#   SVLT_SVMT <- na.omit(SVLT_SVMT)
#   
#   
#   validate_TCGARNA_STvsLT = validateWithTCGARNA(DEres = resSignif_STvsLT, genesMap = genesMap, 
#                                                 RNAdata = rna_zscores_tcga[, c(1,2,match(names(SVLT_Type), colnames(rna_zscores_tcga)))],
#                                                 responseVariable = SVLT_Type)
#   
#   validate_TCGARNA_MGMT_Surv = validateWithTCGARNA(DEres = resSignif_MGMT_Surv, genesMap = genesMap, 
#                                                    RNAdata = rna_zscores_tcga[, c(1,2,match(names(SVLT_Type), colnames(rna_zscores_tcga)))],
#                                                    responseVariable = SVLT_Type)
#   validate_TCGARNA_MGMT_MGMT = validateWithTCGARNA(DEres = resSignif_MGMT_MGMT, genesMap = genesMap, 
#                                                    RNAdata = rna_zscores_tcga[, c(1,2,match(names(SVLT_MGMT), colnames(rna_zscores_tcga)))],
#                                                    responseVariable = SVLT_MGMT)
#   validate_TCGARNA_MGMT_MGMTSurv = validateWithTCGARNA(DEres = resSignif_MGMT_MGMTSurv, genesMap = genesMap, 
#                                                        RNAdata = rna_zscores_tcga[, c(1,2,match(names(SVLT_SVMT), colnames(rna_zscores_tcga)))],
#                                                        responseVariable = SVLT_SVMT)
#   
#   # validate_GTCNV_OS_Surv = validateWithCNV(DEres = resSignif_OS_Surv, genesMap = genesMap, 
#   #                                          CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#   #                                          CNVregions = cnvRegions, responseVariable = OS_all)
#   # validate_GTCNV_OS_MGMTSurv = validateWithCNV(DEres = resSignif_OS_MGMTSurv, genesMap = genesMap, 
#   #                                              CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#   #                                              CNVregions = cnvRegions, responseVariable = OS_all)
#   # validate_GTCNV_OS_MGMT = validateWithCNV(DEres = resSignif_OS_MGMT, genesMap = genesMap, 
#   #                                          CNVdata = wgs_data_cnv[, match(names(OS_all), colnames(wgs_data_cnv))],
#   #                                          CNVregions = cnvRegions, responseVariable = OS_all_MGMT)
#   
#   validate_GTCNV_sets = list(validate_GTCNV_STvsLT,
#                              validate_GTCNV_MGMT_Surv, validate_GTCNV_MGMT_MGMT, validate_GTCNV_MGMT_MGMTSurv,
#                              validate_GTCNV_OS_Surv, validate_GTCNV_OS_MGMTSurv, validate_GTCNV_OS_MGMT)
#   names(validate_GTCNV_sets) = c("validate_GTCNV_STvsLT",
#                                  "validate_GTCNV_MGMT_Surv", "validate_GTCNV_MGMT_MGMT", "validate_GTCNV_MGMT_MGMTSurv",
#                                  "validate_GTCNV_OS_Surv", "validate_GTCNV_OS_MGMTSurv", "validate_GTCNV_OS_MGMT")
#   lapply(validate_GTCNV_sets, function(d) (d[d$kruskal.pvalue != -1,]))
#   lapply(validate_GTCNV_sets, function(d) paste(sum(d$kruskal.pvalue != -1), "/", nrow(d)))
#   
#   
#   resSignif_MGMT_Surv["ENSG00000204666.3",]
#   gloup = wgs_data_cnv["chr19_48225001_54624999", match(names(SVLT_MGMT), colnames(wgs_data_cnv))]
#   cbind(cliniDE_STLT, as.numeric(gloup))[order(cliniDE_STLT$Survivor_Type),]
#   
# }# end of ### Compare to TCGA CNV.
# 
# # ===== II.4 Save results 
# save(dds, res, resLFC, resSignif, file = "./tmp/DESeq2_DE_STvsLT_20200707.RData")
# write.table(x = resSignif, file = "./results/2020.07.07_ST25/DESeq2_DEsignif.tsv", col.names = NA, sep = '\t')
# 
} # deprecated 2020.11.11


# ===== II.1 Run analysis ====
{
  
  ## Prepare objects
  coldata <- clinic_data[clinic_data$Subject_ID %in% colnames(rna_data_raw), c("Subject_ID", "Survivor_Type", "MGMT_Methylation_Status", "rna_batch", "Centers")]
  coldata_de = coldata[coldata$Survivor_Type != "IT", ]
  coldata_de$rna_batch <- as.factor(coldata_de$rna_batch)
  rna_data_de = rna_data_raw[, match(coldata_de$Subject_ID, colnames(rna_data_raw))]
  
  ddsMat_MGMT <- DESeqDataSetFromMatrix(countData = rna_data_de, colData = coldata_de, design = ~ rna_batch + Centers + MGMT_Methylation_Status + Survivor_Type + MGMT_Methylation_Status:Survivor_Type)
  dim(ddsMat_MGMT)
  ddsMat_MGMT <- ddsMat_MGMT[rowSums(counts(ddsMat_MGMT)) > 10, ] # filter out transcripts that have less than 10 count total
  dim(ddsMat_MGMT)
  
  ## Run DE analysis.
  dds_MGMT <- DESeq(ddsMat_MGMT)
  summary(dds_MGMT)
  resSurv <- results(dds_MGMT, contrast = c("Survivor_Type", "ST", "LT"), alpha = 0.05)
  resMGMT <- results(dds_MGMT, contrast = c("MGMT_Methylation_Status", "unmethylated", "methylated"), alpha = 0.05)
  resultsNames(dds_MGMT)
  resMGMTSurv <-results(dds_MGMT, name = "MGMT_Methylation_Statusunmethylated.Survivor_TypeST", alpha = 0.05)
  
  ## Filter significant results.
  hist(resSurv$padj, breaks = 100)
  sum(resSurv$padj < 0.05, na.rm = TRUE)
  summary(resSurv)
  resSignif_MGMT_Surv <- resSurv[which(resSurv$padj < 0.05), ]
  resSignif_MGMT_Surv <- resSignif_MGMT_Surv[order(resSignif_MGMT_Surv$padj),]
  
  hist(resMGMT$padj, breaks = 100)
  sum(resMGMT$padj < 0.05, na.rm = TRUE)
  summary(resMGMT)
  resSignif_MGMT_MGMT <- resMGMT[which(resMGMT$padj < 0.05), ]
  resSignif_MGMT_MGMT <- resSignif_MGMT_MGMT[order(resSignif_MGMT_MGMT$padj),]
  
  hist(resMGMTSurv$padj, breaks = 100)
  sum(resMGMTSurv$padj < 0.05, na.rm = TRUE)
  summary(resMGMTSurv)
  resSignif_MGMT_MGMTSurv <- resMGMTSurv[which(resMGMTSurv$padj < 0.05), ]
  resSignif_MGMT_MGMTSurv <- resSignif_MGMT_MGMTSurv[order(resSignif_MGMT_MGMTSurv$padj),]
  
  resSurv$geneName <- rna2genes[match(rownames(resSurv), rna2genes$Transcript_ID), 2]
  resSignif_MGMT_Surv$geneName <- rna2genes[match(rownames(resSignif_MGMT_Surv), rna2genes$Transcript_ID), 2]
  resSignif_MGMT_MGMT$geneName <- rna2genes[match(rownames(resSignif_MGMT_MGMT), rna2genes$Transcript_ID), 2]
  resSignif_MGMT_MGMTSurv$geneName <- rna2genes[match(rownames(resSignif_MGMT_MGMTSurv), rna2genes$Transcript_ID), 2]
  
  ## Visualize results.
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_Surv)[1], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_Surv)[2], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_Surv)[3], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMT)[1], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMT)[2], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMT)[3], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMTSurv)[1], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMTSurv)[2], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  plotCounts(dds_MGMT, rownames(resSignif_MGMT_MGMTSurv)[3], intgroup = c("Survivor_Type", "MGMT_Methylation_Status"),
             normalized = TRUE, transform = TRUE, pch=20)
  
  plotMA(resSignif_MGMT_Surv, ylim=c(-2,2))
  plotMA(resSignif_MGMT_MGMT, ylim=c(-2,2))
  plotMA(resSignif_MGMT_MGMTSurv, ylim=c(-2,2))
  
  ## Normalize log2foldchange to [-1:1]
  exp = resSurv$log2FoldChange
  n_exp = 2*((exp-min(exp))/(max(exp)-min(exp))) - 1
  hist(exp, breaks = 100)
  hist(n_exp, breaks = 100, ylim=c(0,100))
  
}



# ===== II.2 Explore results =====

### Perform Enrichment Analysis 
{

  geneUniverse = data.frame(probe.name = rna2genes$Gene_Symbol, probe.symbol = rna2genes$Gene_Symbol) # needed for runTopGO()

  runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_Surv$geneName, pvalues = resSignif_MGMT_Surv$padj,
           outPath = "./results/2020.11.12/GT.RNA-Seq_DEA/", title="TopGo__MGMT_Surv_count", 
           ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
  runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_Surv$geneName, pvalues = resSignif_MGMT_Surv$padj,
           outPath = "./results/2020.11.12/GT.RNA-Seq_DEA/", title="TopGo__MGMT_Surv_score", 
           ontology="BP", algorithm = "elim", statistic = "ks", topNods = 20)
  # runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_MGMT$geneName, pvalues = resSignif_MGMT_MGMT$padj,
  #          outPath = "./results/2020.11.12/GT.RNA-Seq_DEA/", title="TopGo_analysis_DESeq2_MGMT_MGMT", 
  #          ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
  # runTopGO(geneUniverse = geneUniverse, investigatedGenes = resSignif_MGMT_MGMTSurv$geneName, pvalues = resSignif_MGMT_MGMTSurv$padj,
  #          outPath = "./results/2020.11.12/GT.RNA-Seq_DEA/", title="TopGo_analysis_DESeq2_MGMT_MGMTSurv", 
  #          ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)

}



# ===== II.3 Compare to other datasets ====

names(resSignif_MGMT_Surv)[2] = "difference"

## With GT data.
{
  clini_validation_STLT = clinic_data[clinic_data$Survivor_Type != "IT",]
  clini_validation_IT = clinic_data[clinic_data$Survivor_Type == "IT",]
  
  # Break down the rownames of cnv data to get loci information.
  cnvRegions <- as.data.frame(t(data.frame(sapply(rownames(wgs_data_cnv), function(reg) strsplit(reg, '_')))))
  colnames(cnvRegions) <- c("chromosome", "start", "end")
  cnvRegions$start <- as.numeric(as.character(cnvRegions$start))
  cnvRegions$end <- as.numeric(as.character(cnvRegions$end))
  
  # CNV data.
  samples_STLT = intersect(clini_validation_STLT$Subject_ID, colnames(wgs_data_cnv))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, clini_validation_STLT$Subject_ID)]
  names(responseVariable) = clini_validation_STLT$Subject_ID[match(samples_STLT, clini_validation_STLT$Subject_ID)]
  validation_GT.CNV.STLT = validate_with_GT.CNV.STLT(DEres = resSignif_MGMT_Surv, genesMap = genesMap, CNVregions = cnvRegions,
                                                     validationData = wgs_data_cnv[, match(samples_STLT, colnames(wgs_data_cnv))],
                                                     responseVariable = responseVariable)
  
  samples_IT = intersect(clini_validation_IT$Subject_ID, colnames(wgs_data_cnv))
  responseVariable = clini_validation_IT$Overall_Survival_Overall_Survival_months[match(samples_IT, clini_validation_IT$Subject_ID)]
  names(responseVariable) = clini_validation_IT$Subject_ID[match(samples_IT, clini_validation_IT$Subject_ID)]
  validation_GT.CNV.IT = validate_with_GT.CNV.IT(DEres = resSignif_MGMT_Surv, genesMap = genesMap, CNVregions = cnvRegions,
                                                 validationData = wgs_data_cnv[, match(samples_IT, colnames(wgs_data_cnv))],
                                                 responseVariable = responseVariable)
  
  
  # RNA-Seq data.
  samples_IT = intersect(clini_validation_IT$Subject_ID, colnames(rna_data_vst))
  responseVariable = clini_validation_IT$Overall_Survival_Overall_Survival_months[match(samples_IT, clini_validation_IT$Subject_ID)]
  names(responseVariable) = clini_validation_IT$Subject_ID[match(samples_IT, clini_validation_IT$Subject_ID)]
  validation_GT.RNASeq.IT = validate_with_GT.RNASeq.IT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                       validationData = rna_data_vst[, match(samples_IT, colnames(rna_data_vst))],
                                                       responseVariable = responseVariable)
  
  
  # See results.
  sum(!is.na(validation_GT.CNV.STLT$matching.region))
  sum(!is.na(validation_GT.CNV.IT$matching.region))
  sapply(list(validation_GT.CNV.STLT, validation_GT.CNV.IT, validation_GT.RNASeq.IT), 
         function(d) nrow(d) - sum(d$test.pvalue %in% c(-1,-2)))
  sapply(list(validation_GT.CNV.STLT, validation_GT.CNV.IT, validation_GT.RNASeq.IT), 
         function(d) sum(!is.na(d$comment)))
  sapply(list(validation_GT.CNV.STLT, validation_GT.CNV.IT, validation_GT.RNASeq.IT), 
         function(d) d[grep("Same",d$comment),])
  sapply(list(validation_GT.CNV.STLT, validation_GT.CNV.IT, validation_GT.RNASeq.IT), 
         function(d) d[grep("Opposite",d$comment),])
  genesMap[genesMap$ensembl_gene_id %in% c("ENSG00000269956", "ENSG00000130957"), "hgnc_symbol"]
  genesMap[genesMap$ensembl_gene_id %in% c("ENSG00000131044", "ENSG00000168658", "ENSG00000215612", "ENSG00000175920", "ENSG00000184471"), "hgnc_symbol"]
  
  
}

## With DASL data.
{
  clini_validation_STLT = clini_dasl[clini_dasl$Survivor_Type != "IT",]
  clini_validation_IT = clini_dasl[clini_dasl$Survivor_Type == "IT",]
  
  samples_STLT = intersect(rownames(clini_validation_STLT), colnames(dasl_data))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, rownames(clini_validation_STLT))]
  names(responseVariable) = rownames(clini_validation_STLT[match(samples_STLT, rownames(clini_validation_STLT)),])
  validation_EMC.DASL.STLT = validate_with_DASL.STLT(DEres = resSignif_MGMT_Surv, genesMap = genesMap, 
                                                     validationData = dasl_data[, match(samples_STLT, colnames(dasl_data))],
                                                     responseVariable = responseVariable)
  
  samples_IT = intersect(rownames(clini_validation_IT), colnames(dasl_data))
  responseVariable = clini_validation_IT$OS[match(samples_IT, rownames(clini_validation_IT))]
  names(responseVariable) = rownames(clini_validation_IT[match(samples_IT, rownames(clini_validation_IT)),])
  validation_EMC.DASL.IT = validate_with_DASL.IT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                 validationData = dasl_data[, match(samples_IT, colnames(dasl_data))],
                                                 responseVariable = responseVariable)

  
  # See results.
  sapply(list(validation_EMC.DASL.STLT, validation_EMC.DASL.IT), 
         function(d) nrow(d) - sum(d$test.pvalue %in% c(-1,-2))) # number of features found in the validation dataset
  sapply(list(validation_EMC.DASL.STLT, validation_EMC.DASL.IT), 
         function(d) sum(!is.na(d$comment))) # number of features significantly DE as well in the validation dataset
  lapply(list(validation_EMC.DASL.STLT, validation_EMC.DASL.IT), 
         function(d) d[grep("Same",d$comment),])
  lapply(list(validation_EMC.DASL.STLT, validation_EMC.DASL.IT), 
         function(d) d[grep("Opposite",d$comment),])
  genesMap[genesMap$ensembl_gene_id %in% c("ENSG00000269956", "ENSG00000130957"), "hgnc_symbol"]
  genesMap[genesMap$ensembl_gene_id %in% c("ENSG00000131044", "ENSG00000168658", "ENSG00000215612", "ENSG00000175920", "ENSG00000184471"), "hgnc_symbol"]
  
  
}

## With TCGA data.
{
  clini_validation_STLT = clini_tcga[clini_tcga$Survivor_Type != "IT",]
  clini_validation_IT = clini_tcga[clini_tcga$Survivor_Type == "IT",]
  
  # CNV data.
  samples_STLT = intersect(clini_validation_STLT$SAMPLE_ID, colnames(cnv_tcga))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_STLT$SAMPLE_ID[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  validation_TCGA.CNV.STLT = validate_with_TCGA.CNV.STLT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                     validationData = cnv_tcga[, c(1,2,3, match(samples_STLT, colnames(cnv_tcga)))],
                                                     responseVariable = responseVariable)
  
  samples_IT = intersect(clini_validation_IT$SAMPLE_ID, colnames(cnv_tcga))
  responseVariable = clini_validation_IT$OS_MONTHS[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_IT$SAMPLE_ID[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  validation_TCGA.CNV.IT = validate_with_TCGA.CNV.IT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                 validationData = cnv_tcga[, c(1,2,3, match(samples_IT, colnames(cnv_tcga)))],
                                                 responseVariable = responseVariable)
  
  # RNA-Seq data.
  samples_STLT = intersect(clini_validation_STLT$SAMPLE_ID, colnames(rna_zscores_tcga))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_STLT$SAMPLE_ID[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  validation_TCGA.RNASeq.STLT = validate_with_TCGA.RNASeq.STLT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                           validationData = rna_zscores_tcga[, c(1,2, match(samples_STLT, colnames(rna_zscores_tcga)))],
                                                           responseVariable = responseVariable)
  
  samples_IT = intersect(clini_validation_IT$SAMPLE_ID, colnames(rna_zscores_tcga))
  responseVariable = clini_validation_IT$OS_MONTHS[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_IT$SAMPLE_ID[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  validation_TCGA.RNASeq.IT = validate_with_TCGA.RNASeq.IT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                       validationData = rna_zscores_tcga[, c(1,2, match(samples_IT, colnames(rna_zscores_tcga)))],
                                                       responseVariable = responseVariable)
  
  
  # See results.
  
}



# ===== II.4 Save results =====

save(ddsMat_MGMT, dds_MGMT, resSurv, resMGMT, resSignif_MGMT_Surv, resSignif_MGMT_MGMT,
     file = "./results/2020.11.12/GT.RNA-Seq_DEA/DEAresults.RData")
write.table(x = resSignif_MGMT_Surv, col.names = NA, sep = '\t', 
            file = "./results/2020.11.12/GT.RNA-Seq_DEA/DESeq2_DEsignif_survGenes.tsv")
save(validation_GT.CNV.STLT, validation_GT.CNV.IT, validation_GT.RNASeq.IT, validation_EMC.DASL.STLT, validation_EMC.DASL.IT,
     validation_TCGA.CNV.STLT, validation_TCGA.CNV.IT, validation_TCGA.RNASeq.STLT, validation_TCGA.RNASeq.IT,
     file = "./results/2020.11.12/GT.RNA-Seq_DEA/crossvalidation_results.RData")




# oooooooooooooooooooooooooooooooooooooooooooooooooo
# III - DE Analysis for WGS GISTIC CNV data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

# ===== III.1 Run analysis =====

# ST vs LT accounting for MGMT ----

### Prepare data.
clini_analysis_STLT = clinic_data[clinic_data$Survivor_Type != "IT",]
samples_STLT = intersect(clini_analysis_STLT$Subject_ID, colnames(wgs_data_cnv))
clini_analysis = clini_analysis_STLT[match(samples_STLT, clini_analysis_STLT$Subject_ID), c("Subject_ID", "Survivor_Type", "MGMT_Methylation_Status")]
rownames(clini_analysis) = clini_analysis$Subject_ID
clini_analysis = clini_analysis[, c("Survivor_Type", "MGMT_Methylation_Status")]
wgs_data_analysis = wgs_data_cnv[, match(samples_STLT, colnames(wgs_data_cnv))]
dea_gt.cnv_results = data.frame(survival.coef = rep(NA, nrow(wgs_data_analysis)),
                                mgmt.coef = rep(NA, nrow(wgs_data_analysis)) ) # will be filled as each region is analysed
rownames(dea_gt.cnv_results) = rownames(wgs_data_analysis)


### Analyse for each gene.
for(region in 1:nrow(wgs_data_analysis)){
  
  ## Define analysis matrix.
  dat = data.frame(Survivor_Type = clini_analysis$Survivor_Type, 
                   MGMT_Methylation_Status = clini_analysis$MGMT_Methylation_Status,
                   Region_cnv = as.numeric(wgs_data_analysis[region,]))
  dat$Survivor_Type = factor(dat$Survivor_Type, levels = c("ST", "LT"), ordered = TRUE)
  dat$MGMT_Methylation_Status = factor(dat$MGMT_Methylation_Status, levels = c("unmethylated", "methylated"), ordered = TRUE)
  dat$Region_cnv = factor(dat$Region_cnv, levels = c(-1, 0, 1, 2), ordered = TRUE)
  
  ## Train the model.
  model = polr(Region_cnv ~ Survivor_Type + MGMT_Methylation_Status, data = dat, Hess = TRUE)
  model.summary = summary(model)
  
  ## Complete the results matrix.
  dea_gt.cnv_results$survival.coef[region] = model.summary$coefficients[1,1]
  dea_gt.cnv_results$mgmt.coef[region] = model.summary$coefficients[2,1]
  
}


### Check out results.
hist(dea_gt.cnv_results$survival.coef, breaks = 100)
hist(dea_gt.cnv_results$mgmt.coef, breaks = 100)



# ===== III.2 Explore results =====

# Curate results ----
{
  ### Check out distributions.
  hist(dea_gt.cnv_results$survival.coef, breaks = 100)
  hist(dea_gt.cnv_results$mgmt.coef, breaks = 100)
  sum(abs(dea_gt.cnv_results$survival.coef) >=1 )
  sum(abs(dea_gt.cnv_results$mgmt.coef) >=1 )
  
  ### Extract genes from shortlisted regions.
  dea_gt.cnv_results.survival = data.frame(geneName = NA, difference = NA, region = NA)
  dea_gt.cnv_results.mgmt = data.frame(geneName = NA, difference = NA, region = NA)
  
  for(reg in 1:nrow(dea_gt.cnv_results)){
    
    if(abs(dea_gt.cnv_results$survival.coef[reg]) >=1){
      region = strsplit(rownames(dea_gt.cnv_results)[reg], '_')[[1]]
      genesList = genesMap$ensembl_gene_id[genesMap$chromosome_name == gsub("chr", "", region[1]) & 
                                             genesMap$start_position >= region[2] & genesMap$end_position >= region[3] ]
      toAdd = data.frame(geneName = genesList, difference = rep(dea_gt.cnv_results$survival.coef[reg], length(genesList)), 
                         region = rep(rownames(dea_gt.cnv_results)[reg], length(genesList)))
      dea_gt.cnv_results.survival = rbind(dea_gt.cnv_results.survival, toAdd)
    }
    
    if(abs(dea_gt.cnv_results$mgmt.coef[reg]) >=1){
      region = strsplit(rownames(dea_gt.cnv_results)[reg], '_')[[1]]
      genesList = genesMap$ensembl_gene_id[genesMap$chromosome_name == gsub("chr", "", region[1]) & 
                                             genesMap$start_position >= region[2] & genesMap$end_position >= region[3] ]
      toAdd = data.frame(geneName = genesList, difference = rep(dea_gt.cnv_results$mgmt.coef[reg], length(genesList)), 
                         region = rep(rownames(dea_gt.cnv_results)[reg], length(genesList)))
      dea_gt.cnv_results.mgmt = rbind(dea_gt.cnv_results.mgmt, toAdd)
    }
  }
  
  dea_gt.cnv_results.survival = dea_gt.cnv_results.survival[-1,] 
  dea_gt.cnv_results.mgmt = dea_gt.cnv_results.mgmt[-1,] 
  
}

# Perform Enrichment Analysis ----
{
  ### Define all genes that were included in the analysis.
  ## Identify genes included in the analysis.
  geneUniverse_probes = c()
  
  for(reg in 1:nrow(dea_gt.cnv_results)){
      region = strsplit(rownames(dea_gt.cnv_results)[reg], '_')[[1]]
      genesList = which(genesMap$chromosome_name == gsub("chr", "", region[1]) & 
                          genesMap$start_position >= region[2] & genesMap$end_position >= region[3] )
      geneUniverse_probes = c(geneUniverse_probes, genesMap$ensembl_gene_id[genesList])
  }
  geneUniverse_probes = unique(geneUniverse_probes)
  
  ## Find out corresponding gene names.
  geneUniverse_symbols = c()
  for(g in geneUniverse_probes){
    posInMart = which(genesMap$ensembl_gene_id == g)
    if(length(posInMart) > 1){  posInMart = min(posInMart)  }
    geneUniverse_symbols = c(geneUniverse_symbols, 
                             ifelse(genesMap$hgnc_symbol[posInMart] != "", 
                                    genesMap$hgnc_symbol[posInMart], 
                                    genesMap$external_gene_name[posInMart]))
  }
  
  geneUniverse = data.frame(probe.name = geneUniverse_probes, probe.symbol = geneUniverse_symbols)
  
  
  ### Run gene enrichment analysis.
  runTopGO(geneUniverse = geneUniverse, investigatedGenes = dea_gt.cnv_results.survival$geneName,
           outPath = "./results/2020.11.12/GT.CNV_DEA/", title="TopGo_analysis_DESeq2_MGMT_Surv", 
           ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)
  runTopGO(geneUniverse = geneUniverse, investigatedGenes = dea_gt.cnv_results.mgmt$geneName,
           outPath = "./results/2020.11.12/GT.CNV_DEA/", title="TopGo_analysis_DESeq2_MGMT_MGMT", 
           ontology="BP", algorithm = "elim", statistic = "fisher", topNods = 20)

}


# ===== III.3 Compare to other datasets =====

## With GT data.
{
  clini_validation_STLT = clinic_data[clinic_data$Survivor_Type != "IT",]
  clini_validation_IT = clinic_data[clinic_data$Survivor_Type == "IT",]
  
  # Break down the rownames of cnv data to get loci information.
  cnvRegions <- as.data.frame(t(data.frame(sapply(rownames(wgs_data_cnv), function(reg) strsplit(reg, '_')))))
  colnames(cnvRegions) <- c("chromosome", "start", "end")
  cnvRegions$start <- as.numeric(as.character(cnvRegions$start))
  cnvRegions$end <- as.numeric(as.character(cnvRegions$end))
  
  # CNV data.
  samples_IT = intersect(clini_validation_IT$Subject_ID, colnames(wgs_data_cnv))
  responseVariable = clini_validation_IT$Overall_Survival_Overall_Survival_months[match(samples_IT, clini_validation_IT$Subject_ID)]
  names(responseVariable) = clini_validation_IT$Subject_ID[match(samples_IT, clini_validation_IT$Subject_ID)]
  validation_GT.CNV.IT = validate_with_GT.CNV.IT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap, CNVregions = cnvRegions,
                                                 validationData = wgs_data_cnv[, match(samples_IT, colnames(wgs_data_cnv))],
                                                 responseVariable = responseVariable)
  sum(!is.na(validation_GT.CNV.IT$matching.region))
  
  # RNA-Seq data.
  samples_STLT = intersect(clini_validation_STLT$Subject_ID, colnames(rna_data_vst))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, clini_validation_STLT$Subject_ID)]
  names(responseVariable) = clini_validation_STLT$Subject_ID[match(samples_STLT, clini_validation_STLT$Subject_ID)]
  validation_GT.RNASeq.STLT = validate_with_GT.RNASeq.STLT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap,
                                                       validationData = rna_data_vst[, match(samples_STLT, colnames(rna_data_vst))],
                                                       responseVariable = responseVariable)
  
  samples_IT = intersect(clini_validation_IT$Subject_ID, colnames(rna_data_vst))
  responseVariable = clini_validation_IT$Overall_Survival_Overall_Survival_months[match(samples_IT, clini_validation_IT$Subject_ID)]
  names(responseVariable) = clini_validation_IT$Subject_ID[match(samples_IT, clini_validation_IT$Subject_ID)]
  validation_GT.RNASeq.IT = validate_with_GT.RNASeq.IT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap,
                                                       validationData = rna_data_vst[, match(samples_IT, colnames(rna_data_vst))],
                                                       responseVariable = responseVariable)
  
  
  # See results.
  
  
}

## With DASL data.
{
  clini_validation_STLT = clini_dasl[clini_dasl$Survivor_Type != "IT",]
  clini_validation_IT = clini_dasl[clini_dasl$Survivor_Type == "IT",]
  
  samples_STLT = intersect(rownames(clini_validation_STLT), colnames(dasl_data))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, rownames(clini_validation_STLT))]
  names(responseVariable) = rownames(clini_validation_STLT[match(samples_STLT, rownames(clini_validation_STLT)),])
  validation_EMC.DASL.STLT = validate_with_DASL.STLT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap, 
                                                     validationData = dasl_data[, match(samples_STLT, colnames(dasl_data))],
                                                     responseVariable = responseVariable)
  
  samples_IT = intersect(rownames(clini_validation_IT), colnames(dasl_data))
  responseVariable = clini_validation_IT$OS[match(samples_IT, rownames(clini_validation_IT))]
  names(responseVariable) = rownames(clini_validation_IT[match(samples_IT, rownames(clini_validation_IT)),])
  validation_EMC.DASL.IT = validate_with_DASL.IT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap,
                                                 validationData = dasl_data[, match(samples_IT, colnames(dasl_data))],
                                                 responseVariable = responseVariable)
  
  
  # See results.
  
}

## With TCGA data.
{
  clini_validation_STLT = clini_tcga[clini_tcga$Survivor_Type != "IT",]
  clini_validation_IT = clini_tcga[clini_tcga$Survivor_Type == "IT",]
  
  # CNV data.
  samples_STLT = intersect(clini_validation_STLT$SAMPLE_ID, colnames(cnv_tcga))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_STLT$SAMPLE_ID[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  validation_TCGA.CNV.STLT = validate_with_TCGA.CNV.STLT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap,
                                                         validationData = cnv_tcga[, c(1,2,3, match(samples_STLT, colnames(cnv_tcga)))],
                                                         responseVariable = responseVariable)
  
  samples_IT = intersect(clini_validation_IT$SAMPLE_ID, colnames(cnv_tcga))
  responseVariable = clini_validation_IT$OS_MONTHS[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_IT$SAMPLE_ID[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  validation_TCGA.CNV.IT = validate_with_TCGA.CNV.IT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap,
                                                     validationData = cnv_tcga[, c(1,2,3, match(samples_IT, colnames(cnv_tcga)))],
                                                     responseVariable = responseVariable)
  
  # RNA-Seq data.
  samples_STLT = intersect(clini_validation_STLT$SAMPLE_ID, colnames(rna_zscores_tcga))
  responseVariable = clini_validation_STLT$Survivor_Type[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_STLT$SAMPLE_ID[match(samples_STLT, clini_validation_STLT$SAMPLE_ID)]
  validation_TCGA.RNASeq.STLT = validate_with_TCGA.RNASeq.STLT(DEres = dea_gt.cnv_results.survival, genesMap = genesMap,
                                                               validationData = rna_zscores_tcga[, c(1,2, match(samples_STLT, colnames(rna_zscores_tcga)))],
                                                               responseVariable = responseVariable)
  
  samples_IT = intersect(clini_validation_IT$SAMPLE_ID, colnames(rna_zscores_tcga))
  responseVariable = clini_validation_IT$OS_MONTHS[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  names(responseVariable) = clini_validation_IT$SAMPLE_ID[match(samples_IT, clini_validation_IT$SAMPLE_ID)]
  validation_TCGA.RNASeq.IT = validate_with_TCGA.RNASeq.IT(DEres = resSignif_MGMT_Surv, genesMap = genesMap,
                                                           validationData = rna_zscores_tcga[, c(1,2, match(samples_IT, colnames(rna_zscores_tcga)))],
                                                           responseVariable = responseVariable)
  
  
  # See results.
  
}


# ===== III.4 Save results =====










