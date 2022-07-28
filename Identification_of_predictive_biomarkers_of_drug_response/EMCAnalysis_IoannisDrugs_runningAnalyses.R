# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#                   Erasmus Data Analysis
#
# Run analysis of the EMC data, using PCA, WGCNA, LASSO and other statistiscal methods. 
#
# Investigated data includes:
#   * DASL transcriptomic microarray
#   * IC50 values for 108 drugs in 45 GBM samples
#   * Clinical data: diagnosis and grade of the tumour and patient's overall survival
#
# Several subsets of the data are investigated:
#   * The complete set of 88 samples from the DASL data (keyword: "all88")
#   * All GBM-diagnosed samples from the complete set of samples (keyword: "allGBM")
#   * All Primary GBM-diagnosed samples from the complete set of samples (keyword: "allPrimaryGBM")
#   * All Recurrent GBM-diagnosed samples from the complete set of samples (keyword: "allRecurrentGBM")
#   * All 45 samples that are present in the IC50 dataset (keyword: "drugExposedAllRecurrent")
#   * All Primary samples that are present in the IC50 dataset (keyword: "drugExposedPrimary")
#   * All Recurrent samples that are present in the IC50 dataset (keyword: "drugExposedRecurrent")
#
# The analysis mostly consists in trying to identify markers in the DASL data 
#    for either the IC50 drug response (keyword: "IC50s"), 
#    or the recurrence status + the overall survival + grade (keyword: "Clinical").
# Several approaches are taken:
#   * An unbiased approach (keyword: "unbiased")
#   * Focusing on known cancer pathways/genes (keyword: "cancerGenes")
#   * Focusing on the ALDH genes and Retinoic Acid Pathway (keyword: "aldhGenes")
#
# Started on the 12th of March, 2020. Copied from previous version
# 
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

# oooooooooooooooooooooooooooooooooooooooooooooooooo
# 0 - Get ready ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # ===== 0.1 Set-up environment =====
  setwd("/path/to/working_directory/")
  
  library(ade4)
  library(WGCNA)
  library(glmnet)
  library(topGO)
  library(Rgraphviz)
  
  # ===== 0.2 Define functions =====
  
  runLASSO <- function(x, y, runNbr = 100, freqRatio = 0.5){
    # x is the independent variables dataset
    # y is the dependent variables dataset
    # runNbr is the number of times LASSO is run
    # freqRatio is the % of times a given gene must be associated with a given trait to be included in final results

    runs = list() # results of the simulations
    y.df = as.data.frame(y)
    
    ## Start the runs
    for(run in 1:runNbr){
      
      cvfits = list()
      coefs_1se = list()
      coefs_min = list()
      
      for(trait in 1:ncol(y.df)){
        
        # Select samples for which we have a value
        notNA = !is.na(y.df[ , trait])
        clini_trait = y.df[notNA, trait]
        transc_trait = as.matrix(x[notNA, ])
        
        # Run glmnet
        cvfit = cv.glmnet(transc_trait, clini_trait, family="gaussian", alpha=1)
        cvfits[[trait]] = cvfit
        coefs_1se[[trait]] = coef(cvfit, s="lambda.1se") 
        coefs_min[[trait]] = coef(cvfit, s="lambda.min") 
        
      }
      
      names(cvfits) = colnames(y.df)
      names(coefs_1se) = colnames(y.df)
      names(coefs_min) = colnames(y.df)
      
      co1se = sapply(coefs_1se, function(d)  sum(as.matrix(d)>0)  )
      co1se = co1se[co1se>1]
      coMin = sapply(coefs_min, function(d)  sum(as.matrix(d)>0)  )
      traitMin = names(coMin)
      
      runResults = list()
      for(d in 1:length(traitMin)){
        coefs = as.matrix(coefs_min[[which(names(coefs_min)==traitMin[d])]])
        non0coefs = coefs[coefs>0,]
        
        runResults[[d]] = names(non0coefs)
        #if(length(runResults[[d]])>0) names(runResults[[d]]) = traitMin[d]
        
      }
      
      if(length(runResults)>0){  names(runResults) = traitMin[1:length(runResults)]  }
      #for(i in 1:length(runResults)){  if(length(runResults[[i]])>0) names(runResults[[i]]) = traitMin[i]  }
      #names(runResults) = traitMin
      
      runs[[run]] = runResults
      
    }
    
    
    ## Extract variables regularly associated to each trait by LASSO.
    LASSO_results = list()
    for(trait in 1:ncol(y.df)){
      traitName = colnames(y.df)[trait]
      traitGenes = c()
      traitFreq = 0
      for(run in 1:length(runs)){
        if(traitName %in% names(runs[[run]])){
          if(length(runs[[run]][[traitName]])>0){
            traitGenes = c(traitGenes, runs[[run]][[traitName]])
            traitFreq = traitFreq + 1
          }
        }
      }
      genesRatio = (table(traitGenes))/traitFreq
      
      results_probes = names(genesRatio[genesRatio >= freqRatio])
      LASSO_results[[trait]] = as.character(results_probes)#(probesGenes[match(results_probes, probesGenes$ProbeID),1])
      
    }
    names(LASSO_results) = colnames(y.df)
    
    ## Put results in a matrix
    LassoResults = matrix(data=NA, nrow=length(LASSO_results), ncol = max(sapply(LASSO_results, length)))
    for(trait in 1:length(LASSO_results)){
      for(var in 1:length(LASSO_results[[trait]])){
        LassoResults[trait, var] = LASSO_results[[trait]][var]
      }
    }
    rownames(LassoResults) = names(LASSO_results)
    
    return(LassoResults)
    
  } # end of runLASSO()
  
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

  runWGCNA <- function(path="./", title="WCNArun", omicData, samplesTraits, flag=1, maxBlockSize=NULL){
    {
      ## Performs Weigthed Gene Correlation Analysis to compare molecular data to other traits such as clinical data.
      # @params <omicData> the quantitative dataset to be analysed, with samples as rows and genes as colums                    
      # @params <samplesTraits> the external traits data, such as clinical features, to be correlated with identified modules;
      #                         Values in this dataset must also be numeric.
      # @params <flag> A value used to determine where to start the analysis; this is helpful if the analysis has already been run partly
      #                 or there is a need to try different parameters. Values correspond to:
      #                    - 1 : start from the beginning
      #                    - 2 : preprocessing and outliers exclusion has been performed. start at the soft thresholding power selection step
      #                    - 3 : soft thresholding power selection step is done. start at the Topology Overlap Matrix calculation step
      #                    - 4 : Topology Overlap Matrix calculation step is done. start at the Genes Clusters Analysis step
      #                    - 5 : Genes Clusters Analysis was performed; compute correlations with traits
      #                    - 6 : everything is done
      # @params <maxBlockSize> maximum size of a block to be analysed at a time, to be chosen depending on the size of 
      #                        the <omicData> dataset and the computer capacities; recommanded value for modest laptops is 5000
      #                        If left as NULL, the whole dataset will be analysed in one go
      # @params <title> name of the run; will be added to all outputed data and files for identification purposes 
      #
      # @return <mergedNetwork> the network resulting from merging elements form <dataList> using the SNF method
    } # runWGCNA header
    
    cat("=====  Launching WGCNA method ===== \n")
    outPath = paste0(path,title)
    
    while(flag < 6){
      if(flag==1){
        
        cat("==  First look at Datasets  == \n")
        gsg = goodSamplesGenes(omicData, verbose = 3)
        print(gsg$allOK)
        
        ## Cluster samples, detect and remove outliers 
        sampleTree = hclust(dist(omicData), method = "average");
        # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
        sizeGrWindow(12,9)
        par(cex = 0.6);
        par(mar = c(0,4,2,0))
        
        sampleCut = FALSE # used to loop until cutheight satisfies user
        while(sampleCut == FALSE){
          plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
          cutHeightSamples = readInteger("Enter cutHeigt:")
          abline(h = cutHeightSamples, col = "red"); # Plot a line to show the cut
          sampleCut = askConfirmation()
        }
        
        clust = cutreeStatic(sampleTree, cutHeight = cutHeightSamples, minSize = 3)
        cat("Samples clusters:\n")
        print(table(clust)) # Determine cluster under the line
        keepSamples = (clust!=0) # clust 1 contains the samples we want to keep.
        omic = omicData[keepSamples, ]
        nGenes = ncol(omic)
        nSamples = nrow(omic)
        # Save the plot to file
        dev.copy2pdf(file = paste0(outPath,"/graphs/samplesClustering.pdf"));
        dev.off()
        
        collectGarbage()
        
        save(omic, nGenes, nSamples, file = paste0(outPath,"/tmpData/flag1_dataWithoutOutliers.RData"))
        save(cutHeightSamples, file = paste0(outPath,"/tmpData/parameters.RData"))
        flag = 2
        
      }else if (flag==2){
        
        ## Select soft-thresholding power
        cat("==  Selection of soft-thresholding power  == \n")
        load(paste0(outPath,"/tmpData/flag1_dataWithoutOutliers.RData"))
        
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        # Call the network topology analysis function
        sft = pickSoftThreshold(omic, powerVector = powers, verbose = 5)
        # Plot the results:
        sizeGrWindow(9, 5)
        par(mfrow = c(1,2));
        cex1 = 0.9;
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence"));
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,cex=cex1,col="red");
        # this line corresponds to using an R^2 cut-off of h
        abline(h=0.90,col="red")
        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
        dev.copy2pdf(file = paste0(outPath,"/graphs/SfThrSelection.pdf"));
        dev.off()
        
        softPower = readInteger("Enter value for Soft Thresholding Power: ");
        
        load(paste0(outPath,"/tmpData/parameters.RData"))
        save(cutHeightSamples, softPower, file = paste0(outPath,"/tmpData/parameters.RData"))
        flag = 3
        
      }else if(flag==3){
        
        ## Adjacency and Topology overlap calculations
        cat("==  Computations on Adjacency and Topology Overlap Matrix  == \n")
        load(paste0(outPath,"/tmpData/flag1_dataWithoutOutliers.RData"))
        load(paste0(outPath,"/tmpData/parameters.RData"))
        
        adjacency = adjacency(omic, power = softPower);
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency);
        dissTOM = 1-TOM
        # Call the hierarchical clustering function
        geneTree = hclust(as.dist(dissTOM), method = "average");
        # Plot the resulting clustering tree (dendrogram)
        sizeGrWindow(12,9)
        plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
             labels = FALSE, hang = 0.04);
        dev.copy2pdf(file = paste0(outPath,"/graphs/geneClustering.pdf"));
        dev.off()
        
        # Module identification using dynamic tree cut:
        dynamicCut = FALSE # used to loop and select appropriate deepsplit and cutHeight
        while(dynamicCut == FALSE){
          deepSplit = readInteger("deepSplit? ([0,4], default is 2) ")
          cutHeightGeneTree = readInteger("cutHeight? (default is 99 (%)) ")
          dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight = cutHeightGeneTree,
                                      deepSplit = deepSplit, pamRespectsDendro = FALSE,
                                      minClusterSize = 30);
          cat("Genes Clusters:\n")
          print(table(dynamicMods))
          dynamicCut = askConfirmation()
        }
        
        save(TOM, dissTOM, dynamicMods, geneTree, file = paste0(outPath,"/tmpData/flag3_AdjacencyAndTOM.RData"))
        save(cutHeightSamples, softPower, deepSplit, cutHeightGeneTree, file = paste0(outPath,"/tmpData/parameters.RData"))
        flag = 4
        
      }else if(flag==4){
        
        ## Treating the results
        cat("==  Analysis of Genes Clusters  == \n")
        load(paste0(outPath,"/tmpData/flag1_dataWithoutOutliers.RData"))
        load(paste0(outPath,"/tmpData/flag3_AdjacencyAndTOM.RData"))
        load(paste0(outPath,"/tmpData/parameters.RData"))
        
        # Convert numeric lables into colors
        dynamicColors = labels2colors(dynamicMods)
        cat("Modules-associated colors:\n")
        print(table(dynamicColors))
        # Plot the dendrogram and colors underneath
        sizeGrWindow(12,6)
        plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Gene dendrogram and module colors")
        dev.copy2pdf(file = paste0(outPath,"/graphs/ModuleColorsDendrogram.pdf"));
        dev.off()
        
        # Calculate eigengenes
        MEList = moduleEigengenes(omic, colors = dynamicColors)
        MEs = MEList$eigengenes
        # Calculate dissimilarity of module eigengenes
        MEDiss = 1-cor(MEs);
        # Cluster module eigengenes
        METree = hclust(as.dist(MEDiss), method = "average");
        # Plot the result
        sizeGrWindow(7, 6)
        
        cutModules = FALSE # used to loop for selection of satisfying cutheight
        while(cutModules == FALSE){
          plot(METree, main = "Clustering of module eigengenes",
               xlab = "", sub = "")
          cutHeigthModules = readDouble("Enter merging threshold: ")
          # Plot the cut line into the dendrogram
          abline(h=cutHeigthModules, col = "red")
          cutModules = askConfirmation()
        }
        
        dev.copy2pdf(file = paste0(outPath,"/graphs/eigengeneClustering.pdf"));
        dev.off()
        # Call an automatic merging function
        merge = mergeCloseModules(omic, dynamicColors, cutHeight = cutHeigthModules, verbose = 3)
        # The merged module colors
        mergedColors = merge$colors;
        # Eigengenes of the new merged modules:
        mergedMEs = merge$newMEs;
        
        sizeGrWindow(30, 9)
        plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                            c("Dynamic Tree Cut", "Merged dynamic"),
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        dev.copy2pdf(file = paste0(outPath,"/graphs/FinalDendroAndColors.pdf"));
        dev.off()
        
        
        # Rename to moduleColors
        moduleColors = mergedColors
        # Construct numerical labels corresponding to the colors
        colorOrder = c("grey", standardColors(50));
        moduleLabels = match(moduleColors, colorOrder)-1;
        MEs = mergedMEs;
        # Recalculate MEs with color labels
        MEs0 = moduleEigengenes(omic, moduleColors)$eigengenes
        MEs = orderMEs(MEs0)
        
        
        save(MEs, moduleLabels, moduleColors, geneTree, file = paste0(outPath,"/tmpData/flag4_modules.RData"))
        save(cutHeightSamples, softPower, deepSplit, cutHeightGeneTree, cutHeigthModules, 
             file = paste0(outPath,"/tmpData/parameters.RData"))
        flag = 5
        
      }else if(flag==5){
        
        ## Visualise correlation with sample traits
        cat("==  Correlation between molecular data and sample traits  == \n")
        load(paste0(outPath,"/tmpData/flag1_dataWithoutOutliers.RData"))
        load(paste0(outPath,"/tmpData/flag4_modules.RData"))
        
        # Visualize how the samples traits relate to the sample dendrogram
        traits = as.data.frame(samplesTraits[match(rownames(omic), rownames(samplesTraits)),])
        sampleTree2 = hclust(dist(omic), method = "average") # Re-cluster samples
        traitColors = numbers2colors(traits, signed = FALSE); # Convert traits to a color representation: white means low, red means high, grey means missing entry
        # Plot the sample dendrogram and the colors underneath.
        plotDendroAndColors(sampleTree2, traitColors,
                            groupLabels = names(traits),
                            main = "Sample dendrogram and trait heatmap")
        dev.copy2pdf(file = paste0(outPath,"/graphs/traitHeatmap.pdf"));
        dev.off()
        
        moduleTraitCor = cor(MEs, traits, use = "p");
        moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
        sizeGrWindow(50,20)
        # Will display correlations and their p-values
        textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                            signif(moduleTraitPvalue, 1), ")", sep = "");
        dim(textMatrix) = dim(moduleTraitCor)
        par(mar = c(6, 8.5, 3, 3));
        # Display the correlation values within a heatmap plot
        labeledHeatmap(Matrix = moduleTraitCor,
                       xLabels = names(traits),
                       yLabels = names(MEs),
                       ySymbols = names(MEs),
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = 0.15, cex.lab=0.4,
                       zlim = c(-1,1),
                       main = paste("Module-trait relationships"))
        dev.copy2pdf(file = paste0(outPath,"/graphs/CorrelationHeatmap.pdf"));
        dev.off()
        
        save(traits, MEs, moduleLabels, moduleColors, geneTree, file = paste0(outPath,"/tmpData/flag5_correlations.RData"))
        write.table(x=moduleTraitCor, file=paste0(outPath,"/outData/WGCNAcorrelations_modulesDrugs.txt"), sep="\t")
        write.table(x=moduleTraitPvalue, file=paste0(outPath,"/outData/WGCNAcorrelationsPvalues_modulesDrugs.txt"), sep="\t")
        flag = 6
        
      }
    } # end of while flag < 6 
    
  } # end of runWGCNA()
  
  outputWGCNAmodules <- function(path="./", title="WGCNArun", traitsModules){
    {
      ## Extracts the genes corresponding to modules strongly associated with traits of interest, from the result of an already-performed WGCNA
      #
      # @params <path> string giving the path to the parent folder of the WGCNA run folder
      # @params <title> string associated with the WGCNA run; corresponds to the folder in which results can be retrieved and where the gene lists will be output
      # @params <traitsModules> list of vectors; each vector holds the names of the modules strongly associated with a given trait; the name of that trait must be the name of the vector in the list
    } # header of outputWGCNAmodules()
    
    cat("=====  Extracting genes of interest  =====\n")
    
    ## Load the needed WGCNA results
    load(paste0(path,title,"/tmpData/flag1_dataWithoutOutliers.RData"))
    load(paste0(path,title,"/tmpData/flag5_correlations.RData"))
    
    
    ## Compute Gene Significance and Module Membership
    cat("==  Compute Gene Significance and Module Membership  ==\n")
    
    # Module Membership values
    modNames = substring(names(MEs), 3) # names (colors) of the modules
    geneModuleMembership = as.data.frame(cor(omic, MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    names(geneModuleMembership) = paste("MM.", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    
    # Gene Significance per user-provided traits of interest
    screenedTraits = names(traitsModules)
    gtSignificances = list()
    GSPvalues = list()
    for(trait in 1:length(screenedTraits)){
      traitValues = traits[, colnames(traits)==screenedTraits[trait]]
      gtSignificances[[trait]] = as.data.frame(cor(omic, traitValues, use = "p"));
      GSPvalues[[trait]] = as.data.frame(corPvalueStudent(as.matrix(gtSignificances[[trait]]), nSamples));
    }
    names(gtSignificances) = paste0("GS.",screenedTraits)
    names(GSPvalues) = paste0("p.GS.",screenedTraits)
    
    ## Save these matrices in independent file
    save(geneModuleMembership, MMPvalue, gtSignificances, GSPvalues, file = paste0(path,title,"/outData/GeneSignificances_ModulesMemberships.RData"))
    
    
    ## Extract Genes per trait/module association
    cat("==  Extract Genes per trait/module association  ==\n")
    
    for(trait in 1:length(traitsModules)){
      
      cat("Extracting genes lists for ", names(traitsModules)[trait], "...\n", sep="")
      
      nMods = length(traitsModules[[trait]])
      gtSignif = gtSignificances[[trait]]
      gsPvals = GSPvalues[[trait]]
      maxSize = 0 # will be used to determine dimension of the output dataset
      extractedGenes = list() # will be filled with the genes names in each module as well as their significance and associed pvalues
      if(nMods>1){  par(mfrow=c((nMods%/%2 + nMods%%2), 2))  }
      
      if(nMods>6){  par(mar=c(1,1,1,1))  }
      for(mod in 1:length(traitsModules[[trait]])){
        
        module = traitsModules[[trait]][mod]
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(gtSignif[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = paste("Gene significance for", names(traitsModules)[trait]),
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, pch=20, col = module)
        
        moduleList = cbind(colnames(omic)[moduleGenes], gtSignif[moduleGenes, 1], gsPvals[moduleGenes,1])
        colnames(moduleList) = c(paste0(module,".Gene"),paste0(module,".Significance"),paste0(module,".SignificancePvalue"))
        extractedGenes[[mod]] = moduleList[order(moduleList[,2], decreasing=TRUE),]
        if(nrow(moduleList)>maxSize){  maxSize = nrow(moduleList)  }
        
      }
      
      
      cat("Saving results for ", names(traitsModules)[trait], "...\n", sep="")
      
      ## Create a matrix and fill it with the genes information
      genesLists = matrix(data=NA, nrow=maxSize, ncol=length(traitsModules[[trait]])*3)
      cNames = c()
      for(mod in 1:length(traitsModules[[trait]])){
        cNames = c(cNames,colnames(extractedGenes[[mod]]))
        for(gene in 1:nrow(extractedGenes[[mod]])){
          genesLists[gene,(mod-1)*3+1] = extractedGenes[[mod]][gene,1]
          genesLists[gene,(mod-1)*3+2] = extractedGenes[[mod]][gene,2]
          genesLists[gene,(mod-1)*3+3] = extractedGenes[[mod]][gene,3]
        }
      }
      colnames(genesLists) = cNames
      
      ## Save to files
      write.csv2(genesLists, file=paste0(path,title,"/outData/Modules_",names(traitsModules)[trait],".csv"),
                 row.names=FALSE, na="")
      dev.copy2pdf(file = paste0(path,title,"/graphs/Modules_",names(traitsModules)[trait],".pdf"));
      dev.off()
      par(mfrow=c(1,1))
      if(nMods>6){  par(mar=c(5.1, 4.1, 4.1, 2.1))  }
      
    } # end of for(trait in 1:length(traitsModule))
    
    
  } # end of outputWGCNAmodules()
}



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# I - Unbiased approach ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # ===== I.1 - Load Data =====
  load("./tmp/samplesLists.RData")
  load("./tmp/transcFiltered_unbiasedApproach.RData")
  #load("./tmp/ic50_datasets.RData")
  AUCs = read.csv2("./tmp/AUCs_drugExposure.csv", h=T, row.names = 1)
  transc = transcFiltered_unbiasedApproach # for easier handling
  
  
  # ===== I.2 - DASL vs. Clinical data =====
  {
    
    ## I.2.a) All samples ----
    {
      
      transc_all88 = transc[ match(samples_all88, rownames(transc)),  ]
      diagnosis_all88 = diagnosis[ match(samples_all88, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_all88, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_all88$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_all88$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Pathological diagnosis
        gcol = rainbow(length(levels(diagnosis_all88$Pathological_diagnosis)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_all88$Pathological_diagnosis), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
      }
      
      ### LASSO
      lasso_all88_unbiased = runLASSO(transc_all88, diagnosis_all88[, c("WHO_Grade","OS", "Recurrent_Num")])
      print(lasso_all88_unbiased)
      
      ### WGCNA
      runWGCNA("WGCNA_all88_unbiased_clinical", transc_all88, diagnosis_all88[, c("WHO_Grade","OS", "Recurrent_Num")])
      modules = list(c("darkred","darkmagenta","bisque4"), 
                     c("darkred","darkmagenta","bisque4","brown","mediumorchid","darkseagreen4","skyblue3","orangered3"))
      names(modules) = c("WHO_Grade","OS")#,"Recurrent_Num")
      outputWGCNAmodules(path="./results/2019.November_DASLvsClinical/", title="WGCNA_all88_unbiased_clinical", traitsModules = modules)
      
      
      ### TopGO with results
      {
        ## LASSO
        # For Grade ==> inconclusive results
        runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                 investigatedGenes = c("ADAMTS14", "FLNC", "FMN1", "GLI2", "HOXC4", "IGF2BP3", "LRRC2", "NAGS", "PALM2", "ROBO2", "SHOX2", "ZNF324"),
                 pvalues = NULL,#as.double(na.omit(modulesOS_all88_unbiased$darkorange2.SignificancePvalue)),
                 algorithm = "elim", statistic = "fisher")
        # For OS ==> inconclusive results
        runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                 investigatedGenes = c("ATOH8", "C14orf21", "C18orf34", "CALB1", "DACH1", "FBN2", "FOXP2", "KIAA1958", "PARD3", "SYTL4", "THRA"),
                 pvalues = NULL,#as.double(na.omit(modulesOS_all88_unbiased$darkorange2.SignificancePvalue)),
                 algorithm = "elim", statistic = "fisher")
        
        
        ## WGCNA
        modulesOS_all88_unbiased = read.csv2(file="./results/2019.November_DASLvsClinical/Unbiased/WGCNA_all88_unbiased_clinical/outData/Modules_OS.csv", na.strings = "")
        dim(modulesOS_all88_unbiased)
        names(modulesOS_all88_unbiased)
        sapply(modulesOS_all88_unbiased[,c("yellow.Gene","green.Gene","darkorange2.Gene")],
               function(gene) nrow(modulesOS_all88_unbiased)-sum(is.na(gene)))
        
        # For darkorange2 ==> inconclusive results: some interesting results in terms but no matching genes
        darkorange2GO = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                                investigatedGenes = as.character(na.omit(modulesOS_all88_unbiased$darkorange2.Gene)),
                                pvalues = NULL,
                                algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0023019","GO:1902229")), 
               function(term) intersect(term,as.character(na.omit(modulesOS_all88_unbiased$darkorange2.Gene))))
        
        # For green ==> shows some interesting results, in particular on angiogenesis and extracellular matrix structure
        greenGO = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                           investigatedGenes = as.character(na.omit(modulesOS_all88_unbiased$green.Gene)),
                           pvalues = NULL,
                           algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0030198","GO:0001525","GO:1900122","GO:0010469")), 
               function(term) intersect(term,as.character(na.omit(modulesOS_all88_unbiased$green.Gene))))
        
        # For yellow ==> shows some interesting results around DNA structure and replication
        yellowGO = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                                 investigatedGenes = as.character(na.omit(modulesOS_all88_unbiased$yellow.Gene)),
                                 pvalues = NULL,
                                 algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0051301","GO:0032200","GO:0006270","GO:0000281")), 
               function(term) intersect(term,as.character(na.omit(modulesOS_all88_unbiased$yellow.Gene))))
        
      } # end of topGO
      
    }
    
    ## I.2.b) All GBM samples ----
    {
      
      transc_allGBM = transc[ match(samples_allGBM, rownames(transc)),  ]
      diagnosis_allGBM = diagnosis[ match(samples_allGBM, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_allGBM, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_allGBM$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allGBM$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Pathological diagnosis
        gcol = rainbow(length(levels(diagnosis_allGBM$Pathological_diagnosis)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allGBM$Pathological_diagnosis), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
      }
      
      ### Comparison of short-term vs long-term survivors ----
      {
        ## Define groups
        transc_STSurvivors = transc_allGBM[intersect(samples_allSTsurvivors, rownames(transc_allGBM)), ]
        transc_LTSurvivors = transc_allGBM[intersect(samples_allLTsurvivors, rownames(transc_allGBM)), ]
        Nprobes = ncol(transc_allGBM)
        
        ## Run the t-test
        ttestResults = lapply(  1:Nprobes, function(probe) t.test( transc_STSurvivors[, probe], transc_LTSurvivors[, probe] )  )
        names(ttestResults) = colnames(transc_allGBM)
        ttestResults_pvals = sapply(ttestResults, function(probeTest) probeTest$p.value)
        ttestResults_Signif = ttestResults_pvals[ttestResults_pvals<0.05]
        length(ttestResults_Signif)
        ttestResults_Signif
        
        ## Bonferroni multi-test correction
        ttestResults_Bonferroni = (p.adjust(ttestResults_pvals, "bonferroni"))
        ttestResults_BonferroniSignif = ttestResults_Bonferroni[ttestResults_Bonferroni<0.05]
        length(ttestResults_BonferroniSignif)
        ttestResults_BonferroniSignif
        
        ## FDR multi-test correction
        ttestResults_FDR = (p.adjust(ttestResults_pvals, "fdr"))
        ttestResults_FDRSignif = ttestResults_FDR[ttestResults_FDR<0.05]
        length(ttestResults_FDRSignif)
        ttestResults_FDRSignif
        
      }
      
      ### LASSO
      lasso_allGBM_unbiased = runLASSO(transc_allGBM, diagnosis_allGBM[, c("OS", "Recurrent_Num")])
      print(lasso_allGBM_unbiased)
      
      ### WGCNA
      runWGCNA("WGCNA_allGBM_unbiased_clinical", transc_allGBM, diagnosis_allGBM[, c("OS", "Recurrent_Num")])
      modules = list(c("lightgreen","darkolivegreen","orange"), 
                     c("sienna3","darkmagenta"))
      names(modules) = c("OS","Recurrent_Num")
      outputWGCNAmodules(path="./results/2019.August_DASLvsClinical/Unbiased/", title="WGCNA_allGBM_unbiased_clinical", traitsModules = modules)
      
      ### TopGO with results
      {
        ## LASSO
        # For Recurrence ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = c("AHSP", "HBG2", "KIAA1958", "PTCH1", "SNORD65"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        # For OS ==> inconclusive results, too few genes
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = c("ACSM5", "CAB39L", "ZC3H12C"),
                               pvalues = NULL,#as.double(na.omit(modulesOS_all88_unbiased$darkorange2.SignificancePvalue)),
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0007050","GO:0045786","GO:0007049")), function(term) intersect(term,c("ACSM5", "CAB39L", "ZC3H12C")))
        
        
        ## WGCNA
        modulesOS = read.csv2(file="./results/2019.August_DASLvsClinical/Unbiased/WGCNA_allGBM_unbiased_clinical/outData/Modules_OS.csv", na.strings = "")
        dim(modulesOS)
        names(modulesOS)
        sapply(modulesOS[,c("lightgreen.Gene","darkolivegreen.Gene","orange.Gene")],
               function(gene) nrow(modulesOS)-sum(is.na(gene)))
        
        # For lightgreen ==> shows some interesting results, in particular on angiogenesis
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$lightgreen.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0001525")), 
               function(term) intersect(term,as.character(na.omit(modulesOS$lightgreen.Gene))))
        
        # For darkolivegreen ==> shows limited results, on negative regulation of ERK1/2 cascade
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$darkolivegreen.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0070373")), 
               function(term) intersect(term,as.character(na.omit(modulesOS$darkolivegreen.Gene))))
        
        # For orange ==> shows limited results, on insulin receptor pathway
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$orange.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0046626")), 
               function(term) intersect(term,as.character(na.omit(modulesOS$orange.Gene))))
        
        modulesRec = read.csv2(file="./results/2019.August_DASLvsClinical/Unbiased/WGCNA_allGBM_unbiased_clinical/outData/Modules_Recurrent_Num.csv", na.strings = "")
        dim(modulesRec)
        names(modulesRec)
        sapply(modulesRec[,c("sienna3.Gene","darkmagenta.Gene")],
               function(gene) nrow(modulesRec)-sum(is.na(gene)))
        
        # For sienna3 ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesRec$sienna3.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        # For darkolivegreen ==> shows interesting results, in paracrine but especially around angiogenesis
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesRec$darkmagenta.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0016525","GO:0035793","GO:0001570","GO:0038001","GO:0045766")), 
               function(term) intersect(term,as.character(na.omit(modulesRec$darkmagenta.Gene))))
        
      } # end of topGO
      
    }
    
    ## I.2.c) All Primary GBM samples ----
    {
      
      transc_allPrimaryGBM = transc[ match(samples_allPrimaryGBM, rownames(transc)),  ]
      diagnosis_allPrimaryGBM = diagnosis[ match(samples_allPrimaryGBM, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_allPrimaryGBM, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_allPrimaryGBM$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allPrimaryGBM$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))

      }
      
      ### Comparison of short-term vs long-term survivors ----
      {
        ## Define groups
        transc_STSurvivors = transc_allPrimaryGBM[intersect(samples_allSTsurvivors, rownames(transc_allPrimaryGBM)), ]
        transc_LTSurvivors = transc_allPrimaryGBM[intersect(samples_allLTsurvivors, rownames(transc_allPrimaryGBM)), ]
        Nprobes = ncol(transc_allPrimaryGBM)
        
        ## Run the t-test
        ttestResults = lapply(  1:Nprobes, function(probe) t.test( transc_STSurvivors[, probe], transc_LTSurvivors[, probe] )  )
        names(ttestResults) = colnames(transc_allPrimaryGBM)
        ttestResults_pvals = sapply(ttestResults, function(probeTest) probeTest$p.value)
        ttestResults_Signif = ttestResults_pvals[ttestResults_pvals<0.05]
        length(ttestResults_Signif)
        ttestResults_Signif
        
        ## Bonferroni multi-test correction
        ttestResults_Bonferroni = (p.adjust(ttestResults_pvals, "bonferroni"))
        ttestResults_BonferroniSignif = ttestResults_Bonferroni[ttestResults_Bonferroni<0.05]
        length(ttestResults_BonferroniSignif)
        ttestResults_BonferroniSignif
        
        ## FDR multi-test correction
        ttestResults_FDR = (p.adjust(ttestResults_pvals, "fdr"))
        ttestResults_FDRSignif = ttestResults_FDR[ttestResults_FDR<0.05]
        length(ttestResults_FDRSignif)
        ttestResults_FDRSignif
        
      }
      
      ### LASSO
      lasso_allPrimaryGBM_unbiased = runLASSO(transc_allPrimaryGBM, diagnosis_allPrimaryGBM[, c("OS")])
      print(lasso_allPrimaryGBM_unbiased)
      
      ### WGCNA
      runWGCNA("WGCNA_allPrimaryGBM_unbiased_clinical", transc_allPrimaryGBM, diagnosis_allPrimaryGBM[, c("OS", "Recurrent_Num")])
      modules = list(c("tan","grey60","saddlebrown","sienna3","lightyellow"))
      names(modules) = c("OS")
      outputWGCNAmodules(path="./results/2019.August_DASLvsClinical/Unbiased/", title="WGCNA_allPrimaryGBM_unbiased_clinical", traitsModules = modules)
      
      ### TopGO with results
      {
        ## LASSO
        # For OS ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = c("EPHA4", "IRX2", "NFYB", "OMD", "SATB1", "TMEM80", "WDR16"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        ## WGCNA
        modulesOS = read.csv2(file="./results/2019.August_DASLvsClinical/Unbiased/WGCNA_allPrimaryGBM_unbiased_clinical/outData/Modules_OS.csv", na.strings = "")
        names(modulesOS)
        dim(modulesOS)
        sapply(modulesOS[,c("tan.Gene","grey60.Gene","saddlebrown.Gene","sienna3.Gene","lightyellow.Gene")],
               function(gene) nrow(modulesOS)-sum(is.na(gene)))
        
        # For tan ==> inconclusive results, related to cilium assembly and motility
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$tan.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        # For grey60 ==> shows some interesting results, on G protein-coupled receptor signaling pathway, ERK1/2 cascade and TOR signalling
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$grey60.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0106071","GO:0106072","GO:0070373","GO:0031929")), 
               function(term) intersect(term,as.character(na.omit(modulesOS$grey60.Gene))))
        
        # For saddlebrown ==> inconclusive, focused on general translational/homeostasis/catabolic processes
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$saddlebrown.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        # For sienna3 ==> shows interesting results around DNA repair
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$sienna3.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0006297","GO:1900264","GO:0006260","GO:0000724","GO:0006270")), 
               function(term) intersect(term,as.character(na.omit(modulesOS$sienna3.Gene))))
        
        # For lightyellow ==> shows interesting results, around GTPase activity, response to drug, and IL6 production
        topGOresult = runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                               investigatedGenes = as.character(na.omit(modulesOS$lightyellow.Gene)),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0043547","GO:0042493","GO:0032715")), 
               function(term) intersect(term,as.character(na.omit(modulesOS$lightyellow.Gene))))
        
      } # end of topGO
      
    }
    
    ## I.2.d) All Recurrent GBM samples ----
    {
      
      transc_allRecurrentGBM = transc[ match(samples_allRecurrentGBM, rownames(transc)),  ]
      diagnosis_allRecurrentGBM = diagnosis[ match(samples_allRecurrentGBM, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_allRecurrentGBM, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_allRecurrentGBM$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allRecurrentGBM$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
      }
      
      ### LASSO
      lasso_allRecurrentGBM_unbiased = runLASSO(transc_allRecurrentGBM, diagnosis_allRecurrentGBM[, c("OS")])
      print(lasso_allRecurrentGBM_unbiased)
      
      ### WGCNA
      runWGCNA("WGCNA_allRecurrentGBM_unbiased_clinical", transc_allRecurrentGBM, diagnosis_allRecurrentGBM[, c("OS", "Recurrent_Num")])
      
    }

    
  }# end of I.2 - DASL vs. Clinical data
  
  
  # ===== I.3 - DASL vs. AUCs =====
  {
    
    ## I.3.a) All 45 samples ----
    {
      transc_allExposed = transc[ match(samples_drugExposedAll, rownames(transc)),  ]
      AUCs_allExposed = AUCs[ match(samples_drugExposedAll, rownames(AUCs)),  ]
      cat(dim(transc_allExposed), dim(AUCs_allExposed))

      ### LASSO
      lasso_allExposed_unbiased = runLASSO(transc_allExposed, AUCs_allExposed)
      dim(lasso_allExposed_unbiased)
      write.csv2(lasso_allExposed_unbiased, "./results/2019.November/allExposed_AUCs/LASSO/results_allExposed_AUCs.csv")
      
      ### WGCNA
      runWGCNA(path="./results/2019.November/allExposed_AUCs/WGCNA/", title="WGCNA_AUCallExposed1_unbiased", transc_allExposed, AUCs_allExposed[,1:44], flag=5)
      runWGCNA(path="./results/2019.November/allExposed_AUCs/WGCNA/", title="WGCNA_AUCallExposed2_unbiased", transc_allExposed, AUCs_allExposed[,45:97], flag=5)
      
      modules1 = list(Hydroxyurea = c("greenyellow"), Fluorouracil = c("mediumpurple3"), Mercaptopurine = c("brown4","darkgrey"),
                      Mechlorethamine.hydrochloride = c("brown", "brown4"), Thiotepa = c("brown4", "darkgrey"), 
                      Busulfan = c("skyblue3"), Methoxsalen = c("floralwhite", "cyan", "mediumpurple3"), Decitabine = c("brown4"),
                      Carmustine = c("skyblue3", "brown4", "darkgrey"), Uracil.mustard = c("floralwhite", "lightyellow", "brown4"),
                      Procarbazine.hydrochloride = c("orange", "brown4"), Tretinoin = c("darkgrey"), 
                      Dexrazoxane = c("floralwhite"), Bortezomib = c("lightcyan1"), Celecoxib = c("greenyellow"),
                      Pemetrexed = c("saddlebrown"), Methotrexate = c("saddlebrown"))
      modules2 = list(Topotecan.hydrochloride = c("lightcyan1", "brown", "mediumpurple3"), Raloxifene = c("green"), 
                      Pralatrexate = c("brown","greenyellow"), Ixabepilone = c("lightyellow", "skyblue3", "darkolivegreen", "lightgreen"), 
                      Valrubicin = c("cyan"), Dactinomycin = c("lightsteelblue1", "darkolivegreen", "green", "darkorange2", "darkgrey"),
                      Bleomycin.sulfate = c("skyblue3", "darkgreen", "mediumpurple3", "greenyellow"),
                      Gemcitabine.hydrochloride = c("lightyellow"),
                      Arsenic.trioxide = c("lightyellow", "plum1", "darkslateblue", "brown4", "darkgrey"),
                      Oxaliplatin = c("lightcyan1"), Amiodarone.hydrochloride = c("darkgreen", "green"))
      
      outputWGCNAmodules("./results/2019.November/allExposed_AUCs/WGCNA/", "WGCNA_AUCallExposed1_unbiased", modules1)
      outputWGCNAmodules("./results/2019.November/allExposed_AUCs/WGCNA/", "WGCNA_AUCallExposed2_unbiased", modules2)
      

      ### TopGO with results
      {
        ## LASSO
        # For OS ==> inconclusive results
        for(drug in rownames(lasso_allExposed_unbiased)){
          if(sum(!is.na(lasso_allExposed_unbiased[drug,]))>10){
            lastGene_idx = sum(!is.na(lasso_allExposed_unbiased[drug,]))
            cat(drug, lastGene_idx, sum(!is.na(lasso_allExposed_unbiased[drug,])), '\n')
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = ifelse(lasso_allExposed_unbiased[drug,1] == "(Intercept)", 
                                                lasso_allExposed_unbiased[drug,2:lastGene_idx],
                                                lasso_allExposed_unbiased[drug,1:lastGene_idx]),
                     outPath = "./results/2019.November/unbiasedApproach/allExposed_AUCs/LASSO/", title = paste0("TopGO_",drug),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
        }
        
        ## WGCNA
        for(drug in names(modules1)){
          
          modules = read.csv2(file=paste0("./results/2019.November/allExposed_AUCs/WGCNA/WGCNA_AUCallExposed1_unbiased/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modules1[[drug]]){
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/allExposed_AUCs/WGCNA/WGCNA_AUCallExposed1_unbiased/outData/", 
                     title = paste0("TopGO_WGCNA_allExposed_AUCs_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
        for(drug in names(modules2)){
          
          modules = read.csv2(file=paste0("./results/2019.November/allExposed_AUCs/WGCNA/WGCNA_AUCallExposed2_unbiased/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modules2[[drug]]){
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/allExposed_AUCs/WGCNA/WGCNA_AUCallExposed2_unbiased/outData/", 
                     title = paste0("TopGO_WGCNA_allExposed_AUCs_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
      } # end of topGO
    }
    
    ## I.3.b) All Primary GBM samples that have AUCs available ----
    {
      
      transc_PrimaryExposed = transc[ match(samples_drugExposedPrimary, rownames(transc)),  ]
      AUCs_PrimaryExposed = AUCs[ match(samples_drugExposedPrimary, rownames(AUCs)),  ]
      cat(dim(transc_PrimaryExposed), dim(AUCs_PrimaryExposed))
      
      ### LASSO
      lasso_PrimaryExposed_unbiased = runLASSO(transc_PrimaryExposed, AUCs_PrimaryExposed)
      dim(lasso_PrimaryExposed_unbiased)
      write.csv2(lasso_PrimaryExposed_unbiased, "./results/2019.November/PrimaryExposed_AUCs/LASSO/results_PrimaryExposed_AUCs.csv")
      
      ### WGCNA
      runWGCNA(path="./results/2019.November/PrimaryExposed_AUCs/WGCNA/", title="WGCNA_AUCPrimaryExposed1_unbiased", transc_PrimaryExposed, AUCs_PrimaryExposed[,1:44])
      runWGCNA(path="./results/2019.November/PrimaryExposed_AUCs/WGCNA/", title="WGCNA_AUCPrimaryExposed2_unbiased", transc_PrimaryExposed, AUCs_PrimaryExposed[,45:97])
      
      modules1 = list(Hydroxyurea = c("indianred4"), Mercaptopurine = c("blue2"), Thiotepa = c("yellow4"), 
                      Busulfan = c("lightpink4", "plum3"), 
                      Methoxsalen = c("darkgreen", "tan", "thistle2", "brown4", "lightcyan1", "violet", "orange"), 
                      Lomustine = c("greenyellow", "yellow4", "skyblue2", "brown2", "darkolivegreen4", "mediumorchid"),
                      Azacitidine = c("magenta", "coral1", "indianred4"), Decitabine = c("orange"),
                      Carmustine = c("antiquewhite4", "navajowhite2", "firebrick4", "indianred4"), 
                      Uracil.mustard = c("darkslateblue", "violet", "palevioletred3"),
                      Procarbazine.hydrochloride = c("darkred", "bisque4", "green", "maroon"), 
                      Dexrazoxane = c("darkgreen", "orangered3", "yellow4", "bisque4", "firebrick4", "magenta", "darkolivegreen4"),
                      Mitomycin = c("darkgreen", "yellow4"), Pipobroman = c("darkgreen", "yellow4", "darkolivegreen"),
                      Fludarabine.phosphate = c("darkturquoise", "lightsteelblue", "skyblue1", "firebrick4", "steelblue"),
                      Bortezomib = c("greenyellow"), Celecoxib = c("pink"), Sunitinib = c("greenyellow", "skyblue1"),
                      Mitoxantrone = c("magenta"), Pemetrexed = c("darkturquoise", "skyblue1"), Gefitinib = c("firebrick4"), 
                      Methotrexate = c("lightsteelblue", "skyblue2", "darkred", "navajowhite2", "brown2", "darkolivegreen4"))
      modules2 = list(Topotecan.hydrochloride = c("darkgreen", "greenyellow", "yellow4", "maroon", "brown2", "darkolivegreen", "mediumorchid"),
                      Pazopanib.hydrochloride = c("antiquewhite4", "coral2"), Imatinib = c("orangered4"),  
                      Sorafenib = c("steelblue"), Raloxifene = c("orangered4"),
                      Pralatrexate = c("coral2", "brown2", "darkolivegreen4", "royalblue", "darkslateblue","yellow4", "skyblue2", 
                                       "navajowhite2", "mediumpurple2", "maroon", "blue2", "palevioletred2", "lightsteelblue1", 
                                       "coral1", "darkolivegreen"), 
                      Enzalutamide = c("orangered4"), Vemurafenib = c("darkviolet"), 
                      Idarubicin.hydrochloride = c("brown2", "darkolivegreen", "mediumorchid"), Nilotinib = c("yellow4", "darkolivegreen"),
                      Ixabepilone = c("mediumorchid"), Ponatinib = c("paleturquoise"), Doxorubicin.hydrochloride = c("lightcyan1"),
                      Tamoxifen.citrate = c("antiquewhite4"), Epirubicin.hydrochloride = c("lightcyan1"), 
                      Lapatinib = c("firebrick4", "violet"), Irinotecan.hydrochloride = c("maroon"),
                      Valrubicin = c("thistle2", "darkolivegreen"), Carfilzomib = c("orangered3", "coral1"), Paclitaxel = c("greenyellow"),
                      Dactinomycin = c("darkturquoise", "skyblue1", "coral2", "orangered4", "lightcyan", "maroon", "darkslateblue",
                                       "tan", "antiquewhite4", "green", "blue2", "brown2", "lightsteelblue1", "royalblue"),
                      Bleomycin.sulfate = c("darkgreen", "darkturquoise", "skyblue1"), Vinorelbine.tartrate = c("palevioletred3"),
                      Temsirolimus = c("coral1", "darkolivegreen4", "mediumorchid"), 
                      Oxaliplatin = c("greenyellow", "yellow4", "bisque4", "green", "maroon", "royalblue", "plum3"), 
                      Melphalan.hydrochloride = c("greenyellow", "yellow4", "maroon", "darkolivegreen"),
                      Amiodarone.hydrochloride = c("magenta", "orangered4", "palevioletred3"))
      
      outputWGCNAmodules("./results/2019.November/PrimaryExposed_AUCs/WGCNA/", "WGCNA_AUCPrimaryExposed1_unbiased", modules1)
      outputWGCNAmodules("./results/2019.November/PrimaryExposed_AUCs/WGCNA/", "WGCNA_AUCPrimaryExposed2_unbiased", modules2)
      
      
      ### TopGO with results
      {
        ## LASSO
        for(drug in rownames(lasso_PrimaryExposed_unbiased)){
          if(sum(!is.na(lasso_PrimaryExposed_unbiased[drug,]))>10){
            lastGene_idx = sum(!is.na(lasso_PrimaryExposed_unbiased[drug,]))
            cat(drug, lastGene_idx, sum(!is.na(lasso_PrimaryExposed_unbiased[drug,])), '\n')
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = ifelse(lasso_PrimaryExposed_unbiased[drug,1] == "(Intercept)", 
                                                lasso_PrimaryExposed_unbiased[drug,2:lastGene_idx],
                                                lasso_PrimaryExposed_unbiased[drug,1:lastGene_idx]),
                     outPath = "./results/2019.November/unbiasedApproach/PrimaryExposed_AUCs/LASSO/", title = paste0("TopGO_LASSO_PrimaryExposed_AUCs_",drug),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
        }
        
        ## WGCNA
        for(drug in names(modules1)){
          
          modules = read.csv2(file=paste0("./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed1_unbiased/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modules1[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed1_unbiased/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
        for(drug in names(modules2)){
          
          modules = read.csv2(file=paste0("./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed2_unbiased/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modules2[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed2_unbiased/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
      } # end of topGO
      
    }
    
    ## I.3.c) All Recurrent GBM samples that have AUCs available ----
    {
      
      transc_RecurrentExposed = transc[ match(samples_drugExposedRecurrent, rownames(transc)),  ]
      AUCs_RecurrentExposed = AUCs[ match(samples_drugExposedRecurrent, rownames(AUCs)),  ]
      cat(dim(transc_RecurrentExposed), dim(AUCs_RecurrentExposed))
      
      ### LASSO
      lasso_RecurrentExposed_unbiased = runLASSO(transc_RecurrentExposed, AUCs_RecurrentExposed[,1:96])
      dim(lasso_RecurrentExposed_unbiased)
      write.csv2(lasso_RecurrentExposed_unbiased, "./results/2019.November/RecurrentExposed_AUCs/LASSO/results_RecurrentExposed_AUCs.csv")
      
      ### WGCNA
      runWGCNA(path="./results/2019.November/RecurrentExposed_AUCs/WGCNA/", title="WGCNA_AUCRecurrentExposed1_unbiased", transc_RecurrentExposed, AUCs_RecurrentExposed[,1:44])
      runWGCNA(path="./results/2019.November/RecurrentExposed_AUCs/WGCNA/", title="WGCNA_AUCRecurrentExposed2_unbiased", transc_RecurrentExposed, AUCs_RecurrentExposed[,45:97])
      
      modules1 = list(Hydroxyurea = c("whitesmoke"), Fluorouracil = c("royalblue3"), Thioguanine = c("royalblue", "darkgoldenrod4", "floralwhite"),
                      Mechlorethamine.hydrochloride = c("deeppink", "aliceblue", "coral1"), Thiotepa = c("aliceblue", "chocolate3"), 
                      Aminolevulinic.acid.hydrochloride = c("darkmagenta", "tomato", "mediumorchid4", "magenta2"),
                      Temozolomide = c("moccasin", "plum4", "deeppink", "darkseagreen3", "bisque4", "darkseagreen4", "darkolivegreen4"), 
                      Busulfan = c("aliceblue", "yellow2"), Carmustine = c("chocolate3"), Uracil.mustard = c("aliceblue", "coral1"),
                      Procarbazine.hydrochloride = c("lightskyblue3", "plum4", "orangered", "darkseagreen3", "lightpink1",  
                                                     "indianred1", "brown2", "firebrick2", "navajowhite1", "floralwhite", 
                                                     "lightskyblue4", "blue1", "darkolivegreen1", "darkolivegreen"), 
                      Streptozocin = c("tan3", "indianred2"), Cisplatin = c("mistyrose", "darkseagreen3", "coral2", "brown3", "magenta2"),
                      Tretinoin = c("coral1", "green"), Nelarabine = c("bisque4", "plum3", "royalblue", "lightskyblue4"),
                      Mitomycin = c("lightskyblue3"), Mitotane = c("thistle4"), 
                      Bendamustine.hydrochloride = c("darkslateblue", "tan3", "blueviolet", "darkgoldenrod4", "floralwhite"),
                      Fludarabine.phosphate = c("ivory"), Bortezomib = c("indianred4"), Axitinib = c("antiquewhite1", "darkolivegreen"),
                      Mitoxantrone = c("darkolivegreen1"), Pemetrexed = c("palevioletred", "whitesmoke"), 
                      Vismodegib = c("darkslateblue", "aliceblue", "coral1", "yellow2"), Methotrexate = c("palevioletred", "whitesmoke"))
      modules2 = list(Dasatinib = c("yellow2"), Imatinib = c("lightslateblue"), Raloxifene = c("darkseagreen3"), Afatinib = c("firebrick2"),
                      Pazopanib.hydrochloride = c("darkseagreen3", "brown2", "magenta2", "orangered4", "navajowhite1"), 
                      Pralatrexate = c("cornflowerblue", "palevioletred", "tomato", "chocolate3", "brown1", "blue", "indianred1", 
                                       "indianred3", "blue3"), 
                      Vandetanib = c("firebrick2"), Vemurafenib = c("tan3", "chocolate4", "darkolivegreen"), 
                      Romidepsin = c("lightskyblue3"), Omacetaxine.mepesuccinate = c("lightskyblue3"), Ponatinib = c("tan3", "thistle4"), 
                      Cabozantinib = c("antiquewhite2", "mediumpurple4", "magenta2"), Doxorubicin.hydrochloride = c("slateblue1", "thistle"),
                      Etoposide = c("plum3", "yellow2"), Epirubicin.hydrochloride = c("coral1"), Estramustine.phosphate.sodium = c("antiquewhite1"),
                      Lapatinib = c("blue4"), Irinotecan.hydrochloride = c("darkolivegreen"), Carfilzomib = c("antiquewhite2", "firebrick2"), 
                      Fulvestrant = c("lightskyblue3", "royalblue", "darkseagreen1", "floralwhite", "indianred2"),
                      Dabrafenib.mesylate = c("mistyrose", "royalblue3"), Teniposide = c("darkolivegreen"),
                      Paclitaxel = c("tomato2", "skyblue2", "tan4", "indianred3", "yellow2"),
                      Vinblastine.sulfate = c("brown2", "mediumpurple4", "blue4", "yellow2"),
                      Vincristine.sulfate = c("darkslateblue", "mistyrose", "royalblue3", "indianred4", "plum4", "blueviolet", "yellow2", "floralwhite"),
                      Sirolimus = c("palevioletred3", "antiquewhite2"), Plicamycin = c("lightskyblue3", "firebrick2"),
                      Bleomycin.sulfate = c("royalblue3", "lightslateblue", "blue4"), Vinorelbine.tartrate = c("tomato", "lightskyblue4", "mediumpurple2"),
                      Temsirolimus = c("lightskyblue3"), Gemcitabine.hydrochloride = c("indianred4"),
                      Dacarbazine = c("palevioletred", "salmon1", "thistle", "royalblue3", "blueviolet", "blue3", "blue4", "antiquewhite1"),
                      Arsenic.trioxide = c("palevioletred", "thistle", "chocolate3", "lavenderblush"), Carboplatin = c("coral1"),
                      Erlotinib.hydrochloride = c("thistle4"), Melphalan.hydrochloride = c("aliceblue"))
      
      modules = list()
      WGCNAcorrelation = read.table("./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed1_unbiased/outData/WGCNAcorrelations_modulesDrugs.txt", h=T)
      for(drug in colnames(WGCNAcorrelation)){
        mods = rownames(WGCNAcorrelation)[abs(WGCNAcorrelation[,drug]) > 0.7]
        if(length(mods) > 0){  modules[[drug]] = gsub("ME", "", mods)  }
      }
      
      outputWGCNAmodules("./results/2019.November/RecurrentExposed_AUCs/WGCNA/", "WGCNA_AUCRecurrentExposed1_unbiased", modules1)
      outputWGCNAmodules("./results/2019.November/RecurrentExposed_AUCs/WGCNA/", "WGCNA_AUCRecurrentExposed2_unbiased", modules2)
      
      
      ### TopGO with results
      {
        ## LASSO
        # For OS ==> inconclusive results
        for(drug in rownames(lasso_RecurrentExposed_unbiased)){
          if(sum(!is.na(lasso_RecurrentExposed_unbiased[drug,]))>10){
            lastGene_idx = sum(!is.na(lasso_RecurrentExposed_unbiased[drug,]))
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = lasso_RecurrentExposed_unbiased[drug,2:lastGene_idx],
                     outPath = "./results/2019.November/RecurrentExposed_AUCs/LASSO/", title = paste0("TopGO_LASSO_RecurrentExposed_AUCs_",drug),
                     pvalues = NULL,
                     algorithm = "elim", statistic = "fisher")
          }
        }
        
        ## WGCNA
        for(drug in names(modules1)){
          
          modules = read.csv2(file=paste0("./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed1_unbiased/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modules1[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed1_unbiased/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
        for(drug in names(modules2)){
          
          modules = read.csv2(file=paste0("./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed2_unbiased/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modules2[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_unbiasedApproach,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed2_unbiased/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
      } # end of topGO
    }
    
    
  }# end of I.3 - DASL vs. AUCs
  
  
} # end of I - Unbiased approach




# oooooooooooooooooooooooooooooooooooooooooooooooooo
# II - Cancer Genes - focused approach ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # ===== II.1 - Load Data =====
  load("./tmp/samplesLists.RData")
  load("./tmp/transcFiltered_cancerGenes.RData")
  #load("./tmp/ic50_datasets.RData")
  AUCs = read.csv2("./tmp/AUCs_drugExposure.csv", h=T, row.names = 1)
  transc = transcFiltered_cancerGenes # for easier handling
  
  
  # ===== II.2 - DASL vs. Clinical data =====
  {
    
    ## II.2.a) All samples ----
    {
      
      transc_all88 = transc[ match(samples_all88, rownames(transc)),  ]
      diagnosis_all88 = diagnosis[ match(samples_all88, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_all88, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_all88$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_all88$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Pathological diagnosis
        gcol = rainbow(length(levels(diagnosis_all88$Pathological_diagnosis)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_all88$Pathological_diagnosis), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
      }
      
      ### LASSO
      lasso_all88_cancerGenes = runLASSO(transc_all88, diagnosis_all88[, c("WHO_Grade","OS", "Recurrent_Num")])
      print(lasso_all88_cancerGenes)
      
      ### WGCNA
      runWGCNA("WGCNA_all88_cancerGenes_clinical", transc_all88, diagnosis_all88[, c("WHO_Grade","OS", "Recurrent_Num")])
      modules = list(c("red","blue","magenta","black","brown","lightcyan"))
      names(modules) = c("OS")
      outputWGCNAmodules(path="./results/2019.August_DASLvsClinical/CancerGenes/", title="WGCNA_all88_CancerGenes_clinical", traitsModules = modules)
      
      ### TopGO with results
      {
        ## LASSO
        # For OS ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = c("AMH", "BCR", "CDKN2B", "FLJ46838", "IGF1R", "PRKX", "SMAD9", "TSHR"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0046777","GO:0051389")), 
               function(term) intersect(term,c("AMH", "BCR", "CDKN2B", "FLJ46838", "IGF1R", "PRKX", "SMAD9", "TSHR")))
        
        # For Recurrence ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = c("BCOR", "C14orf173", "CD28", "DICER1", "HBB", "HBD", "MAX", "UTX"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0015671","GO:0098869")), 
               function(term) intersect(term,c("BCOR", "C14orf173", "CD28", "DICER1", "HBB", "HBD", "MAX", "UTX")))
        
        # For Grade ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = c("BCL11A", "CDC6", "FLNC", "ITGA2", "SOX11", "SP8"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        
        ## WGCNA
        modulesOS = read.csv2(file="./results/2019.August_DASLvsClinical/CancerGenes/WGCNA_all88_CancerGenes_clinical/outData/Modules_OS.csv", na.strings = "")
        names(modulesOS)
        dim(modulesOS)
        sapply(modulesOS[,c("red.Gene","blue.Gene","magenta.Gene","black.Gene","brown.Gene","lightcyan.Gene")],
               function(gene) nrow(modulesOS)-sum(is.na(gene)))
        
        # For red ==> shows some interesting results, around cellular response to growth factor, differentiation, proliferation
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$red.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0071363","GO:0030154","GO:0048665","GO:0008284")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        # For blue ==> shows some interesting results, on insulin like growth factor pathway and inflammatory response
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$blue.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0045944","GO:0006954","GO:0048009")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        # For magenta ==> shows very interesting results, focusing on signalling pathways and protein kinase activity
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$magenta.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0070373","GO:0071901","GO:0045725","GO:0007169","GO:0051387","GO:0032147")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        # For black ==> shows some interesting results, on TOR signaling, response to DNA damage and growth factor stimulus 
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$black.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0006974","GO:0032008","GO:1990090","GO:0048589")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        # For brown ==> shows very interesting results, focusing on cell cycle and DNA replication and repair
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$brown.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL, topNods=20,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0006270","GO:0000082","GO:0051301","GO:0000086","GO:0070317","GO:0000723","GO:0036297","GO:0000724","GO:0006977","GO:0006260","GO:1901796")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        # For lightcyan ==> shows some interesting results, focusing on angiogenesis in particular
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$lightcyan.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0045766","GO:0072104","GO:0043536","GO:0001974","GO:0070374","GO:0048663","GO:0010594")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        
      } # end of topGO
      
    }
    
    ## II.2.b) All GBM samples ----
    {
      
      transc_allGBM = transc[ match(samples_allGBM, rownames(transc)),  ]
      diagnosis_allGBM = diagnosis[ match(samples_allGBM, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_allGBM, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_allGBM$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allGBM$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
      }
      
      ### Comparison of short-term vs long-term survivors ----
      {
        ## Define groups
        transc_STSurvivors = transc_allGBM[intersect(samples_allSTsurvivors, rownames(transc_allGBM)), ]
        transc_LTSurvivors = transc_allGBM[intersect(samples_allLTsurvivors, rownames(transc_allGBM)), ]
        Nprobes = ncol(transc_allGBM)
        
        ## Run the t-test
        ttestResults = lapply(  1:Nprobes, function(probe) t.test( transc_STSurvivors[, probe], transc_LTSurvivors[, probe] )  )
        names(ttestResults) = colnames(transc_allGBM)
        ttestResults_pvals = sapply(ttestResults, function(probeTest) probeTest$p.value)
        ttestResults_Signif = ttestResults_pvals[ttestResults_pvals<0.05]
        length(ttestResults_Signif)
        ttestResults_Signif
        
        ## Bonferroni multi-test correction
        ttestResults_Bonferroni = (p.adjust(ttestResults_pvals, "bonferroni"))
        ttestResults_BonferroniSignif = ttestResults_Bonferroni[ttestResults_Bonferroni<0.05]
        length(ttestResults_BonferroniSignif)
        ttestResults_BonferroniSignif
        
        ## FDR multi-test correction
        ttestResults_FDR = (p.adjust(ttestResults_pvals, "fdr"))
        ttestResults_FDRSignif = ttestResults_FDR[ttestResults_FDR<0.05]
        length(ttestResults_FDRSignif)
        ttestResults_FDRSignif
        
      }
      
      ### LASSO
      lasso_allGBM_cancerGenes = runLASSO(transc_allGBM, diagnosis_allGBM[, c("OS", "Recurrent_Num")])
      print(lasso_allGBM_cancerGenes)
      
      ### WGCNA
      runWGCNA("WGCNA_allGBM_cancerGenes_clinical", transc_allGBM, diagnosis_allGBM[, c("OS", "Recurrent_Num")])
      modules = list(c("cyan","black"), c("purple"))
      names(modules) = c("OS","Recurrent_Num")
      outputWGCNAmodules(path="./results/2019.August_DASLvsClinical/CancerGenes/", title="WGCNA_allGBM_CancerGenes_clinical", traitsModules = modules)
      
      ### TopGO with results
      {
        ## LASSO
        # For OS ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = c("CR1", "P11", "TSHR"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        # For Recurrence ==> shows some very limited results around autophagy, Notch receptor, PI3K; due to EP300, PIK3CA, PBX1, ERCC4, CD28, IGFR1
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = c("BCOR", "CACNG6", "CD28", "CNTF", "DICER1", "EFNA2", "EP300", "ERCC4", "HBB", "HBD", "IGF1R", "PBX1", "PIK3CA", "PTCH1", "TSHR", "UTX"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0007221","GO:0010506","GO:0014065")), 
               function(term) intersect(term,c("BCOR", "CACNG6", "CD28", "CNTF", "DICER1", "EFNA2", "EP300", "ERCC4", "HBB", "HBD", "IGF1R", "PBX1", "PIK3CA", "PTCH1", "TSHR", "UTX")))
        
        ## WGCNA
        modulesOS = read.csv2(file="./results/2019.August_DASLvsClinical/CancerGenes/WGCNA_allGBM_CancerGenes_clinical/outData/Modules_OS.csv", na.strings = "")
        names(modulesOS)
        dim(modulesOS)
        sapply(modulesOS[,c("cyan.Gene","black.Gene")],
               function(gene) nrow(modulesOS)-sum(is.na(gene)))
        
        # For cyan ==> shows some interesting results but all over the place: signalling pathways, DNA damage response...
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$cyan.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0045944","GO:0051387","GO:0034260","GO:0070373","GO:0009967","GO:0000077","GO:0071479")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        # For black ==> shows some interesting results, on cytokine signalling pathway, cell proliferation and immune response
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$black.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL, topNods=20,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0019221","GO:0008285","GO:0006955","GO:0071364","GO:0150077")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        modulesRec = read.csv2(file="./results/2019.August_DASLvsClinical/CancerGenes/WGCNA_allGBM_CancerGenes_clinical/outData/Modules_Recurrent_Num.csv", na.strings = "")
        names(modulesRec)
        dim(modulesRec)
        sapply(modulesRec[,c("purple.Gene")],
               function(gene) nrow(modulesRec)-sum(is.na(gene)))
        
        # For purple ==> shows limited results with TOR signalling
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesRec$purple.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0032007")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
        
      } # end of topGO
      
    }
    
    ## II.2.c) All Primary GBM samples ----
    {
      
      transc_allPrimaryGBM = transc[ match(samples_allPrimaryGBM, rownames(transc)),  ]
      diagnosis_allPrimaryGBM = diagnosis[ match(samples_allPrimaryGBM, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_allPrimaryGBM, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_allPrimaryGBM$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allPrimaryGBM$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
      }
      
      ### Comparison of short-term vs long-term survivors ----
      {
        ## Define groups
        transc_STSurvivors = transc_allPrimaryGBM[intersect(samples_allSTsurvivors, rownames(transc_allPrimaryGBM)), ]
        transc_LTSurvivors = transc_allPrimaryGBM[intersect(samples_allLTsurvivors, rownames(transc_allPrimaryGBM)), ]
        Nprobes = ncol(transc_allPrimaryGBM)
        
        ## Run the t-test
        ttestResults = lapply(  1:Nprobes, function(probe) t.test( transc_STSurvivors[, probe], transc_LTSurvivors[, probe] )  )
        names(ttestResults) = colnames(transc_allPrimaryGBM)
        ttestResults_pvals = sapply(ttestResults, function(probeTest) probeTest$p.value)
        ttestResults_Signif = ttestResults_pvals[ttestResults_pvals<0.05]
        length(ttestResults_Signif)
        ttestResults_Signif
        
        ## Bonferroni multi-test correction
        ttestResults_Bonferroni = (p.adjust(ttestResults_pvals, "bonferroni"))
        ttestResults_BonferroniSignif = ttestResults_Bonferroni[ttestResults_Bonferroni<0.05]
        length(ttestResults_BonferroniSignif)
        ttestResults_BonferroniSignif
        
        ## FDR multi-test correction
        ttestResults_FDR = (p.adjust(ttestResults_pvals, "fdr"))
        ttestResults_FDRSignif = ttestResults_FDR[ttestResults_FDR<0.05]
        length(ttestResults_FDRSignif)
        ttestResults_FDRSignif
        
      }
      
      ### LASSO
      lasso_allPrimaryGBM_cancerGenes = runLASSO(transc_allPrimaryGBM, diagnosis_allPrimaryGBM[, c("OS")])
      print(lasso_allPrimaryGBM_cancerGenes)
      
      ### WGCNA
      runWGCNA("WGCNA_allPrimaryGBM_cancerGenes_clinical", transc_allPrimaryGBM, diagnosis_allPrimaryGBM[, c("OS", "Recurrent_Num")])
      modules = list(c("cyan","salmon"))
      names(modules) = c("OS")
      outputWGCNAmodules(path="./results/2019.August_DASLvsClinical/CancerGenes/", title="WGCNA_allPrimaryGBM_CancerGenes_clinical", traitsModules = modules)
      
      ### TopGO with results
      {
        ## LASSO
        # For OS ==> inconclusive results
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = c("CR1", "CXorf22", "FAM23B", "TSHR"),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        
        ## WGCNA
        modulesOS = read.csv2(file="./results/2019.August_DASLvsClinical/CancerGenes/WGCNA_allPrimaryGBM_CancerGenes_clinical/outData/Modules_OS.csv", na.strings = "")
        names(modulesOS)
        dim(modulesOS)
        sapply(modulesOS[,c("cyan.Gene","salmon.Gene")],
               function(gene) nrow(modulesOS)-sum(is.na(gene)))
        
        # For cyan ==> inconclusive results
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$cyan.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        
        # For salmon ==> shows some interesting results with only a few genes per interesting term, on signalling pathways, stem cell division, GTPase activity, cell proliferation, glycolitic process
        unnumberedGenes = na.omit(gsub("\\.\\d+$", "", modulesOS$salmon.Gene))
        topGOresult = runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                               investigatedGenes = as.character(unnumberedGenes),
                               pvalues = NULL,
                               algorithm = "elim", statistic = "fisher")
        lapply(genesInTerm(topGOresult, c("GO:0051387","GO:0043409","GO:0070373","GO:0017145","GO:0090630","GO:0008285","GO:0006110")), 
               function(term) intersect(term,as.character(unnumberedGenes)))
        
      } # end of topGO
    }
    
    ## II.2.d) All Recurrent GBM samples ----
    {
      
      transc_allRecurrentGBM = transc[ match(samples_allRecurrentGBM, rownames(transc)),  ]
      diagnosis_allRecurrentGBM = diagnosis[ match(samples_allRecurrentGBM, rownames(diagnosis)),  ]
      
      ### PCA
      {
        acp = dudi.pca(transc_allRecurrentGBM, scannf = FALSE, nf = 10)
        pve <- 100*acp$eig/sum(acp$eig)
        cumsum(pve)
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.label(acp$li, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2)) 
        }
        par(mfrow=c(1,1))
        
        # Visualize distribution for Prim/Rec samples
        gcol = rainbow(length(levels(diagnosis_allRecurrentGBM$Recurrent)))
        par(mfrow=c(2,2))
        for(ax in 1:4){  
          s.class(dfxy = acp$li, fac=as.factor(diagnosis_allRecurrentGBM$Recurrent), col=gcol, xax=(2*(ax-1)+1), yax=(2*(ax-1)+2))
        }
        par(mfrow=c(1,1))
        
      }
      
      ### LASSO
      lasso_allRecurrentGBM_cancerGenes = runLASSO(transc_allRecurrentGBM, diagnosis_allRecurrentGBM[, c("OS")])
      print(lasso_allRecurrentGBM_cancerGenes)
      
      ### WGCNA
      runWGCNA("WGCNA_allRecurrentGBM_cancerGenes_clinical", transc_allRecurrentGBM, diagnosis_allRecurrentGBM[, c("OS", "Recurrent_Num")])
      
    }
    
    
  }# end of II.2 - DASL vs. Clinical data
  
  
  # ===== II.3 - DASL vs. AUCs =====
  {
    
    ## II.3.a) All 45 samples ----
    {
      
      transc_AllExposed = transc[ match(samples_drugExposedAll, rownames(transc)),  ]
      AUCs_AllExposed = AUCs[ match(samples_drugExposedAll, rownames(AUCs)),  ]
      cat(dim(transc_AllExposed), dim(AUCs_AllExposed))
      
      ### LASSO
      lasso_AllExposed_cancerGenes = runLASSO(transc_AllExposed, AUCs_AllExposed)
      dim(lasso_AllExposed_cancerGenes)
      write.csv2(lasso_AllExposed_cancerGenes, "./results/2019.November/AllExposed_AUCs/LASSO/results_AllExposed_AUCs.csv")
      
      ### WGCNA
      runWGCNA(path="./results/2019.November/AllExposed_AUCs/WGCNA/", title="WGCNA_AUCAllExposed_cancerGenes", 
               transc_AllExposed, AUCs_AllExposed)
      WGCNAcorrelation = read.table("./results/2019.November/AllExposed_AUCs/WGCNA/WGCNA_AUCAllExposed_cancerGenes/outData/WGCNAcorrelations_modulesDrugs.txt", h=T)
      modulesRun = list()
      for(drug in colnames(WGCNAcorrelation)[1:96]){
        mods = rownames(WGCNAcorrelation)[abs(WGCNAcorrelation[,drug]) > 0.4]
        if(length(mods) > 0){  modulesRun[[drug]] = gsub("ME", "", mods)  }
      }
      modulesRun
      
      outputWGCNAmodules("./results/2019.November/AllExposed_AUCs/WGCNA/", "WGCNA_AUCAllExposed_cancerGenes", modulesRun)
      
      
      ### TopGO with results
      {
        ## LASSO
        #lasso_AllExposed_cancerGenes = read.csv2("./results/2019.November/cancerGenesApproach/AllExposed_AUCs/LASSO/results_AllExposed_AUCs.csv", 
         #                                        row.names = 1, stringsAsFactors = F)
        for(drug in rownames(lasso_AllExposed_cancerGenes)){
          if(sum(!is.na(lasso_AllExposed_cancerGenes[drug,]))>10){
            lastGene_idx = sum(!is.na(lasso_AllExposed_cancerGenes[drug,]))
            cat(drug, lastGene_idx, sum(!is.na(lasso_AllExposed_cancerGenes[drug,])), '\n')
            #investigatedList = probesGenesMap_cancerGenes$probe.symbol[match(lasso_AllExposed_cancerGenes[drug,2:lastGene_idx], probesGenesMap_cancerGenes$probe.name)]
            runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                     investigatedGenes = ifelse(lasso_AllExposed_cancerGenes[drug,1] == "(Intercept)", 
                                                lasso_AllExposed_cancerGenes[drug,2:lastGene_idx],
                                                lasso_AllExposed_cancerGenes[drug,1:lastGene_idx]),
                     outPath = "./results/2019.November/cancerGenesApproach/AllExposed_AUCs/LASSO/", 
                     title = paste0("TopGO_LASSO_AllExposed_AUCs_",drug),
                     pvalues = NULL, topNods = 20, algorithm = "elim", statistic = "fisher")
          }
        }
        
        
        ## WGCNA
        for(drug in names(modulesRun)){
          
          modules = read.csv2(file=paste0("./results/2019.November/AllExposed_AUCs/WGCNA/WGCNA_AUCAllExposed_cancerGenes/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modulesRun[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleAllGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleAllGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleAllGenesSignif) >=0.3 & moduleAllGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/AllExposed_AUCs/WGCNA/WGCNA_AUCAllExposed_cancerGenes/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
      } # end of topGO
    }
    
    ## II.3.b) All Primary GBM samples that have IC50s available ----
    {
      transc_PrimaryExposed = transc[ match(samples_drugExposedPrimary, rownames(transc)),  ]
      AUCs_PrimaryExposed = AUCs[ match(samples_drugExposedPrimary, rownames(AUCs)),  ]
      cat(dim(transc_PrimaryExposed), dim(AUCs_PrimaryExposed))
      
      ### LASSO
      lasso_PrimaryExposed_cancerGenes = runLASSO(transc_PrimaryExposed, AUCs_PrimaryExposed)
      dim(lasso_PrimaryExposed_cancerGenes)
      write.csv2(lasso_PrimaryExposed_cancerGenes, "./results/2019.November/PrimaryExposed_AUCs/LASSO/results_PrimaryExposed_AUCs.csv")
      
      ### WGCNA
      runWGCNA(path="./results/2019.November/PrimaryExposed_AUCs/WGCNA/", title="WGCNA_AUCPrimaryExposed_cancerGenes", 
               transc_PrimaryExposed, AUCs_PrimaryExposed)
      WGCNAcorrelation = read.table("./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed_cancerGenes/outData/WGCNAcorrelations_modulesDrugs.txt", h=T)
      modulesRun = list()
      for(drug in colnames(WGCNAcorrelation)[1:96]){
        mods = rownames(WGCNAcorrelation)[abs(WGCNAcorrelation[,drug]) > 0.4]
        if(length(mods) > 0){  modulesRun[[drug]] = gsub("ME", "", mods)  }
      }
      modulesRun
      outputWGCNAmodules("./results/2019.November/PrimaryExposed_AUCs/WGCNA/", "WGCNA_AUCPrimaryExposed_cancerGenes", modulesRun)
      
      
      ### TopGO with results
      {
        ## LASSO
        lasso_PrimaryExposed_cancerGenes = read.csv2("./results/2019.November/cancerGenesApproach/PrimaryExposed_AUCs/LASSO/results_PrimaryExposed_AUCs.csv", 
                                                     row.names = 1, stringsAsFactors = F)
        for(drug in rownames(lasso_PrimaryExposed_cancerGenes)){
          if(sum(!is.na(lasso_PrimaryExposed_cancerGenes[drug,]))>10){
            lastGene_idx = sum(!is.na(lasso_PrimaryExposed_cancerGenes[drug,]))
            cat(drug, lastGene_idx, sum(!is.na(lasso_PrimaryExposed_cancerGenes[drug,])), '\n')
            runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                     investigatedGenes = ifelse(lasso_PrimaryExposed_cancerGenes[drug,1] == "(Intercept)", 
                                                lasso_PrimaryExposed_cancerGenes[drug,2:lastGene_idx],
                                                lasso_PrimaryExposed_cancerGenes[drug,1:lastGene_idx]),
                     outPath = "./results/2019.November/cancerGenesApproach/PrimaryExposed_AUCs/LASSO/", title = paste0("TopGO_LASSO_PrimaryExposed_AUCs_",drug),
                     pvalues = NULL, topNods = 20, algorithm = "elim", statistic = "fisher")
          }
        }

        ## WGCNA
        for(drug in names(modulesRun)){
          
          modules = read.csv2(file=paste0("./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed_cancerGenes/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modulesRun[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            modulePrimaryGenesSignif = as.numeric(as.character(moduleData[,2]))
            modulePrimaryGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(modulePrimaryGenesSignif) >=0.3 & modulePrimaryGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/PrimaryExposed_AUCs/WGCNA/WGCNA_AUCPrimaryExposed_cancerGenes/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
      } # end of topGO
    }
    
    ## II.3.c) All Recurrent GBM samples that have AUCs available ----
    {
      
      transc_RecurrentExposed = transc[ match(samples_drugExposedRecurrent, rownames(transc)),  ]
      AUCs_RecurrentExposed = AUCs[ match(samples_drugExposedRecurrent, rownames(AUCs)),  ]
      cat(dim(transc_RecurrentExposed), dim(AUCs_RecurrentExposed))
      
      ### LASSO
      lasso_RecurrentExposed_cancerGenes = runLASSO(transc_RecurrentExposed, AUCs_RecurrentExposed[1:96])
      dim(lasso_RecurrentExposed_cancerGenes)
      write.csv2(lasso_RecurrentExposed_cancerGenes, "./results/2019.November/RecurrentExposed_AUCs/LASSO/results_RecurrentExposed_AUCs.csv")
      
      ### WGCNA
      runWGCNA(path="./results/2019.November/RecurrentExposed_AUCs/WGCNA/", title="WGCNA_AUCRecurrentExposed_cancerGenes", 
               transc_RecurrentExposed, AUCs_RecurrentExposed)
      WGCNAcorrelation = read.table("./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed_cancerGenes/outData/WGCNAcorrelations_modulesDrugs.txt", h=T)
      modulesRun = list()
      for(drug in colnames(WGCNAcorrelation)[1:96]){
        mods = rownames(WGCNAcorrelation)[abs(WGCNAcorrelation[,drug]) > 0.7]
        if(length(mods) > 0){  modulesRun[[drug]] = gsub("ME", "", mods)  }
      }
      modulesRun
      outputWGCNAmodules("./results/2019.November/RecurrentExposed_AUCs/WGCNA/", "WGCNA_AUCRecurrentExposed_cancerGenes", modulesRun)
      
      
      ### TopGO with results
      {
        ## LASSO
        for(drug in rownames(lasso_RecurrentExposed_cancerGenes)){
          if(sum(!is.na(lasso_RecurrentExposed_cancerGenes[drug,]))>10){
            lastGene_idx = sum(!is.na(lasso_RecurrentExposed_cancerGenes[drug,]))
            #investigatedList = probesGenesMap_cancerGenes$probe.symbol[match(lasso_RecurrentExposed_cancerGenes[drug,2:lastGene_idx], probesGenesMap_cancerGenes$probe.name)]
            runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                     investigatedGenes = lasso_RecurrentExposed_cancerGenes[drug,2:lastGene_idx],
                     outPath = "./results/2019.November/RecurrentExposed_AUCs/LASSO/", title = paste0("TopGO_LASSO_RecurrentExposed_AUCs_",drug),
                     pvalues = NULL, algorithm = "elim", statistic = "fisher")
          }
        }
        
        ## WGCNA
        for(drug in names(modulesRun)){
          
          modules = read.csv2(file=paste0("./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed_cancerGenes/outData/Modules_", drug,".csv"), na.strings = "")
          
          for(module in modulesRun[[drug]]){
            cat(drug, module, '\n')
            moduleData = modules[, c(paste0(module, ".Gene"), paste0(module, ".Significance"),paste0(module, ".SignificancePvalue"))]
            moduleData = moduleData[!is.na(moduleData[,1]),]
            moduleRecurrentGenesSignif = as.numeric(as.character(moduleData[,2]))
            moduleRecurrentGenesSignifPval = as.numeric(as.character(moduleData[,3]))
            moduleGenes = moduleData[(abs(moduleRecurrentGenesSignif) >=0.3 & moduleRecurrentGenesSignifPval < 0.05), 1]
            runTopGO(geneUniverse = probesGenesMap_cancerGenes,
                     investigatedGenes = moduleGenes,
                     outPath = "./results/2019.November/RecurrentExposed_AUCs/WGCNA/WGCNA_AUCRecurrentExposed_cancerGenes/outData/", 
                     title = paste0("TopGO_",drug, "_", module),
                     pvalues = NULL, topNods = 20,
                     algorithm = "elim", statistic = "fisher")
          }
          
        }
        
      } # end of topGO
    }
    
    
  }# end of II.3 - DASL vs. IC50s
  
  
} # end of II - Cancer Genes - focused approach




# oooooooooooooooooooooooooooooooooooooooooooooooooo
# III - ALDH Genes and Retinoic Acid - focused approach ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # ===== III.1 - Load Data =====
  load("./tmp/samplesLists.RData")
  load("./tmp/transcFiltered_aldhGenes.RData")
  load("./tmp/ic50_datasets.RData")
  transc = transcFiltered_aldhGenes # for easier handling
  
  
  # ===== III.2 - With all 88 samples =====
  {
    
    ## Extract data
    transc_all88 = transc[ match(samples_all88, rownames(transc)),  ]
    diagnosis_all88 = diagnosis[ match(samples_all88, rownames(diagnosis)),  ]
    
    ALDH1A1av = (transc_all88[,1] + transc_all88[,2])/2 # Since the probes for ALDH1A1 (respecitvely, ALDH1A3) are relatively well correlated 
    ALDH1A3av = (transc_all88[,3] + transc_all88[,4])/2 #   (.78 between the probes for ALDH1A1 and .85 for ALDH1A3 probes, cf below), they are averaged
    SOX2.1 = transc_all88[,5]
    SOX2.2 = transc_all88[,6]
    OLIG2 = transc_all88[,7]
    GFAP = transc_all88[,8]
    survival = diagnosis_all88$OS
    names(survival) = rownames(diagnosis_all88)
    
    
    ## III.2.a) Simple vizualisation ----
    {
      
      boxplot(transc_all88)
      
      ## See if probes for a same gene have the same profile
      {
        probesDiff = cbind(transc_all88[,1] - transc_all88[,2], # ALDH1A1.1 - ALDH1A1.2
                           transc_all88[,3] - transc_all88[,4], # ALDH1A3.1 - ALDH1A3.2 
                           transc_all88[,5] - transc_all88[,6]) # SOX2.1 - SOX2.2
        colnames(probesDiff) = c("ALDH1A1", "ALDH1A3", "SOX2")
        head(probesDiff)
        
        boxplot(probesDiff)
        abline(h=0, lty=2)
        
        plot(probesDiff[,1], main="Difference of 2 probes for one gene", xlab="Samples", ylab="Difference", pch=20, col="blue")
        points(probesDiff[,2], pch=20, col="green")
        points(probesDiff[,3], pch=20, col="red")
        abline(h=0, lty=2)
        legend("topright", legend=colnames(probesDiff), col=c("blue","green","red"), pch=20)
      }
      
      ## Look at correlations
      {
        gradeColors = c("grey","green","yellow","orange", "red")
        
        # Correlation panel
        panel.cor <- function(x, y){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          corTest = cor.test(x, y, method="pearson")
          corVal <- round(corTest$estimate, digits=2)
          corPVal <- round(corTest$p.value, digits=3)
          txt <- paste0(corVal, "\n(", corPVal, ")")
          cex.cor <- 1
          text(0.5, 0.5, txt, cex = cex.cor, col=ifelse(corPVal<0.05, yes="red", no="black") )#* (0.0+abs(r)))
        }
        # Customize upper panel
        upper.panel<-function(x, y){
          points(x,y, pch = c(1,16)[diagnosis_all88$Recurrent], col = gradeColors[as.factor(diagnosis_all88$WHO_Grade)])
        }
        # Create the plots
        pairs(cbind(transc_all88, diagnosis_all88[, c("WHO_Grade", "OS")]),
              lower.panel = panel.cor,
              upper.panel = upper.panel)
      }
      
      
    }
    
    
    ## III.2.b) ANOVAs between ALDH Genes ----
    {
      
      aovALDH1A1 <- aov(ALDH1A1av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A1)
      
      aovALDH1A3 <- aov(ALDH1A3av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A3)
      
    }
    
    
    ## III.2.c) ANOVAs between ALDH Genes and Overall Survival ----
    {
      
      aovOS <- aov(survival ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS)
      
      aovOS1A1 <- aov(survival ~ ALDH1A1av*GFAP)
      summary(aovOS1A1)
      
      aovOS1A3 <- aov(survival ~ ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS1A3)
      
    }
    
    
  } # end of III.2 - With all 88 samples
  
  
  # ===== III.3 - With all GBM samples =====
  {
    
    ## Extract data
    transc_allGBM = transc[ match(samples_allGBM, rownames(transc)),  ]
    diagnosis_allGBM = diagnosis[ match(samples_allGBM, rownames(diagnosis)),  ]
    
    ALDH1A1av = (transc_allGBM[,1] + transc_allGBM[,2])/2 # Since the probes for ALDH1A1 (respecitvely, ALDH1A3) are relatively well correlated 
    ALDH1A3av = (transc_allGBM[,3] + transc_allGBM[,4])/2 #   (.78 between the probes for ALDH1A1 and .85 for ALDH1A3 probes, cf below), they are averaged
    SOX2.1 = transc_allGBM[,5]
    SOX2.2 = transc_allGBM[,6]
    OLIG2 = transc_allGBM[,7]
    GFAP = transc_allGBM[,8]
    survival = diagnosis_allGBM$OS
    names(survival) = rownames(diagnosis_allGBM)
    
    
    ## III.3.a) Simple vizualisation ----
    {
      
      boxplot(transc_allGBM)
      
      ## Look at correlations
      {
        
        # Correlation panel
        panel.cor <- function(x, y){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          corTest = cor.test(x, y, method="pearson")
          corVal <- round(corTest$estimate, digits=2)
          corPVal <- round(corTest$p.value, digits=3)
          txt <- paste0(corVal, "\n(", corPVal, ")")
          cex.cor <- 1
          text(0.5, 0.5, txt, cex = cex.cor, col=ifelse(corPVal<0.05, yes="red", no="black") )#* (0.0+abs(r)))
        }
        # Customize upper panel
        upper.panel<-function(x, y){
          points(x,y, pch = 20, col = c("blue","red")[as.factor(diagnosis_allGBM$Recurrent)])
        }
        # Create the plots
        pairs(cbind(transc_allGBM, diagnosis_allGBM[, "OS"]),
              lower.panel = panel.cor,
              upper.panel = upper.panel)
      }
      
      
    }
    
    
    ## III.3.b) ANOVAs between ALDH Genes ----
    {
      
      aovALDH1A1 <- aov(ALDH1A1av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A1)
      
      aovALDH1A3 <- aov(ALDH1A3av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A3)
      
    }
    
    
    ## III.3.c) ANOVAs between ALDH Genes and Overall Survival ----
    {
      
      aovOS <- aov(survival ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS)
      
      aovOS1A1 <- aov(survival ~ ALDH1A1av*GFAP)
      summary(aovOS1A1)
      
      aovOS1A3 <- aov(survival ~ ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS1A3)
      
    }
    
    
  } # end of III.3 - With all GBM samples
  
  
  # ===== III.4 - With all Primary GBM samples =====
  {
    
    ## Extract data
    transc_allPrimaryGBM = transc[ match(samples_allPrimaryGBM, rownames(transc)),  ]
    diagnosis_allPrimaryGBM = diagnosis[ match(samples_allPrimaryGBM, rownames(diagnosis)),  ]
    
    ALDH1A1av = (transc_allPrimaryGBM[,1] + transc_allPrimaryGBM[,2])/2 # Since the probes for ALDH1A1 (respecitvely, ALDH1A3) are relatively well correlated 
    ALDH1A3av = (transc_allPrimaryGBM[,3] + transc_allPrimaryGBM[,4])/2 #   (.78 between the probes for ALDH1A1 and .85 for ALDH1A3 probes, cf below), they are averaged
    SOX2.1 = transc_allPrimaryGBM[,5]
    SOX2.2 = transc_allPrimaryGBM[,6]
    OLIG2 = transc_allPrimaryGBM[,7]
    GFAP = transc_allPrimaryGBM[,8]
    survival = diagnosis_allPrimaryGBM$OS
    names(survival) = rownames(diagnosis_allPrimaryGBM)
    
    
    ## III.4.a) Simple vizualisation ----
    {
      
      boxplot(transc_allPrimaryGBM)
      
      ## Look at correlations
      {
        
        # Correlation panel
        panel.cor <- function(x, y){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          corTest = cor.test(x, y, method="pearson")
          corVal <- round(corTest$estimate, digits=2)
          corPVal <- round(corTest$p.value, digits=3)
          txt <- paste0(corVal, "\n(", corPVal, ")")
          cex.cor <- 1
          text(0.5, 0.5, txt, cex = cex.cor, col=ifelse(corPVal<0.05, yes="red", no="black") )#* (0.0+abs(r)))
        }
        # Customize upper panel
        upper.panel<-function(x, y){
          points(x,y, pch = 20, col = "black")
        }
        # Create the plots
        pairs(cbind(transc_allPrimaryGBM, diagnosis_allPrimaryGBM[, "OS"]),
              lower.panel = panel.cor,
              upper.panel = upper.panel)
      }
      
      
    }
    
    
    ## III.4.b) ANOVAs between ALDH Genes ----
    {
      
      aovALDH1A1 <- aov(ALDH1A1av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A1)
      
      aovALDH1A3 <- aov(ALDH1A3av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A3)
      
    }
    
    
    ## III.4.c) ANOVAs between ALDH Genes and Overall Survival ----
    {
      
      aovOS <- aov(survival ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS)
      
      aovOS1A1 <- aov(survival ~ ALDH1A1av*GFAP)
      summary(aovOS1A1)
      
      aovOS1A3 <- aov(survival ~ ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS1A3)
      
    }
    
    
  } # end of III.4 - With all GBM samples
  
  
  # ===== III.5 - With all samples that have Tretinoin IC50s available =====
  {
    
    ## Extract data
    tretinoinBounded = ic50s_bounded[!is.na(ic50s_bounded[,"Drug_random_70"]), "Drug_random_70"]
    samples_tretinoin = names(tretinoinBounded)
    tretinoinNormed = ic50s_normed[match(samples_tretinoin, rownames(ic50s_normed)), "Drug_random_70"]
    transc_tretinoin = transc[ match(samples_tretinoin, rownames(transc)),  ]
    diagnosis_tretinoin = diagnosis[ match(samples_tretinoin, rownames(diagnosis)),  ]
    
    ALDH1A1av = (transc_tretinoin[,1] + transc_tretinoin[,2])/2 # Since the probes for ALDH1A1 (respecitvely, ALDH1A3) are relatively well correlated 
    ALDH1A3av = (transc_tretinoin[,3] + transc_tretinoin[,4])/2 #   (.78 between the probes for ALDH1A1 and .85 for ALDH1A3 probes, cf below), they are averaged
    SOX2.1 = transc_tretinoin[,5]
    SOX2.2 = transc_tretinoin[,6]
    OLIG2 = transc_tretinoin[,7]
    GFAP = transc_tretinoin[,8]
    survival = diagnosis_tretinoin$OS
    names(survival) = rownames(diagnosis_tretinoin)
    
    
    ## III.5.a) Simple vizualisation ----
    {
      
      boxplot(transc_tretinoin)
      
      ## Look at correlations
      {
        
        # Correlation panel
        panel.cor <- function(x, y){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          corTest = cor.test(x, y, method="pearson")
          corVal <- round(corTest$estimate, digits=2)
          corPVal <- round(corTest$p.value, digits=3)
          txt <- paste0(corVal, "\n(", corPVal, ")")
          cex.cor <- 1
          text(0.5, 0.5, txt, cex = cex.cor, col=ifelse(corPVal<0.05, yes="red", no="black") )#* (0.0+abs(r)))
        }
        # Customize upper panel
        upper.panel<-function(x, y){
          points(x,y, pch = 20, col = "black")
        }
        # Create the plots
        transc_tretinoinIC50s = cbind(transc_tretinoin, tretinoinBounded, tretinoinNormed, survival)
        colnames(transc_tretinoinIC50s) = c(colnames(transc_tretinoin), "TretBounded","TretNormed", "OS")
        pairs(transc_tretinoinIC50s,
              lower.panel = panel.cor,
              upper.panel = upper.panel)
      }
      
      
    }
    
    
    ## III.5.b) ANOVAs between ALDH Genes and Overall Survival and drug response ----
    {
      
      aovOS <- aov(survival ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS)
      
      aovTretBounded <- aov(tretinoinBounded ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovTretBounded)
      
      aovTretNormed <- aov(tretinoinNormed ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovTretNormed)
      
    }
    
    
  } # end of III.5 - With all GBM samples
  
  
  # ===== III.6 - With Primary samples that have Tretinoin IC50s available =====
  {
    
    ## Extract data
    tretinoinBounded = ic50s_bounded[!is.na(ic50s_bounded[,"Drug_random_70"]), "Drug_random_70"]
    tretinoinBounded = tretinoinBounded[ names(tretinoinBounded) %in% samples_drugExposedPrimary  ]
    samples_tretinoin = names(tretinoinBounded)
    tretinoinNormed = ic50s_normed[match(samples_tretinoin, rownames(ic50s_normed)), "Drug_random_70"]
    transc_tretinoin = transc[ match(samples_tretinoin, rownames(transc)),  ]
    diagnosis_tretinoin = diagnosis[ match(samples_tretinoin, rownames(diagnosis)),  ]
    
    ALDH1A1av = (transc_tretinoin[,1] + transc_tretinoin[,2])/2 # Since the probes for ALDH1A1 (respecitvely, ALDH1A3) are relatively well correlated 
    ALDH1A3av = (transc_tretinoin[,3] + transc_tretinoin[,4])/2 #   (.78 between the probes for ALDH1A1 and .85 for ALDH1A3 probes, cf below), they are averaged
    SOX2.1 = transc_tretinoin[,5]
    SOX2.2 = transc_tretinoin[,6]
    OLIG2 = transc_tretinoin[,7]
    GFAP = transc_tretinoin[,8]
    survival = diagnosis_tretinoin$OS
    names(survival) = rownames(diagnosis_tretinoin)
    
    
    ## III.6.a) Simple vizualisation ----
    {
      
      boxplot(transc_tretinoin)
      
      ## Look at correlations
      {
        
        # Correlation panel
        panel.cor <- function(x, y){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          corTest = cor.test(x, y, method="pearson")
          corVal <- round(corTest$estimate, digits=2)
          corPVal <- round(corTest$p.value, digits=3)
          txt <- paste0(corVal, "\n(", corPVal, ")")
          cex.cor <- 1
          text(0.5, 0.5, txt, cex = cex.cor, col=ifelse(corPVal<0.05, yes="red", no="black") )#* (0.0+abs(r)))
        }
        # Customize upper panel
        upper.panel<-function(x, y){
          points(x,y, pch = 20, col = "black")
        }
        # Create the plots
        transc_tretinoinIC50s = cbind(transc_tretinoin, tretinoinBounded, tretinoinNormed, survival)
        colnames(transc_tretinoinIC50s) = c(colnames(transc_tretinoin), "TretBounded","TretNormed", "OS")
        pairs(transc_tretinoinIC50s,
              lower.panel = panel.cor,
              upper.panel = upper.panel)
      }
      
      
    }
    
    
    ## III.6.b) ANOVAs between ALDH Genes and Overall Survival and drug response ----
    {
      
      aovOS <- aov(survival ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS)
      
      aovTretBounded <- aov(tretinoinBounded ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovTretBounded)
      
      aovTretNormed <- aov(tretinoinNormed ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovTretNormed)
      
    }
    
    
  } # end of III.6 - With all GBM samples
  
  
  # ===== III.7 - With Recurrent samples that have Tretinoin IC50s available =====
  {
    
    ## Extract data
    tretinoinBounded = ic50s_bounded[!is.na(ic50s_bounded[,"Drug_random_70"]), "Drug_random_70"]
    tretinoinBounded = tretinoinBounded[ names(tretinoinBounded) %in% samples_drugExposedRecurrent  ]
    samples_tretinoin = names(tretinoinBounded)
    tretinoinNormed = ic50s_normed[match(samples_tretinoin, rownames(ic50s_normed)), "Drug_random_70"]
    transc_tretinoin = transc[ match(samples_tretinoin, rownames(transc)),  ]
    diagnosis_tretinoin = diagnosis[ match(samples_tretinoin, rownames(diagnosis)),  ]
    
    ALDH1A1av = (transc_tretinoin[,1] + transc_tretinoin[,2])/2 # Since the probes for ALDH1A1 (respecitvely, ALDH1A3) are relatively well correlated 
    ALDH1A3av = (transc_tretinoin[,3] + transc_tretinoin[,4])/2 #   (.78 between the probes for ALDH1A1 and .85 for ALDH1A3 probes, cf below), they are averaged
    SOX2.1 = transc_tretinoin[,5]
    SOX2.2 = transc_tretinoin[,6]
    OLIG2 = transc_tretinoin[,7]
    GFAP = transc_tretinoin[,8]
    survival = diagnosis_tretinoin$OS
    names(survival) = rownames(diagnosis_tretinoin)
    
    
    ## III.7.a) Simple vizualisation ----
    {
      
      boxplot(transc_tretinoin)
      
      ## Look at correlations
      {
        
        # Correlation panel
        panel.cor <- function(x, y){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          corTest = cor.test(x, y, method="pearson")
          corVal <- round(corTest$estimate, digits=2)
          corPVal <- round(corTest$p.value, digits=3)
          txt <- paste0(corVal, "\n(", corPVal, ")")
          cex.cor <- 1
          text(0.5, 0.5, txt, cex = cex.cor, col=ifelse(corPVal<0.05, yes="red", no="black") )#* (0.0+abs(r)))
        }
        # Customize upper panel
        upper.panel<-function(x, y){
          points(x,y, pch = 20, col = "black")
        }
        # Create the plots
        transc_tretinoinIC50s = cbind(transc_tretinoin, tretinoinBounded, tretinoinNormed, survival)
        colnames(transc_tretinoinIC50s) = c(colnames(transc_tretinoin), "TretBounded","TretNormed", "OS")
        pairs(transc_tretinoinIC50s,
              lower.panel = panel.cor,
              upper.panel = upper.panel)
      }
      
      
    }
    
    
    ## III.7.b) ANOVAs between ALDH Genes and Overall Survival and drug response ----
    {
      
      aovOS <- aov(survival ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovOS)
      
      aovTretBounded <- aov(tretinoinBounded ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovTretBounded)
      
      aovTretNormed <- aov(tretinoinNormed ~ ALDH1A1av*GFAP + ALDH1A3av*SOX2.1*SOX2.2*OLIG2)
      summary(aovTretNormed)
      
    }
    
    ## III.7.c) ANOVAs between ALDH Genes ----
    {
      
      aovALDH1A1 <- aov(ALDH1A1av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A1)
      
      aovALDH1A3 <- aov(ALDH1A3av ~ SOX2.1 + SOX2.2 + OLIG2 + GFAP +
                          SOX2.1:SOX2.2 + SOX2.1:OLIG2 + SOX2.1:GFAP +
                          SOX2.2:OLIG2 + SOX2.2:GFAP + OLIG2:GFAP)
      summary(aovALDH1A3)
      
    }
    
    
  } # end of III.7 - With all GBM samples
  
  
} # end of III - ALDH Genes and Retinoic Acid - focused approach




# oooooooooooooooooooooooooooooooooooooooooooooooooo
# IV - Analysis involving the complete dataset ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # ===== IV.1 - Load Data =====
  
  setwd("C:/Users/Romain Tching/Documents/GLIOTRAIN/GBMAnalysis/IoannisResistance/Data")
  load("./tmp/samplesLists.RData")
  load("./tmp/transc_preprocessed.RData")
  load("./tmp/ic50_datasets.RData")
  transc = transc_preprocessed # for easier handling
  
  
  # ===== IV.2 - Comparison of short-term vs long-term survivors =====

  ## Using all GBM samples ----
  {
    transc_allGBM = transc[ match(samples_allGBM, rownames(transc)),  ]
    
    ## Define groups
    transc_STSurvivors = transc_preprocessed[intersect(samples_allSTsurvivors, rownames(transc_allGBM)), ]
    transc_LTSurvivors = transc_preprocessed[intersect(samples_allLTsurvivors, rownames(transc_allGBM)), ]
    Nprobes = ncol(transc_allGBM)
    
    ## Run the t-test
    ttestResults = lapply(  1:Nprobes, function(probe) t.test( transc_STSurvivors[, probe], transc_LTSurvivors[, probe] )  )
    names(ttestResults) = colnames(transc_allGBM)
    ttestResults_pvals = sapply(ttestResults, function(probeTest) probeTest$p.value)
    ttestResults_Signif = ttestResults_pvals[ttestResults_pvals<0.05]
    length(ttestResults_Signif)
    ttestResults_Signif
    
    ## Bonferroni multi-test correction
    ttestResults_Bonferroni = (p.adjust(ttestResults_pvals, "bonferroni"))
    ttestResults_BonferroniSignif = ttestResults_Bonferroni[ttestResults_Bonferroni<0.05]
    length(ttestResults_BonferroniSignif)
    ttestResults_BonferroniSignif
    
    ## FDR multi-test correction
    ttestResults_FDR = (p.adjust(ttestResults_pvals, "fdr"))
    ttestResults_FDRSignif = ttestResults_FDR[ttestResults_FDR<0.05]
    length(ttestResults_FDRSignif)
    ttestResults_FDRSignif
    
  }
  
  ## Using all Primary GBM samples ----
  {
    transc_allPrimaryGBM = transc[ match(samples_allPrimaryGBM, rownames(transc)),  ]
    
    ## Define groups
    transc_STSurvivors = transc_preprocessed[intersect(samples_allSTsurvivors, rownames(transc_allPrimaryGBM)), ]
    transc_LTSurvivors = transc_preprocessed[intersect(samples_allLTsurvivors, rownames(transc_allPrimaryGBM)), ]
    Nprobes = ncol(transc_allPrimaryGBM)
    
    ## Run the t-test
    ttestResults = lapply(  1:Nprobes, function(probe) t.test( transc_STSurvivors[, probe], transc_LTSurvivors[, probe] )  )
    names(ttestResults) = colnames(transc_allPrimaryGBM)
    ttestResults_pvals = sapply(ttestResults, function(probeTest) probeTest$p.value)
    ttestResults_Signif = ttestResults_pvals[ttestResults_pvals<0.05]
    length(ttestResults_Signif)
    ttestResults_Signif
    
    ## Bonferroni multi-test correction
    ttestResults_Bonferroni = (p.adjust(ttestResults_pvals, "bonferroni"))
    ttestResults_BonferroniSignif = ttestResults_Bonferroni[ttestResults_Bonferroni<0.05]
    length(ttestResults_BonferroniSignif)
    ttestResults_BonferroniSignif
    
    ## FDR multi-test correction
    ttestResults_FDR = (p.adjust(ttestResults_pvals, "fdr"))
    ttestResults_FDRSignif = ttestResults_FDR[ttestResults_FDR<0.05]
    length(ttestResults_FDRSignif)
    ttestResults_FDRSignif
    
  }

  ## Looking at SPRY4 and PRKAR1B distributions ----
  {
    diagnosis_OS = diagnosis[!is.na(diagnosis$OS),]
    dim(diagnosis_OS)
    diagnosis_OSGBM = diagnosis_OS[diagnosis_OS$Pathological_diagnosis=="GBM", ]
    dim(diagnosis_OSGBM)
    diagnosis_OSPGBM = diagnosis_OSGBM[diagnosis_OSGBM$Recurrent=="No", ]
    dim(diagnosis_OSPGBM)
    
    transc_OS = transc[match(rownames(diagnosis_OS),rownames(transc)),]
    transc_OSGBM = transc_OS[match(rownames(diagnosis_OSGBM),rownames(transc_OS)),]
    transc_OSPGBM = transc_OSGBM[match(rownames(diagnosis_OSPGBM),rownames(transc_OSGBM)),]
    
    ## PRKAR1B
    # With all samples having data
    cor.test(diagnosis_OS$OS, transc_OS[,"PRKAR1B.1"])
    plot(diagnosis_OS$OS, transc_OS[,"PRKAR1B.1"], xlab = "Survival (months)", ylab = "PRKAR1B.1 expression",
         pch = ifelse(diagnosis_OS$Recurrent=="Yes", 1, 20), 
         col = ifelse(diagnosis_OS$Pathological_diagnosis=="GBM","red","black"))
    abline(v=c(9, 36), col="red", lty=2)
    
    # Only GBM Samples
    cor.test(diagnosis_OSGBM$OS, transc_OSGBM[,"PRKAR1B.1"])
    plot(diagnosis_OSGBM$OS, transc_OSGBM[,"PRKAR1B.1"], xlab = "Survival (months)", ylab = "PRKAR1B.1 expression",
         pch = ifelse(diagnosis_OSGBM$Recurrent=="Yes", 1, 20))
    abline(v=c(9, 36), col="red", lty=2)
    
    # Only Primary GBM Samples
    cor.test(diagnosis_OSPGBM$OS, transc_OSPGBM[,"PRKAR1B.1"])
    plot(diagnosis_OSPGBM$OS, transc_OSPGBM[,"PRKAR1B.1"], xlab = "Survival (months)", ylab = "PRKAR1B.1 expression",
         pch = ifelse(diagnosis_OSPGBM$Recurrent=="Yes", 1, 20))
    abline(v=c(9, 36), col="red", lty=2)
    
    
    ## SPRY4
    # With all samples having data
    cor.test(diagnosis_OS$OS, transc_OS[,"SPRY4.3"])
    plot(diagnosis_OS$OS, transc_OS[,"SPRY4.3"], xlab = "Survival (months)", ylab = "SPRY4.3 expression",
         pch = ifelse(diagnosis_OS$Recurrent=="Yes", 1, 20), 
         col = ifelse(diagnosis_OS$Pathological_diagnosis=="GBM","red","black"))
    abline(v=c(9, 36), col="red", lty=2)
    
    # Only GBM Samples
    glop = cor.test(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.3"])
    corVal <- round(glop$estimate, digits=2) # extract the value of the correlation
    corPVal <- round(glop$p.value, digits=3) # extract the pvalue from the test
    plot(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.3"], xlab = "Survival (months)", ylab = "SPRY4 expression",
         pch = ifelse(diagnosis_OSGBM$Recurrent=="Yes", 1, 20))
    legend("topright", legend=c("Primary","Recurrent"), pch=c(20,1))
    abline(a = transc_OSGBM[1,"SPRY4.3"] - corVal*diagnosis_OSGBM$OS[1], # intercept of the line; just make sure that OS[1] and olig2[1] are not NA
           b = corVal) # slope of the line
    text(70, 10, # coordinates for the text
         paste0("cor = ", corVal, "\npval = ", corPVal),
         col="black", cex=0.75)
    
    abline(v=c(9, 36), col="red", lty=2)
    
    # Only Primary GBM Samples
    cor.test(diagnosis_OSPGBM$OS, transc_OSPGBM[,"SPRY4.3"])
    plot(diagnosis_OSPGBM$OS, transc_OSPGBM[,"SPRY4.3"], xlab = "Survival (months)", ylab = "SPRY4.3 expression",
         pch = ifelse(diagnosis_OSPGBM$Recurrent=="Yes", 1, 20))
    abline(v=c(9, 36), col="red", lty=2)
    
    # All SPRY4 probes with GBM samples
    plot(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.3"], xlab = "Survival (months)", ylab = "SPRY4 expression", ylim=c(7.8,15.5),
         pch = ifelse(diagnosis_OSGBM$Recurrent=="Yes", 1, 20))
    abline(v=c(9, 36), col="red", lty=2)
    points(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.1"], col = "blue",
         pch = ifelse(diagnosis_OSGBM$Recurrent=="Yes", 1, 20))
    points(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.2"], col = "green",
           pch = ifelse(diagnosis_OSGBM$Recurrent=="Yes", 1, 20))
    points(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.4"], col = "red",
           pch = ifelse(diagnosis_OSGBM$Recurrent=="Yes", 1, 20))
    cor.test(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.1"])
    cor.test(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.2"])
    cor.test(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY4.4"])
    
    ## SPRY1,2,3
    plot(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY1.1"], xlab = "Survival (months)", ylab = "SPRY expression",
         pch = 20, col = "blue")
    abline(v=c(9, 36), col="red", lty=2)
    points(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY1.2"], xlab = "Survival (months)", ylab = "SPRY expression",
           pch = 20, col = "cyan")
    points(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY1.3"], xlab = "Survival (months)", ylab = "SPRY expression",
           pch = 20, col = "purple")
    plot(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY2.1"], xlab = "Survival (months)", ylab = "SPRY expression",
         pch = 20, col = "green")
    plot(diagnosis_OSGBM$OS, transc_OSGBM[,"SPRY3.1"], xlab = "Survival (months)", ylab = "SPRY expression",
           pch = 20, col = "red")
  }


}



# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# V - Investigation of the Multivariate Analyses Results ----
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
{
  
  # ===== V.1 - Unbiased approach =====
  load("./tmp/samplesLists.RData")
  load("./tmp/transcFiltered_unbiasedApproach.RData")

  ## V.1.a) LASSO Results ----


} # end of Section V-





# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# VI - Pulp analysis ----
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

load("./tmp/transcFiltered_unbiasedApproach.RData")
load("./tmp/transcFiltered_cancerGenes.RData")
load("./tmp/transc_preprocessed.RData")
colnames(transc_preprocessed)[grep("ABCB1\\.", colnames(transc_preprocessed))]
colnames(transcFiltered_unbiasedApproach)[grep("ABCB1", colnames(transcFiltered_unbiasedApproach))]
colnames(transcFiltered_cancerGenes)[grep("ABCB1\\.", colnames(transcFiltered_cancerGenes))]
dim(transcFiltered_cancerGenes)
ABCB1.1 = transcFiltered_cancerGenes[,"ABCB1.1"]
ABCB1.2 = transcFiltered_cancerGenes[,"ABCB1.2"]
updatedIC50s = read.csv("./tmp/ic50updated_raw.csv", row.names = 1)
updatedIC50s[1:5,1:5]

drugsTested = c("Drug_random_89", "Drug_random_84", "Drug_random_12", 
                "Drug_random_63", "Drug_random_66", "Drug_random_79")

for(drug in drugsTested){
  
  extracted = updatedIC50s[,drug]
  extrNum = as.numeric(as.character(extracted))
  names(extrNum) = rownames(updatedIC50s)
  drug_ic50s = extrNum[!is.na(extrNum)]
  
  drug_ABCB1.1 = ABCB1.1[match(names(drug_ic50s), names(ABCB1.1))]
  drug_ABCB1.2 = ABCB1.2[match(names(drug_ic50s), names(ABCB1.2))]
  
  ## Run correlation test
  corTest1 = cor.test(drug_ic50s, drug_ABCB1.1, method="spearman") # Correlation test
  corVal1 <- round(corTest1$estimate, digits=2) # extract the value of the correlation
  corPVal1 <- round(corTest1$p.value, digits=3) # extract the pvalue from the test
  
  corTest2 = cor.test(drug_ic50s, drug_ABCB1.2, method="spearman") # Correlation test
  corVal2 <- round(corTest2$estimate, digits=2) # extract the value of the correlation
  corPVal2 <- round(corTest2$p.value, digits=3) # extract the pvalue from the test
  
  ## Plot corresonding graph
  png(filename = paste0("./results/SpearmanCorrelation_ABCB1_",drug,".png"))
  par(mfrow=c(1,2))
  plot(drug_ic50s, drug_ABCB1.1, pch=20, col="black",
       main = paste0(drug," ~ ABCB1.1\n ","cor = ", corVal1, " pval = ", corPVal1),
       xlab="IC50", ylab="ABCB1.1 expression")
  abline(a = drug_ABCB1.1[1] - corVal1*drug_ic50s[1], # intercept of the line; just make sure that OS[1] and olig2[1] are not NA
         b = corVal1) # slope of the line
  
  plot(drug_ic50s, drug_ABCB1.2, pch=20, col="black",
       main = paste0(drug," ~ ABCB1.2\n ","cor = ", corVal2, " pval = ", corPVal2),
       xlab="IC50", ylab="ABCB1.2 expression")
  abline(a = drug_ABCB1.2[1] - corVal2*drug_ic50s[1], # intercept of the line; just make sure that OS[1] and olig2[1] are not NA
         b = corVal2) # slope of the line
  par(mfrow=c(1,1))
  dev.off()
}

d19 = updatedIC50s$Drug_random_19
names(d19) = rownames(updatedIC50s)
d19_sensitive = c()
d19_resistant = c()
for(s in 1:length(d19)){
  if(d19[s] == ">160"){  d19_resistant = c(d19_resistant, names(d19)[s])  }
  else if(!is.na(as.numeric(as.character(d19[s])))){
    val = as.numeric(as.character(d19[s]))
    if(val<160){  d19_sensitive = c(d19_sensitive, names(d19)[s])  }
    else{  d19_resistant = c(d19_resistant, names(d19)[s])  }
  }
}
cat("Sensitive: ", length(d19_sensitive), " ; Resistant: ", length(d19_resistant), "\n")

ABCB1.1_res = ABCB1.1[names(ABCB1.1) %in% d19_resistant]
ABCB1.1_sen = ABCB1.1[names(ABCB1.1) %in% d19_sensitive]
ABCB1.2_res = ABCB1.2[names(ABCB1.2) %in% d19_resistant]
ABCB1.2_sen = ABCB1.2[names(ABCB1.2) %in% d19_sensitive]
d19_groups = list(ABCB1.1_res,ABCB1.1_sen,ABCB1.2_res,ABCB1.2_sen)
names(d19_groups) = c("ABCB1.1_res","ABCB1.1_sen","ABCB1.2_res","ABCB1.2_sen")

ttest1 = t.test(ABCB1.1_res, ABCB1.1_sen)
ttest2 = t.test(ABCB1.2_res, ABCB1.2_sen)
pval1 = round(ttest1$p.value, 3)
pval2 = round(ttest2$p.value, 3)

boxplot(d19_groups, ylab = "ABCB1 expression", xlab="",
        main=paste0("ABCB1 ~ Drug_19 responsiveness\npval_1 = ", pval1, ", pval_2 = ", pval2))


ABCB1.1_mean = mean(ABCB1.1)
ABCB1.1_sdev = sd(ABCB1.1)
ABCB1.2_mean = mean(ABCB1.2)
ABCB1.2_sdev = sd(ABCB1.2)

ABCB1.1_extrems = ABCB1.1[abs(ABCB1.1 - ABCB1.1_mean) > (2*ABCB1.1_sdev)]
ABCB1.2_extrems = ABCB1.2[abs(ABCB1.2 - ABCB1.2_mean) > (2*ABCB1.2_sdev)]
