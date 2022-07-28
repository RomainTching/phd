# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#                   GLIOTRAIN Data Preprocessing and Extraction
#
# Process the data from GLIOTRAIN to have it ready for tranSMART loading. That includes:
#   * For Genomics (low coverage WGS) data:
#     - Identifying the overlap of the identified regions between samples
#   * For RNA Seq data:
#     - 
#   * For RPPA data
#     -  
#   * For Clinical data
#     -  
#   * For Methylation data
#     - 
#
# Started on the 27th of August, 2019
# 
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

### Environment set-up ----
setwd("/path/to/working_directory/")

library(readxl)
library(biomaRt)
library(org.Hs.eg.db)
library(ChromHeatMap)
library(ggplot2)
library(missMDA)
library(ade4)
library(Rtsne)
library(umap)
library(glmpca)
library(DESeq2)
#library(FactoMineR)

# gliotrain_id = gsub("\\.", "-", wgs_samplesMap$`Sample ID`)
# wgs_id = wgs_samplesMap$DILA_ID_DNA
# rnaSeq_id = gsub("s*(FLO\\d+)_\\w+(_\\w+)", "r\\1\\2", wgs_samplesMap$DILA_ID_DNA)
# sampleIDs_mapping = as.data.frame(cbind(gliotrain_id, wgs_id, rnaSeq_id))
# write.csv2(sampleIDs_mapping, row.names=F, file = "./tmp/samplesIDs_mapping.csv")


# oooooooooooooooooooooooooooooooooooooooooooooooooo
# I - Clinical Data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

clinical = read.csv("./tmp/clinical_data.csv")
sample_id = gsub("\\.", "\\-", clinical[,1])
patient_id = gsub("\\.", "\\-", clinical[,1])
clinical = cbind(sample_id, patient_id, clinical)
clinical = clinical[, -3]
clinical[,"sex"] = ifelse(clinical[,"sex"] == 'F', "Female", "Male")
clinical[1:5, 1:5]

write.table(clinical, file = "./readyForUpload/parentalTumors/transmart/data_clinical.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = '')


# oooooooooooooooooooooooooooooooooooooooooooooooooo
# II - WGS Data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

wgs_raw = read.table("./rawData/seg_Baf_logR_sorted_WGS.txt", h=T, sep='\t', stringsAsFactors = T)
wgs_samplesMap = read.csv("./tmp/samplesLabels_WGS.csv", h=T)
gistic_amps = (read_excel("./tmp/gistic_ampDels.xlsx", na="NA", sheet = 1, col_names = T))
gistic_dels = (read_excel("./tmp/gistic_ampDels.xlsx", na="NA", sheet = 2, col_names = T))

# ===== II.1 Get a sense of the data =====

### II.1.a) General distribution ----
head(wgs_raw)
head(wgs_samplesMap)
length(levels(wgs_raw$Sample)) # number of samples
length(levels(wgs_raw$Chromosome)) # number of chromosomes appearing
levels(wgs_raw$Chromosome)  # list of chromosomes; ==> no chromosome Y?
table(wgs_raw$Chromosome)[c(paste("chr",1:22, sep=""),"chrX")] # distribution of the number of regions referenced per chromosome
barplot(table(wgs_raw$Chromosome)[c(paste("chr",1:22, sep=""),"chrX")])

hist(wgs_raw$log_R, breaks= seq(min(wgs_raw$log_R)-0.1, max(wgs_raw$log_R)+0.1, 0.1))
hist(wgs_raw$log_R, breaks= seq(min(wgs_raw$log_R)-0.1, max(wgs_raw$log_R)+0.1, 0.1), ylim=c(0,650))


### II.1.b) Looking at the overlap between samples ----
test = wgs_raw #[(wgs_raw$log_R < -0.5 | wgs_raw$log_R > 0.5) , ]
lowestVal  = min(test$log_R)
highestVal = max(test$log_R)

for(chr in levels(test$Chromosome)){
  
  dataChr = test[test$Chromosome==chr,]
  samplesInt = match(dataChr$Sample,levels(test$Sample))
  dataChr = as.data.frame(cbind(dataChr, samplesInt))
  
  #png(filename = paste0("~/Events/2019.10_GT.AnnualMeeting.SystemsBiologyWorkshop_Stuttgart/Presentation/WGSraw_", chr,"contigs.png"))
  print(
    ggplot(data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="darkgrey")) + # create empty plot
      labs(title = chr, x = "Position", y = "Samples", colour = "log_R") +
      xlim(  (min(dataChr$Start) - 1000), (max(dataChr$End) + 1000)  ) + ylim( c( 0, length(levels(dataChr$Sample))+1) )  +
      scale_colour_gradient2(limits=c(lowestVal,highestVal), low="blue",mid="white",high="red") +
      geom_segment(data=dataChr, aes(x=dataChr$Start, xend=dataChr$End, 
                                     y=dataChr$samplesInt, yend=dataChr$samplesInt, 
                                     color=dataChr$log_R)) 
  )
  #dev.off()
  
}


# ===== II.2 Fix segments size =====

### II.2.a) By breaking down segments into common start/end position ----

remodelled_chr = list()

for(chr in levels(wgs_raw$Chromosome)){
  
  dataChr = wgs_raw[ wgs_raw$Chromosome==chr, ]
  positions = sort( unique( c( dataChr$Start, dataChr$End ) ) )
  remodelled = matrix(data=NA, nrow=length(levels(wgs_raw$Sample)), ncol=(length(positions)-1))
  rownames(remodelled) = levels(wgs_raw$Sample)
  colnames(remodelled) = sapply(  1:(length(positions)-1), function(i)  paste(chr, positions[i], positions[i+1], sep='_')  )
  
  for(samp in 1:nrow(remodelled)){
    
    dataChrSamp = dataChr[  dataChr$Sample == rownames(remodelled)[samp],  ]
    
    for(p in 1:ncol(remodelled)){
      if(any( dataChrSamp$Start<=positions[p] & dataChrSamp$End>=positions[p+1] )){
        remodelled[samp, p] = dataChrSamp$log_R[  which( dataChrSamp$Start<=positions[p] & dataChrSamp$End>=positions[p+1] ) ]
      }
    }
  }
  
  remodelled_chr[[chr]] = remodelled
  
}

wgs_matrix = NULL
for(chr in 1:(length(remodelled_chr))){  wgs_matrix = cbind(wgs_matrix, remodelled_chr[[chr]])  }
dim(wgs_matrix)
sum(is.na(wgs_matrix))
colnames(wgs_matrix[, which(apply(wgs_matrix,2, function(seg) sum(is.na(seg))==nrow(wgs_matrix)))])
wgs_matrixTransmart = wgs_matrix[ , !apply(wgs_matrix,2, function(seg) sum(is.na(seg))==nrow(wgs_matrix))  ]
dim(wgs_matrixTransmart)
sum(is.na(wgs_matrixTransmart))
wgs_matrixTransmart[1:5,1:5]

## Characterisation of the produced matrix
sapply(remodelled_chr, ncol)[c(paste("chr",1:22, sep=""),"chrX")]

## Map genes to chromosomal regions 
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genesMap <- getBM( attributes = c("ensembl_gene_id","hgnc_symbol", "external_gene_name",
                                  "chromosome_name", "start_position", "end_position", "band"),
                   mart = ensembl )
genesMap = genesMap[genesMap$chromosome_name %in% c('X', as.character(1:22)),]
results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                "start_position","end_position"),
                 mart=ensembl)

region2genes_map = data.frame(region = colnames(wgs_matrixTransmart), stringsAsFactors = F)
for(r in 1:nrow(region2genes_map)){
  pos = strsplit(region2genes_map$region[r], '_')[[1]]
  region2genes_map$chr[r] = gsub('chr', '', pos[1])
  region2genes_map$start[r] = as.numeric(pos[2])
  region2genes_map$end[r] = as.numeric(pos[3])
}

for(i in 1:nrow(region2genes_map)){
  bmChr = genesMap[genesMap$chromosome_name == region2genes_map$chr[i], ]
  geneInDataChr = c()
  for(j in 1:nrow(bmChr)){
    if(region2genes_map$start[i] < bmChr$end_position[j] & region2genes_map$end[i] > bmChr$start_position[j]){
      geneInDataChr = c(geneInDataChr, bmChr$external_gene_name[j])
    }
  }
  if(length(geneInDataChr) > 0){  region2genes_map$genes[i] = paste(geneInDataChr, collapse = ', ')  }
}

head(region2genes_map)
write.table(x = region2genes_map, file = "./tmp/logR2genes_map.tsv",  sep = "\t", row.names = F)

### II.2.b) Based on GISTIC-identified focal events ----

gistic_ampDels = as.data.frame(rbind((gistic_amps[,2:ncol(gistic_amps)]), (-1 * (gistic_dels[,2:ncol(gistic_dels)]))))
rownames_gistic = gsub("\\(.*\\)", "", c(gistic_amps[,1][[1]], gistic_dels[,1][[1]]))
rownames_gistic = gsub("[:-]", "_", rownames_gistic)
rownames(gistic_ampDels) = rownames_gistic
gistic_ampDels[1:5,1:5]

{
# 
# amplifiedBands = c("1q32.1", "1q44", "3q26.33", "4q12", "5p15.33", "7p14.3", "7p11.2", "7p21.2", "7q31.2", 
#                    "8p23.3", "8q24.21", "11q25", "12q13.3", "12q15")
# deletedBands = c("1p36.23", "1p36.11", "1p32.3", "1q32.3", "2p16.1", "2q33.1", "3q29", "4p14", "4q35.2", "5p15.33", "5p15.1", 
#                  "6q12", "6q23.3", "6q26", "9p23", "9p21.3", "9p21.1", "10q21.3", "10q22.1", "10q23.31", "10q26.2", "11q13.3", 
#                  "13q21.2", "13q34", "14q13.1", "15q21.2", "16p13.3", "16p12.2", "16q12.1", "16q22.1", "17q11.2", "17q21.1", 
#                  "18q11.1", "19q13.33")
# eventsBands = c(amplifiedBands, deletedBands)
# 
# eventsChrStr = gsub("[p,q].*$", "", eventsBands )
# eventsArms = gsub("(\\d+)([p,q])(.*$)", "\\2", eventsBands )
# eventsCytobands = gsub("(\\d+)([p,q])(.*$)", "\\3", eventsBands )
# regionsOfInterest = data.frame(eventsChrStr, eventsArms, eventsCytobands, start=NA, end=NA, amplified=NA)
# colnames(regionsOfInterest) = c("Chromosome", "Arm", "Band", "Start", "End", "Amplified")
# 
# data(cytobands)
# cytobandsAll = cytobands$`Homo sapiens`
# 
# 
# # ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "89")
# # genesMap <- getBM( attributes = c("ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_name", 
# #                                   "chromosome_name", "start_position", "end_position", "band"),
# #                    mart = ensembl )
# # regionsMap = genesMap[genesMap$chromosome_name==regionsOfInterest$Chromosome & genesMap$band==regionsOfInterest$Band,]
# 
# 
# 
# 
# for(band in 1:nrow(regionsOfInterest)){
# 
#   regionsOfInterest$Start[band] = cytobandsAll$start[  cytobandsAll$chr==paste0("chr",regionsOfInterest$Chromosome[band]) & 
#                                                          cytobandsAll$arm==regionsOfInterest$Arm[band] &
#                                                          cytobandsAll$band==regionsOfInterest$Band[band]  ]
#   regionsOfInterest$End[band] = cytobandsAll$end[  cytobandsAll$chr==paste0("chr",regionsOfInterest$Chromosome[band]) & 
#                                                      cytobandsAll$arm==regionsOfInterest$Arm[band] &
#                                                      cytobandsAll$band==regionsOfInterest$Band[band]  ]
#   
#   regionsOfInterest$Amplified[band] = ifelse(paste0(regionsOfInterest$Chromosome[band], regionsOfInterest$Arm[band], regionsOfInterest$Band[band]) %in% amplifiedBands, "red", "blue")
#   
# }
# 
# for(chr in levels(test$Chromosome)){
#   
#   dataChr = test[test$Chromosome==chr,]
#   samplesInt = match(dataChr$Sample,levels(test$Sample))
#   dataChr = as.data.frame(cbind(dataChr, samplesInt))
#   
#   showBands = regionsOfInterest[paste0("chr",regionsOfInterest$Chromosome)==chr,]
#   
#   print(
#     ggplot(data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="darkgrey")) + # create empty plot
#       labs(title = chr, x = "Position", y = "Samples", colour = "log_R") +
#       xlim(  (min(dataChr$Start) - 1000), (max(dataChr$End) + 1000)  ) + ylim( c( 0, length(levels(dataChr$Sample))+1) )  +
#       scale_colour_gradient2(limits=c(lowestVal,highestVal), low="blue",mid="white",high="red") +
#       geom_segment(data=dataChr, aes(x=dataChr$Start, xend=dataChr$End, 
#                                      y=dataChr$samplesInt, yend=dataChr$samplesInt, 
#                                      color=dataChr$log_R)) +
#       geom_vline(xintercept=showBands$Start, color=showBands$Amplified) +
#       geom_vline(xintercept=showBands$End, linetype="dashed", color=showBands$Amplified)
#   )
#   
#   # showBands = regionsOfInterest[regionsOfInterest$Chromosome==chr,]
#   # for(band in 1:nrow(showBands)){
#   #   abline(v=showBands$Start[band], lty=1, col=c("red","blue")[showBands$Amplified[band]], lwd=20000000)
#   #   abline(v=showBands$End[band], lty=2, col=c("red","blue")[showBands$Amplified[band]])
#   # }
#   #abline(v=showBands$Start, lty=1, col=c("red","blue")[showBands$Amplified])
#   #abline(v=showBands$End, lty=2, col=c("red","blue")[showBands$Amplified])
#   # abline(v=32566789, lty=1, col="blue")
#   # abline(v=33039907, lty=2, col="blue")
#   # geom_vline(xintercept=32566789, color="blue")
# }
# 
# 
# ### II.2.b) By defining fixed-size segments through merging and averaging 
# 
# for(chr in levels(test$Chromosome)){
#   dataChr = test[test$Chromosome==chr,]
#   print(chr)
#   print("Start")
#   tabStart = table(dataChr$Start)[table(dataChr$Start)>50]
#   print(tabStart)
#   print("End")
#   tabEnd = table(dataChr$End)[table(dataChr$End)>50]
#   print(tabEnd)
# }
# emptyRegions = data.frame(chromosome = c("chr1", "chr1", "chr2",   "chr3",   "chr4",   "chr5",   "chr5",     "chr6",   "chr7", "chr8",  "chr9",   "chr10",   "chr11",   "chr12", "chr16", "chr17",   "chr18",  "chr19", "chr20", "chrX", "chrX"),
#                           start = c(12825000, 121325000, 89025000, 90175000, 49025000, 45875000, 68775000, 57675000, 57625000, 43275000, 38775000, 38425000, 50175000, 34275000, 35075000, 22125000, 14775000, 24325000, 26125000, 58275000, 154913804),
#                           end   = c(13825000, 145425000, 95625000, 93575000, 52725000, 49625000, 70775000, 62025000, 63425000, 47625000, 71025000, 42875000, 55075000, 38475000, 46575000, 25325000, 18575000, 28375000, 29875000, 62075000, 154925000)
#                           )
# 
# newSegments = data.frame(Sample=character(), Chromosome=character(), Start=integer(), End=integer(), log_R=double())
# for(chr in levels(test$Chromosome)){
#   dataChr = test[test$Chromosome==chr, ]
#   #print(table(dataChr$Sample))
#   missingRegion = emptyRegions[emptyRegions$chromosome==chr, ]
#   for(samp in levels(dataChr$Sample)){
#     dataChrSamp = dataChr[dataChr$Sample==samp, ]
#     dataChrSamp = dataChrSamp[order(dataChrSamp$Start),] # order to have contiguous segments next to each other
#     #cat(chr, samp,"\n")
#     #print(dataChrSamp)
#     nextSegment = 2
#     newSegment = dataChrSamp[1, ]
#       #c(dataChrSamp$Sample[1], dataChrSamp$Chromosome[1], dataChrSamp$Start[1], dataChrSamp$End[1], dataChrSamp$log_R[1])
#     #names(newSegment) = colnames(dataChrSamp)
#     while(nextSegment <= nrow(dataChrSamp)+1){
#       
#       if(nextSegment == nrow(dataChrSamp)+1){ 
#         newSegments = rbind(newSegments, newSegment) 
#         nextSegment = nextSegment + 1
#         
#       }else if( (newSegment$End %in% missingRegion$start) | (dataChrSamp$Start[nextSegment] %in% missingRegion$end)){
#         regionIdx = max( which(newSegment$End == missingRegion$start), which(dataChrSamp$Start[nextSegment] == missingRegion$end))
#         newSegment$End = missingRegion$start[regionIdx]
#         newSegments = rbind(newSegments, newSegment)
#         
#         newSegment = dataChrSamp[nextSegment, ]
#         newSegment$Start = missingRegion$end[regionIdx]
#         nextSegment = nextSegment + 1
#         
#       }else if( abs(newSegment$log_R - dataChrSamp$log_R[nextSegment]) <= 1 ){
#         newSegment$End = dataChrSamp$End[nextSegment]
#         newSegment$log_R = ( (newSegment$End - newSegment$Start)*newSegment$log_R + (dataChrSamp$End[nextSegment] - dataChrSamp$Start[nextSegment])*dataChrSamp$log_R[nextSegment])/
#           ((newSegment$End - newSegment$Start) + (dataChrSamp$End[nextSegment] - dataChrSamp$Start[nextSegment]))
#         nextSegment = nextSegment + 1
#         
#       }else{ 
#         newSegments = rbind(newSegments, newSegment) 
#         newSegment = dataChrSamp[nextSegment, ]
#         nextSegment = nextSegment + 1
#       }
#     }
#   }
# }
# 
# dim(test)
# dim(newSegments)
# 
# for(chr in levels(newSegments$Chromosome)){
#   
#   dataChr = newSegments[newSegments$Chromosome==chr,]
#   samplesInt = match(dataChr$Sample,levels(newSegments$Sample))
#   dataChr = as.data.frame(cbind(dataChr, samplesInt))
#   
#   print(
#     ggplot(data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="darkgrey")) + # create empty plot
#       labs(title = chr, x = "Position", y = "Samples", colour = "log_R") +
#       xlim(  (min(dataChr$Start) - 1000), (max(dataChr$End) + 1000)  ) + ylim( c( 0, length(levels(dataChr$Sample))+1) )  +
#       scale_colour_gradient2(limits=c(lowestVal,highestVal), low="blue",mid="white",high="red") +
#       geom_segment(data=dataChr, aes(x=dataChr$Start, xend=dataChr$End, 
#                                      y=dataChr$samplesInt, yend=dataChr$samplesInt, 
#                                      color=dataChr$log_R)) 
#   )
#   
# }
# 



### II.2.d) By extracting all known genes 

# ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "89")
# genesMap <- getBM( attributes = c("ensembl_gene_id","chromosome_name", "start_position", "end_position","external_gene_name"),
#                    mart = ensembl )
# length(unique(genesMap$ensembl_gene_id))
# dim(genesMap)
# genesMap = genesMap[genesMap$chromosome_name %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13",
#                                                     "14","15","16","17","18","19","20","21","22","X"), ]
# length(unique(genesMap$ensembl_gene_id))
# dim(genesMap)
# 
# positions_chr = list()
# for(chr in levels(wgs_raw$Chromosome)){
#   dataChr = wgs_raw[ wgs_raw$Chromosome==chr, ]
#   positions_chr[[chr]] = sort(  c( unique(dataChr$Start), unique(dataChr$End) )  )
# }
# 
# wgs_genes = matrix(data=NA, nrow=length(levels(wgs_raw$Sample)), ncol=nrow(genesMap))
# colnames(wgs_genes) = genesMap$ensembl_gene_id
# rownames(wgs_genes) = rownames(wgs_matrix)
# wgs_genes[1:5, 1:5]
# 
# for(g in 1:nrow(genesMap)){
#   
#   chrName = paste0("chr", genesMap$chromosome_name[g])
#   dataChr = remodelled_chr[[chrName]]
#   positionsChr = positions_chr[[chrName]]
#   
#   pos_start = genesMap$start_position[g]
#   if( any(positionsChr <= pos_start) ){  posIdx_start = tail(which(positionsChr <= pos_start), n=1)  
#   }else{  posIdx_start = NA  }
#   
#   pos_end = genesMap$end_position[g]
#   if( any(positionsChr >= pos_end) ){  posIdx_end = which(positionsChr >= pos_end)[1]  
#   }else{  posIdx_end = NA  }
#   
#   if((is.na(posIdx_start) & posIdx_end==1) | (is.na(posIdx_end) & posIdx_start==length(positionsChr))){  
#     wgs_genes[,g] = NA
#   }else if(is.na(posIdx_start)){  
#     wgs_genes[,g] = dataChr[  , paste(chrName, positionsChr[1], positionsChr[posIdx_end], sep='_')  ]
#   }else if(is.na(posIdx_end)){  
#     wgs_genes[,g] = dataChr[  , paste(chrName, positionsChr[posIdx_start], positionsChr[posIdx_start+1], sep='_')  ]
#   }else if(posIdx_end == (posIdx_start+1)){
#     wgs_genes[,g] = dataChr[  , paste(chrName, positionsChr[posIdx_start], positionsChr[posIdx_end], sep='_')  ]
#   }else if(posIdx_end == (posIdx_start+2)){
#     
#     seg1 = dataChr[, paste(chrName, positionsChr[posIdx_start], positionsChr[posIdx_start+1], sep='_')  ]
#     seg2 = dataChr[, paste(chrName, positionsChr[posIdx_start+1], positionsChr[posIdx_end], sep='_')  ]
#     diffAbs = abs(seg1-seg2)
#     
#     for(samp in 1:nrow(wgs_genes)){
#       
#       if(is.na(diffAbs[samp])){
#         if(any(!is.na(c(seg1[samp],seg2[samp])))){  wgs_genes[samp,g] = c(seg1[samp], seg2[samp])[!is.na(c(seg1[samp], seg2[samp]))]  
#         }else{  wgs_gene[samp, g] = NA  }
#       }else if(diffAbs[samp] < 1){  wgs_genes[samp,g] = mean(c(seg1[samp], seg2[samp]))  
#       }else{  wgs_genes[samp,g] = c(seg1[samp], seg2[samp])[abs(c(seg1[samp], seg2[samp]))==max(abs(c(seg1[samp], seg2[samp])))]  }
#       
#     }
#   }
# }
# 
# sum(is.na(wgs_genes))
# length(colnames(wgs_genes[, which(apply(wgs_genes,2, function(seg) sum(is.na(seg))==nrow(wgs_genes)))]))
# colnames(wgs_genes[, which(apply(wgs_genes,2, function(seg) sum(is.na(seg))==nrow(wgs_genes)))])
# dim(wgs_genes)
# wgs_genesNoNAs = wgs_genes[ , !apply(wgs_genes,2, function(g) sum(is.na(g))==nrow(wgs_genes))  ]
# dim(wgs_genesNoNAs)
} # other unused methods including breaking down by genes, or infering NAs by averaging neighboring segments, or previous GISTIC files


# ===== II.3 Investigate transformed data =====
{
  ## Replace GILA sample IDs with GT labels
  rownames(wgs_matrix) = wgs_samplesMap[match(rownames(wgs_matrix), wgs_samplesMap$DILA_ID_DNA), 1]
  rownames(wgs_matrixNoNAs) = wgs_samplesMap[match(rownames(wgs_matrixNoNAs), wgs_samplesMap$DILA_ID_DNA), 1]
  wgs_matrix[1:5,1:5]
  rownames(wgs_genes) = wgs_samplesMap[match(rownames(wgs_genes), wgs_samplesMap$DILA_ID_DNA), 1]
  rownames(wgs_genesNoNAs) = wgs_samplesMap[match(rownames(wgs_genesNoNAs), wgs_samplesMap$DILA_ID_DNA), 1]
  wgs_genes[1:5,1:5]
  
  ### II.3.a) Compare overlap between transformed and raw data ----
  sampleExample = "GT-03-30-01-WT-02-D"
  dataSamp_raw = wgs_raw[wgs_raw$Sample==wgs_samplesMap[wgs_samplesMap$Gliotrain_ID==sampleExample,2], ]
  
  for(chr in levels(wgs_raw$Chromosome)){
    
    dataChrSamp = dataSamp_raw[dataSamp_raw$Chromosome==chr,]
    chrSegments = positions_chr[[chr]]
    chrGenes = genesMap[paste0("chr",genesMap$chromosome_name)==chr,]
    
    print(
      ggplot(data.frame()) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="darkgrey")) + # create empty plot
        labs(title = chr, x = "Position", colour = "log_R") +
        xlim(  (min(dataChrSamp$Start) - 1000), (max(dataChrSamp$End) + 1000)  ) + 
        scale_colour_gradient2(limits=c(lowestVal,highestVal), low="blue",mid="white",high="red") +
        geom_segment(data=dataChrSamp, aes(x=dataChrSamp$Start, xend=dataChrSamp$End, 
                                           y=rep(0,nrow(dataChrSamp)), yend=rep(0,nrow(dataChrSamp)), 
                                           size=10,
                                           color=dataChrSamp$log_R))+
        geom_vline(xintercept = chrSegments, color="green")+
        geom_segment(data=chrGenes, color="black", aes(x=chrGenes$start_position, xend=chrGenes$end_position, 
                                                       y=rep(0,nrow(chrGenes)), yend=rep(0,nrow(chrGenes)), 
                                                       size=2))
    )
    
  }
  
  
  ### II.3.b) Batch effect ----
  ## Define groups
  wgs_labels = rownames(wgs_genes)
  source.site = c("RCSI","ICM","EMC")[as.integer(sapply(rownames(wgs_genes), substr, 4,5))]
  biocontents = c("DNA","WT","CP")[sapply(rownames(wgs_genes), function(samp) which(substr(samp,13,13)==c("D","W","C")))]
  
  ## wgs_matrix
  # Run PCA
  nb <- estim_ncpPCA(wgs_matrixNoNAs, ncp.max = 10, method.cv = "Kfold", verbose = TRUE)
  nb$ncp
  plot(0:10, nb$criterion, xlab = "nb dim", ylab = "MSEP")
  res.comp_matrix <- imputePCA(wgs_matrixNoNAs, ncp = nb$ncp) # iterativePCA algorithm
  res.comp_matrix$completeObs[1:3,1:3] # the imputed data set
  
  pca_wgs_matrix = dudi.pca(res.comp_matrix, scannf = FALSE, nf = 10)
  pve <- 100*pca_wgs_matrix$eig/sum(pca_wgs_matrix$eig)
  cumsum(pve)
  par(mfrow=c(2,2))
  s.label(pca_wgs_matrix$li, xax=1, yax=2)
  s.label(pca_wgs_matrix$li, xax=3, yax=4)
  s.label(pca_wgs_matrix$li, xax=5, yax=6)
  s.label(pca_wgs_matrix$li, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  
  # Visualize distribution for source sites
  gcol = c("red","blue","green")
  par(mfrow=c(2,2))
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(source.site), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(source.site), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(source.site), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(source.site), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for shipped content
  par(mfrow=c(2,2))
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(biocontents), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(biocontents), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(biocontents), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_wgs_matrix$li, fac=as.factor(biocontents), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  
  ## wgs_genes
  # Run PCA
  nb <- estim_ncpPCA(wgs_genesNoNAs, ncp.max = 10, method.cv = "Kfold", verbose = TRUE)
  nb$ncp
  plot(0:10, nb$criterion, xlab = "nb dim", ylab = "MSEP")
  res.comp_genes <- imputePCA(wgs_genesNoNAs, ncp = nb$ncp) # iterativePCA algorithm
  res.comp_genes$completeObs[1:3,1:3] # the imputed data set
  
  pca_wgs_genes = dudi.pca(res.comp_genes, scannf = FALSE, nf = 10)
  pve <- 100*pca_wgs_genes$eig/sum(pca_wgs_genes$eig)
  cumsum(pve)
  par(mfrow=c(2,2))
  s.label(pca_wgs_genes$li, xax=1, yax=2)
  s.label(pca_wgs_genes$li, xax=3, yax=4)
  s.label(pca_wgs_genes$li, xax=5, yax=6)
  s.label(pca_wgs_genes$li, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for source sites
  gcol = c("red","blue","green")
  par(mfrow=c(2,2))
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(source.site), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(source.site), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(source.site), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(source.site), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for shipped content
  par(mfrow=c(2,2))
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(biocontents), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(biocontents), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(biocontents), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_wgs_genes$li, fac=as.factor(biocontents), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
} # end of II.3 Investigate transformed data


# ===== II.4 Produce tranSMART ready files =====

### II.4.a) "Raw" data ----
# Raw WGS data will be stored with fixed-size segments broken down and common to all samples, 
#    as if it was microarray data.

### Data files
rawWGS_transmart = t(wgs_matrixTransmart)
colnames(rawWGS_transmart) = wgs_samplesMap$Gliotrain_ID[match(colnames(rawWGS_transmart), wgs_samplesMap$DILA_ID_DNA)]
rawWGS_transmart = cbind(ID_REF = rownames(rawWGS_transmart), as.data.frame(rawWGS_transmart))
rawWGS_transmart[1:5,1:5]
dim(rawWGS_transmart)

# Divide between parental tumors and cell cultures
rawWGS_transmart_CP = rawWGS_transmart[, c(1, grep("CP",colnames(rawWGS_transmart)))]
rawWGS_transmart_PT = rawWGS_transmart[, -grep("CP",colnames(rawWGS_transmart))]
dim(rawWGS_transmart_CP)
dim(rawWGS_transmart_PT)


### Annotation file
annotations_WGSraw = data.frame(GPL_ID = rep("WGSrawGTdata", nrow(rawWGS_transmart)),
                                PROBE_ID = rawWGS_transmart[,1],
                                GENE_SYMBOL = rep("", nrow(rawWGS_transmart)),
                                GENE_ID = rep("", nrow(rawWGS_transmart)),
                                Organism = rep("Homo_sapiens", nrow(rawWGS_transmart)))
head(annotations_WGSraw)


### Subject Samples Mapping file
samples_WGSraw = colnames(rawWGS_transmart)[2:ncol(rawWGS_transmart)]
SubjSampMap_WGSraw = data.frame(STUDY_ID = rep("GTdata", length(samples_WGSraw)),
                                SITE_ID = rep("", length(samples_WGSraw)),
                                SUBJECT_ID = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", samples_WGSraw),
                                SAMPLE_ID = samples_WGSraw,
                                PLATFORM = rep("WGSrawGTdata", length(samples_WGSraw)),
                                TISSUETYPE = rep("Homo_sapiens", length(samples_WGSraw)),
                                ATTR1 = rep("", length(samples_WGSraw)),
                                ATTR2 = rep("", length(samples_WGSraw)),
                                CATEGORY_CD = rep("Biomarker_Data+Genomics+WGS+logR_values", length(samples_WGSraw)),
                                SOURCE_CD= rep("STD", length(samples_WGSraw)))
head(SubjSampMap_WGSraw)

# Divide between parental tumors and cell cultures
SubjSampMap_WGSraw_CP = SubjSampMap_WGSraw[ grep("CP", SubjSampMap_WGSraw$SAMPLE_ID), ]
SubjSampMap_WGSraw_CP$CATEGORY_CD = rep("Biomarker_Data+Genomics+WGS+logR_values_Cell_cultures", nrow(SubjSampMap_WGSraw_CP))
SubjSampMap_WGSraw_PT = SubjSampMap_WGSraw[ -grep("CP", SubjSampMap_WGSraw$SAMPLE_ID), ]
dim(SubjSampMap_WGSraw_CP)
dim(SubjSampMap_WGSraw_PT)


### Batches file
wgs_samplesBatch = read.csv("./tmp/samplesBatches_WGS.csv", h=T)
batches = wgs_samplesBatch$FCID[match(wgs_samplesMap[match(samples_WGSraw, wgs_samplesMap$Gliotrain_ID), 2], wgs_samplesBatch$SampleID)]
batches = paste0("Batch",as.numeric(batches))
lanes = wgs_samplesBatch$Lane[match(wgs_samplesMap[match(samples_WGSraw, wgs_samplesMap$Gliotrain_ID), 2], wgs_samplesBatch$SampleID)]
lanes = paste0("lane", as.numeric(lanes))
batches_wgs = data.frame(patient_id = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", samples_WGSraw),
                         sample_id = samples_WGSraw,
                         WGS_batches = batches, 
                         WGS_lanes_within_batches = lanes,
                         stringsAsFactors = FALSE)

# Divide between parental tumors and cell cultures
batches_wgs_CP = batches_wgs[grep("CP", samples_WGSraw),]
batches_wgs_PT = batches_wgs[-grep("CP", samples_WGSraw),]
dim(batches_wgs_CP)
dim(batches_wgs_PT)


### Save files
# Parental tumors
write.table(rawWGS_transmart_PT, file = "./readyForUpload/parentalTumors/mRNA/raw_gex_matrix_mRNA_L.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_WGSraw, file = "./readyForUpload/parentalTumors/mRNA/WGSrawGTdata.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(SubjSampMap_WGSraw_PT, file = "./readyForUpload/parentalTumors/mRNA/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(batches_wgs_PT, file = "./readyForUpload/parentalTumors/transmart/data_WGSBatches_PT.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")

# Cell cultures
write.table(rawWGS_transmart_CP, file = "./readyForUpload/cellCultures/mRNA/raw_gex_matrix_mRNA_L.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_WGSraw, file = "./readyForUpload/cellCultures/mRNA/WGSrawGTdata.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(SubjSampMap_WGSraw_CP, file = "./readyForUpload/cellCultures/mRNA/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(batches_wgs_CP, file = "./readyForUpload/parentalTumors/transmart/data_WGSBatches_CP.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")


### II.4.b) GISTIC-identified events data ----
# GISTIC output CNV will be sotred as cnv, with the amplified/deleted category {-1,0,1,2} in .flag columns
#   and the raw log_R value in .chip columns.

### Data files
## Copy number categories
gisticWGS.flag_transmart = gistic_ampDels
colnames(gisticWGS.flag_transmart) = paste0(wgs_samplesMap$Gliotrain_ID[match(colnames(gisticWGS.flag_transmart), wgs_samplesMap$DILA_ID_DNA)], ".flag")
gisticWGS.flag_transmart = cbind(region_name = rownames(gisticWGS.flag_transmart), as.data.frame(gisticWGS.flag_transmart))
gisticWGS.flag_transmart[1:5,1:5]
dim(gisticWGS.flag_transmart)

# Divide between parental tumors and cell cultures
gisticWGS.flag_transmart_CP = gisticWGS.flag_transmart[, c(1, grep("CP",colnames(gisticWGS.flag_transmart)))]
gisticWGS.flag_transmart_PT = gisticWGS.flag_transmart[, -grep("CP",colnames(gisticWGS.flag_transmart))]
dim(gisticWGS.flag_transmart_CP)
dim(gisticWGS.flag_transmart_PT)

## Identify the corresponding log_R values from the raw data
## given up, that would require some weighed averaging since it concerns focal events, meaning very heterogenously segmented regions.
{# Breakdown available regions in the raw data.
# rawReg = sapply(rownames(rawWGS_transmart), 
#                 function(region) c(strsplit(region, '_')[[1]][1], 
#                                    strsplit(region, '_')[[1]][2], 
#                                    strsplit(region, '_')[[1]][3]))
# rawRegions = data.frame(chrom = as.character(rawReg[1,]),
#                         start = as.numeric(as.character(rawReg[2,])),
#                         end = as.numeric(as.character(rawReg[3,])))
# head(rawRegions)
# 
# # Map gistic region to raw region
# for(region in rownames(gisticWGS.flag_transmart)){
#   cat("testing: ", region, '\n')
#   region.chr = strsplit(region, '_')[[1]][1]
#   region.start = as.numeric(strsplit(region, '_')[[1]][2])
#   region.end = as.numeric(strsplit(region, '_')[[1]][3])
#   
#   rawRegs = rawRegions[rawRegions$chrom == region.chr &
#                          rawRegions$end > region.start &
#                          rawRegions$start < region.end,]
#   dim(rawRegs)
#   
#   
#   if(any(as.character(rawRegions[,1]) == region.chr & rawRegions[,2] <= region.start & rawRegions[,3] >= region.end)){
#     print(rownames(rawRegions)[(as.character(rawRegions[,1]) == region.chr & rawRegions[,2] <= region.start & rawRegions[,3] >= region.end)])
#   }
#   
# }
}


### Annotation file
annotations_WGSgistic = data.frame(GPL_ID = rep("WGSGisticCNVGTdata", nrow(gisticWGS.flag_transmart)),
                                   REGION_NAME = gisticWGS.flag_transmart[,1],
                                   CHROMOSOME = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   START_BP = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   END_BP = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   NUM_PROBES = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   CYTOBAND = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   GENE_SYMBOL = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   GENE_ID = rep(NA, nrow(gisticWGS.flag_transmart)),
                                   ORGANISM = rep("Homo sapiens", nrow(gisticWGS.flag_transmart)))
for(region in 1:nrow(annotations_WGSgistic)){
  annotations_WGSgistic$CHROMOSOME[region] = gsub("chr", "", strsplit(as.character(annotations_WGSgistic$REGION_NAME[region]), '_')[[1]][1])
  annotations_WGSgistic$START_BP[region] = strsplit(as.character(annotations_WGSgistic$REGION_NAME[region]), '_')[[1]][2]
  annotations_WGSgistic$END_BP[region] = strsplit(as.character(annotations_WGSgistic$REGION_NAME[region]), '_')[[1]][3]
}
head(annotations_WGSgistic)


### Subject Samples Mapping files
samples_WGSgistic = colnames(gisticWGS.flag_transmart)[2:ncol(gisticWGS.flag_transmart)]
SubjSampMap_WGSgistic = data.frame(trial_name = rep("GTdata", length(samples_WGSgistic)),
                                   site_id = rep("", length(samples_WGSgistic)),
                                   subject_id = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+\\.flag$", "", samples_WGSgistic),
                                   sample_cd = gsub("\\.flag$", "", samples_WGSgistic),
                                   platform = rep("WGSGisticCNVGTdata", length(samples_WGSgistic)),
                                   sample_type = rep("", length(samples_WGSgistic)),
                                   tissuetype = rep("Homo sapiens", length(samples_WGSgistic)),
                                   time_point = rep("", length(samples_WGSgistic)),
                                   cat_cd = rep("Biomarker_Data+Genomics+WGS+focal_events", length(samples_WGSgistic)),
                                   src_cd= rep("STD", length(samples_WGSgistic)))
head(SubjSampMap_WGSgistic)

# Divide between parental tumors and cell cultures
SubjSampMap_WGSgistic_CP = SubjSampMap_WGSgistic[ grep("CP", SubjSampMap_WGSgistic$sample_cd), ]
SubjSampMap_WGSgistic_CP$cat_cd = rep("Biomarker_Data+Genomics+WGS+focal_events_Cell_cultures", nrow(SubjSampMap_WGSgistic_CP))
SubjSampMap_WGSgistic_PT = SubjSampMap_WGSgistic[ -grep("CP", SubjSampMap_WGSgistic$sample_cd), ]
dim(SubjSampMap_WGSgistic_CP)
dim(SubjSampMap_WGSgistic_PT)


### Save files
# Parental tumors
write.table(gisticWGS.flag_transmart_PT, file = "./readyForUpload/parentalTumors/CNA/cnv/cnv_data.tsv", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_WGSgistic, file = "./readyForUpload/parentalTumors/CNA/CNV_ANNOT/chromosomal_regions.tsv", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")
write.table(SubjSampMap_WGSgistic_PT, file = "./readyForUpload/parentalTumors/CNA/cnv/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")

# Cell cultures
write.table(gisticWGS.flag_transmart_CP, file = "./readyForUpload/cellCultures/CNA/cnv/cnv_data.tsv", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_WGSgistic, file = "./readyForUpload/cellCultures/CNA/CNV_ANNOT/chromosomal_regions.tsv", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")
write.table(SubjSampMap_WGSgistic_CP, file = "./readyForUpload/cellCultures/CNA/cnv/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")




# oooooooooooooooooooooooooooooooooooooooooooooooooo
# III - RNA Seq Data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

rnaseq_raw = read.csv("./tmp/rnaSeq_raw_resequencedUpdated_200110.csv", header = T, row.names = 1)
rnaseq_samplesMap = read.csv("./tmp/samplesLabels_RNASeq.csv", h=T)
rnaseq_samplesBatch = read.csv("./tmp/samplesBatches_RNASeq.csv", h=T)
rnaseq_annots = read.table("./tmp/conversion_rnaSeqIDs.txt", h=T, sep='\t')

# ===== III.1 Preliminary transformations =====
rnaseq_raw[1:5, 1:5]

## Transpose
rnaseq_rawt = t(rnaseq_raw)
rnaseq_rawt[1:5, 1:5]


## Replace platform ID with GLIOTRAIN ID
rownames(rnaseq_rawt) = as.character(rnaseq_samplesMap$Gliotrain_ID)[  match(rownames(rnaseq_rawt), rnaseq_samplesMap$DILA_ID_RNA)  ]
rnaseq_rawt[1:5, 1:5]


# ===== III.2 Map ENSEMBL IDs to HGNC and ENTREZ standards =====
## Deprecated, a file mapping ENSMBLIDs to "gene symbols" was provided by Gonca
{
# ### Replace ENSEMBL IDs with ENTREZ IDs
# ensmblIDs = gsub("\\.\\d+", "", colnames(rnaseq_rawt))
# ensmblIDs[1:5]
# 
# 
# ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "89")
# genesMap <- getBM( filters = "ensembl_gene_id", values = ensmblIDs, mart = ensembl, 
#                    attributes = c("ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_name", "chromosome_name", "start_position", "end_position"))  #, "version"
# 
# length(ensmblIDs)
# dim(genesMap)
# head(genesMap)
# ensmblMinMap = ensmblIDs[ ! ( ensmblIDs %in% genesMap$ensembl_gene_id ) ]
# length(ensmblMinMap)
# ensmblMinMap
# genesMissingInfo = genesMap[(  is.na(genesMap$entrezgene) | genesMap$hgnc_symbol == ""),  ]
# cat((nrow(genesMap)-nrow(genesMissingInfo)), (nrow(genesMissingInfo)/nrow(genesMap)), "\n")
# 
# 
# ### Use org.HS.db to complete missing values
# ## For HGNC
# noHGNC = genesMap[genesMap$hgnc_symbol=="",]
# dim(noHGNC)
# length(unique(noHGNC$ensembl_gene_id))
# head(noHGNC)
# 
# HGNCmissing = mapIds(org.Hs.eg.db, noHGNC$ensembl_gene_id, "SYMBOL", "ENSEMBL")
# head(HGNCmissing)
# HGNCmissing = HGNCmissing[!is.na(HGNCmissing)]
# head(HGNCmissing)
# genesMapHGNC = genesMap
# genesMapHGNC[match(names(HGNCmissing), genesMapHGNC$ensembl_gene_id), "hgnc_symbol"] = HGNCmissing
# 
# noHGNC = genesMapHGNC[genesMapHGNC$hgnc_symbol == "",]
# dim(noHGNC)
# length(unique(noHGNC$ensembl_gene_id))
# head(noHGNC)
# 
# HGNCmissing = mapIds(org.Hs.eg.db, noHGNC$external_gene_name, "SYMBOL", "ALIAS")
# head(HGNCmissing)
# HGNCmissing = HGNCmissing[!is.na(HGNCmissing)]
# head(HGNCmissing)
# genesMapHGNC[match(names(HGNCmissing), genesMapHGNC$external_gene_name), "hgnc_symbol"] = HGNCmissing
# 
# noHGNC = genesMapHGNC[genesMapHGNC$hgnc_symbol == "",]
# dim(noHGNC)
# length(unique(noHGNC$ensembl_gene_id))
# head(noHGNC)
# 
# ## For ENTREZ
# noENTREZID = genesMapHGNC[is.na(genesMapHGNC$entrezgene),]
# dim(noENTREZID)
# length(unique(noENTREZID$ensembl_gene_id))
# head(noENTREZID)
# 
# ENTREZmissing = mapIds(org.Hs.eg.db, noENTREZID$hgnc_symbol, "ENTREZID", "SYMBOL")
# ENTREZmissing = ENTREZmissing[!is.na(names(ENTREZmissing))]
# head(ENTREZmissing)
# genesMap_HGNC_ENTREZ = genesMapHGNC
# genesMap_HGNC_ENTREZ[match(names(ENTREZmissing), genesMapHGNC$hgnc_symbol), "entrezgene"] = as.numeric(ENTREZmissing)
# 
# noENTREZID = genesMap_HGNC_ENTREZ[is.na(genesMap_HGNC_ENTREZ$entrezgene),]
# dim(noENTREZID)
# length(unique(noENTREZID$ensembl_gene_id))
# head(noENTREZID)
# 
# ENTREZmissing = mapIds(org.Hs.eg.db, noENTREZID$external_gene_name, "ENTREZID", "ALIAS")
# ENTREZmissing = ENTREZmissing[!is.na(names(ENTREZmissing))]
# head(ENTREZmissing)
# genesMap_HGNC_ENTREZ[match(names(ENTREZmissing), genesMap_HGNC_ENTREZ$external_gene_name), "entrezgene"] = as.numeric(ENTREZmissing)
# 
# noENTREZID = genesMap_HGNC_ENTREZ[is.na(genesMap_HGNC_ENTREZ$entrezgene),]
# dim(noENTREZID)
# length(unique(noENTREZID$ensembl_gene_id))
# head(noENTREZID)
# 
# genesMapFinal = genesMap_HGNC_ENTREZ
# genesMissingInfoFinal = genesMapFinal[(  is.na(genesMapFinal$entrezgene) | genesMapFinal$hgnc_symbol == ""),  ]
# dim(genesMissingInfoFinal)
# completeMappingLength = length(unique(genesMapFinal$ensembl_gene_id))
# missingEntrezOrHGNCLength = length(unique(genesMissingInfoFinal$ensembl_gene_id))
# cat("Total Number of transcripts = ", completeMappingLength, "\n",
#     "Number of successfully and unsuccessfully mapped IDs: ", completeMappingLength - missingEntrezOrHGNCLength, missingEntrezOrHGNCLength, "\n",
#     "Percentage of successfully and unsuccessfully mapped IDs: ", (completeMappingLength - missingEntrezOrHGNCLength)/completeMappingLength, missingEntrezOrHGNCLength/completeMappingLength, "\n")
# length(unique(genesMapFinal$hgnc_symbol[!(genesMapFinal$hgnc_symbol=="")]))
# length(unique(genesMapFinal$entrezgene[!is.na(genesMapFinal$entrezgene)]))
}


# ===== III.3 Investigate transformed data =====

### PCA
{
  
  # Define groups
  rna_labels = rownames(rnaseq_rawt)
  source.site = c("RCSI","ICM","EMC")[as.integer(sapply(rna_labels, substr, 4,5))]
  biocontents = c("RNA","WT","CP")[sapply(rna_labels, function(samp) which(substr(samp,13,13)==c("R","W","C")))]
  batches = rnaseq_samplesBatch$fcid[match(rnaseq_samplesMap[match(rna_labels, rnaseq_samplesMap$Gliotrain_ID), 2], rnaseq_samplesBatch$sampleid)]
  batches = as.factor(paste0("Batch", as.numeric(batches)))
  
  # Log dataset to reduce bias
  rnaseq_pcaset = rnaseq_rawt
  rnaseq_pcaset[1:5, 1:5]
  for(i in 1:nrow(rnaseq_pcaset)){ for(j in 1:ncol(rnaseq_pcaset)){  if(rnaseq_pcaset[i, j] == 0){   rnaseq_pcaset[i, j] = 1   }  }}
  rnaseq_pcaset = log(rnaseq_pcaset)
  rnaseq_pcaset[1:5, 1:5]
  
  # Run PCA
  pca_rnaSeq = dudi.pca(rnaseq_pcaset, scannf = FALSE, nf = 10)
  pve <- 100*pca_rnaSeq$eig/sum(pca_rnaSeq$eig)
  cumsum(pve)
  par(mfrow=c(2,2))
  s.label(pca_rnaSeq$li, xax=1, yax=2)
  s.label(pca_rnaSeq$li, xax=3, yax=4)
  s.label(pca_rnaSeq$li, xax=5, yax=6)
  s.label(pca_rnaSeq$li, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  
  # Visualize distribution for source sites
  gcol = c("red","blue","green")
  par(mfrow=c(2,2))
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(source.site), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(source.site), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(source.site), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(source.site), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for shipped content
  par(mfrow=c(2,2))
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(biocontents), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(biocontents), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(biocontents), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_rnaSeq$li, fac=as.factor(biocontents), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for batches
  gcol = c("red","blue","green", "purple")
  par(mfrow=c(2,2))
  s.class(dfxy = pca_rnaSeq$li, fac=batches, col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_rnaSeq$li, fac=batches, col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_rnaSeq$li, fac=batches, col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_rnaSeq$li, fac=batches, col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  
  # Run PCA without outliers
  rnaSeq_noOutlier = rnaseq_pcaset[!(rownames(rnaseq_pcaset) %in% c("GT-02-39-01-R-02-R")),]#, "GT-01-15-01-R-02-R", "GT-03-40-02-CP-02-R")),] # commented out are outliers for non-logged dataset
  pca_rnaSeqNoOut = dudi.pca(rnaSeq_noOutlier, scannf = FALSE, nf = 10)
  pve <- 100*pca_rnaSeqNoOut$eig/sum(pca_rnaSeqNoOut$eig)
  cumsum(pve)
  par(mfrow=c(2,2))
  s.label(pca_rnaSeqNoOut$li, xax=1, yax=2)
  s.label(pca_rnaSeqNoOut$li, xax=3, yax=4)
  s.label(pca_rnaSeqNoOut$li, xax=5, yax=6)
  s.label(pca_rnaSeqNoOut$li, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  #Define groups
  wgs_labelsNoOut = rownames(rnaSeq_noOutlier)
  source.siteNoOut = c("RCSI","ICM","EMC")[as.integer(sapply(wgs_labelsNoOut, substr, 4,5))]
  biocontentsNoOut = c("RNA","WT","CP")[sapply(wgs_labelsNoOut, function(samp) which(substr(samp,13,13)==c("R","W","C")))]
  batchesNoOut = rnaseq_samplesBatch$fcid[match(rnaseq_samplesMap[match(wgs_labelsNoOut, rnaseq_samplesMap$Gliotrain_ID), 2], rnaseq_samplesBatch$sampleid)]
  batchesNoOut = as.factor(paste0("Batch", as.numeric(batchesNoOut)))
  
  # Visualize distribution for source sites
  gcol = c("red","blue","green")
  par(mfrow=c(2,2))
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(source.siteNoOut), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(source.siteNoOut), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(source.siteNoOut), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(source.siteNoOut), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for shipped content
  par(mfrow=c(2,2))
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(biocontentsNoOut), col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(biocontentsNoOut), col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(biocontentsNoOut), col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=as.factor(biocontentsNoOut), col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
  # Visualize distribution for batches
  gcol = c("red","blue","green", "purple")
  par(mfrow=c(2,2))
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=batchesNoOut, col=gcol, xax=1, yax=2)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=batchesNoOut, col=gcol, xax=3, yax=4)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=batchesNoOut, col=gcol, xax=5, yax=6)
  s.class(dfxy = pca_rnaSeqNoOut$li, fac=batchesNoOut, col=gcol, xax=7, yax=8)
  par(mfrow=c(1,1))
  
} # end of PCA exploration


### t-SNE
{

rnaseq_noConstVar = rnaseq_pcaset[, which(apply(rnaseq_pcaset, 2, var) != 0)]
cat(ncol(rnaseq_noConstVar), ncol(rnaseq_rawt))
tsne_rawt = Rtsne(rnaseq_noConstVar, dims=3, pca_scale=T)
d_tsne = as.data.frame(tsne_rawt$Y)
rownames(d_tsne) = rownames(rnaseq_rawt)
dim(d_tsne)

## Vizualise for source site
Site = as.factor(source.site)

ggplot(cbind(d_tsne,Site), aes(x=V1, y=V2, colour=Site)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("t-SNE v1 v2")

ggplot(cbind(d_tsne,Site), aes(x=V1, y=V3, colour=Site)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("t-SNE v1 v3")

ggplot(cbind(d_tsne,Site), aes(x=V2, y=V3, colour=Site)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("t-SNE v2 v3")

## Vizualise for biological content
Contents = as.factor(biocontents)

ggplot(cbind(d_tsne,Contents), aes(x=V1, y=V2, colour=Contents)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("t-SNE v1 v2")

ggplot(cbind(d_tsne,Contents), aes(x=V1, y=V3, colour=Contents)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("t-SNE v1 v3")

ggplot(cbind(d_tsne,Contents), aes(x=V2, y=V3, colour=Contents)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("t-SNE v2 v3")


## Vizualise for batches 
Batch = as.factor(batches)

ggplot(cbind(d_tsne,Batch), aes(x=V1, y=V2, colour=Batch)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green", "purple", "black")) +
  xlab("") + ylab("") + ggtitle("t-SNE v1 v2")

ggplot(cbind(d_tsne,Batch), aes(x=V1, y=V3, colour=Batch)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green", "purple", "black")) +
  xlab("") + ylab("") + ggtitle("t-SNE v1 v3")

ggplot(cbind(d_tsne,Batch), aes(x=V2, y=V3, colour=Batch)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green", "purple", "black")) +
  xlab("") + ylab("") + ggtitle("t-SNE v2 v3")

} # end of t-SNE eploration


### UMAP
{
umap_rawt = umap(rnaseq_pcaset, n_components = 4)
names(umap_rawt)
head(umap_rawt$layout)
umap_comps = as.data.frame(umap_rawt$layout)
colnames(umap_comps) = c("V1","V2","V3","V4")

## Vizualise for source site
Site = as.factor(source.site)

ggplot(cbind(umap_comps,Site), aes(x=V1, y=V2, colour=Site)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("UMAP v1 v2")

ggplot(cbind(umap_comps,Site), aes(x=V3, y=V4, colour=Site)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("UMAP v3 v4")


## Vizualise for biological content
Contents = as.factor(biocontents)

ggplot(cbind(umap_comps,Contents), aes(x=V1, y=V2, colour=Contents)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("UMAP v1 v2")

ggplot(cbind(umap_comps,Contents), aes(x=V3, y=V4, colour=Contents)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green")) +
  xlab("") + ylab("") + ggtitle("UMAP v3 v4")

## Vizualise for batches 
Batch = as.factor(batches)

ggplot(cbind(umap_comps,Batch), aes(x=V1, y=V2, colour=Batch)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green", "purple", "black")) +
  xlab("") + ylab("") + ggtitle("UMAP v1 v2")

ggplot(cbind(umap_comps,Batch), aes(x=V3, y=V4, colour=Batch)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"), legend.position = "bottom") +
  geom_point() + scale_color_manual(values = c("red","blue","green", "purple", "black")) +
  xlab("") + ylab("") + ggtitle("UMAP v3 v4")

} # end of UMAP exmploration


### GLM-PCA
{
  samp_labels = rownames(rnaseq_rawt)
  source.site = c("RCSI","ICM","EMC")[as.integer(sapply(rownames(rnaseq_rawt), substr, 4,5))]
  biocontents = c("RNA","WT","CP")[sapply(rownames(rnaseq_rawt), function(samp) which(substr(samp,13,13)==c("R","W","C")))]
  batches = rnaseq_samplesBatch$fcid[match(rnaseq_samplesMap[match(samp_labels, rnaseq_samplesMap$Gliotrain_ID), 2], rnaseq_samplesBatch$sampleid)]
  batches = as.factor(paste0("Batch",as.numeric(batches)))
  names(source.site) <- names(biocontents) <- names(batches) <- samp_labels
  
  rnaseq_glmPCAset = t(rnaseq_pcaset)
  dim(rnaseq_glmPCAset)
  rnaseq_glmPCAset = rnaseq_glmPCAset[!apply(rnaseq_glmPCAset,1, function(t) sum(t) == 0),]
  dim(rnaseq_glmPCAset)
  rnaseq_glmPCAset[1:5, 1:5]
  
  gpca <- glmpca(Y = rnaseq_glmPCAset, L = 10, verbose = T)
  gpca.dat <- gpca$factors
  gpca.dat$site <- source.site
  gpca.dat$contents <- biocontents
  gpca.dat$batches <- batches
  
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = site, shape = contents)) +
    geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
  ggplot(gpca.dat, aes(x = dim3, y = dim4, color = site, shape = contents)) +
    geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
  
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = site, shape = batches)) +
    geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
  ggplot(gpca.dat, aes(x = dim3, y = dim4, color = site, shape = batches)) +
    geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
  
} # end of GLM-PCA exploration


# ===== III.4 Produce transMART files =====

### Data files
rnaseq_transmart = rnaseq_raw
colnames(rnaseq_transmart) = rnaseq_samplesMap$Gliotrain_ID[match(colnames(rnaseq_transmart), rnaseq_samplesMap$DILA_ID_RNA)]
rnaseq_transmart = cbind(ID_REF = rownames(rnaseq_transmart), as.data.frame(rnaseq_transmart))
rnaseq_transmart[1:5,1:5]

# Divide between parental tumors and cell cultures
rnaseq_transmart_CP = rnaseq_transmart[, c(1, grep("CP",colnames(rnaseq_transmart)))]
rnaseq_transmart_PT = rnaseq_transmart[, -grep("CP",colnames(rnaseq_transmart))]
dim(rnaseq_transmart_CP)
dim(rnaseq_transmart_PT)

### Annotation file
annotations_rnaseq = data.frame(Transcript_ID = rnaseq_annots$id,
                                Gene_Symbol = rnaseq_annots$geneID,
                                Organism = rep("Homo_sapiens", nrow(rnaseq_annots)))
head(annotations_rnaseq)

### Subject Samples Mapping files
samples_rnaseq = colnames(rnaseq_transmart)[2:ncol(rnaseq_transmart)]
SubjSampMap_rnaseq = data.frame(STUDY_ID = rep("GTdata", length(samples_rnaseq)),
                                SITE_ID = rep("", length(samples_rnaseq)),
                                SUBJECT_ID = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", samples_rnaseq),
                                SAMPLE_ID = samples_rnaseq,
                                PLATFORM = rep("MRNASeqGTdata", length(samples_rnaseq)),
                                TISSUETYPE = rep("Homo_sapiens", length(samples_rnaseq)),
                                ATTR1 = rep("", length(samples_rnaseq)),
                                ATTR2 = rep("", length(samples_rnaseq)),
                                CATEGORY_CD = rep("Biomarker_Data+Gene_Expression+RNA-Seq", length(samples_rnaseq)),
                                SOURCE_CD = rep("STD", length(samples_rnaseq)))
head(SubjSampMap_rnaseq)

# Divide between parental tumors and cell cultures
SubjSampMap_rnaseq_CP = SubjSampMap_rnaseq[ grep("CP", SubjSampMap_rnaseq$SAMPLE_ID), ]
SubjSampMap_rnaseq_CP$CATEGORY_CD = rep("Biomarker_Data+Gene_Expression+RNA-Seq_Cell_cultures", nrow(SubjSampMap_rnaseq_CP))
SubjSampMap_rnaseq_PT = SubjSampMap_rnaseq[ -grep("CP", SubjSampMap_rnaseq$SAMPLE_ID), ]
dim(SubjSampMap_rnaseq_CP)
dim(SubjSampMap_rnaseq_PT)


### Batch file
batches = rnaseq_samplesBatch$fcid[match(rnaseq_samplesMap[match(samples_rnaseq, rnaseq_samplesMap$Gliotrain_ID), 2], rnaseq_samplesBatch$sampleid)]
batches = paste0("Batch", as.numeric(batches))
lanes = rnaseq_samplesBatch$lane[match(rnaseq_samplesMap[match(samples_rnaseq, rnaseq_samplesMap$Gliotrain_ID), 2], rnaseq_samplesBatch$sampleid)]
lanes = paste0("lane", as.numeric(lanes))
batches_rnaseq = data.frame(patient_id = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", samples_rnaseq),
                            sample_id = samples_rnaseq,
                            RNASeq_batches = batches, 
                            RNASeq_lanes_within_batches = lanes,
                            stringsAsFactors = FALSE)

# Divide between parental tumors and cell cultures
batches_rnaseq_CP = batches_rnaseq[grep("CP", samples_rnaseq),]
batches_rnaseq_PT = batches_rnaseq[-grep("CP", samples_rnaseq),]
dim(batches_rnaseq_CP)
dim(batches_rnaseq_PT)


### Save files
# Parental tumors
write.table(rnaseq_transmart_PT, file = "./readyForUpload/parentalTumors/mRNA_RNASeq/raw_gex_matrix_mRNA_RNASeq_R.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_rnaseq, file = "./readyForUpload/parentalTumors/mRNA_RNASeq/MRNASeqGTdata_annotation.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(SubjSampMap_rnaseq_PT, file = "./readyForUpload/parentalTumors/mRNA_RNASeq/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(batches_rnaseq_PT, file = "./readyForUpload/parentalTumors/transmart/data_RNASeqBatches_PT.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")

# Cell Cultures
write.table(rnaseq_transmart_CP, file = "./readyForUpload/cellCultures/mRNA_RNASeq/raw_gex_matrix_mRNA_RNASeq_R.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_rnaseq, file = "./readyForUpload/cellCultures/mRNA_RNASeq/MRNASeqGTdata_annotation.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(SubjSampMap_rnaseq_CP, file = "./readyForUpload/cellCultures/mRNA_RNASeq/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(batches_rnaseq_CP, file = "./readyForUpload/parentalTumors/transmart/data_RNASeqBatches_CP.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")



# oooooooooooooooooooooooooooooooooooooooooooooooooo
# IV - Methylation data ----
# oooooooooooooooooooooooooooooooooooooooooooooooooo

mt_raw = read.csv("./rawData/CpGs_GLIOTRAIN.csv", h = T, sep = '\t', stringsAsFactors = F)
mt_samplesMap = read.csv("./tmp/samplesLabels_Methylation.csv", h = T, stringsAsFactors = F)
mt_batches = read.csv("./tmp/samplesBatches_Methylation.csv", h = T, stringsAsFactors = F)

# ===== IV.1 Get a sense of the data =====

dim(mt_raw)
mt_raw[1:5, 1:5]
sum(mt_raw$start != mt_raw$end)
sum(is.na(mt_raw))
na_per_position = apply(mt_raw, 1, function(pos) sum(is.na(pos)))
hist(na_per_position)
na_per_sample = apply(mt_raw, 2, function(samp) sum(is.na(samp)))
hist(na_per_sample, 12)

head(mt_samplesMap)
head(mt_batches)


# ===== IV.2 Make transformations =====

mt_data = data.frame(ID_REF = paste0(mt_raw$chr, '_', mt_raw$start),
                     mt_raw[, 4:18])
colnames(mt_data) = c("ID_REF", mt_samplesMap$Gliotrain_ID[match(names(mt_data)[2:16], mt_samplesMap$DILA_ID_Methylation)])
mt_data[1:5, 1:5]

mt_batches$SampleID = mt_samplesMap$Gliotrain_ID[match(mt_batches$SampleID, mt_samplesMap$DILA_ID_Methylation)]


# ===== IV.3 Produce tranSMART ready files =====

### Annotation file
annotations_mt = data.frame(GPL_ID = rep("MethylationGTdata", nrow(mt_data)),
                            PROBE_ID = mt_data$ID_REF,
                            GENE_SYMBOL = rep("", nrow(mt_data)),
                            GENE_ID = rep("", nrow(mt_data)),
                            Organism = rep("Homo_sapiens", nrow(mt_data)))
head(annotations_mt)


### Subject Samples Mapping file
samples_mt = colnames(mt_data)[2:ncol(mt_data)]
SubjSampMap_mt = data.frame(STUDY_ID = rep("GTdata", length(samples_mt)),
                            SITE_ID = rep("", length(samples_mt)),
                            SUBJECT_ID = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", samples_mt),
                            SAMPLE_ID = samples_mt,
                            PLATFORM = rep("MethylationGTdata", length(samples_mt)),
                            TISSUETYPE = rep("Homo_sapiens", length(samples_mt)),
                            ATTR1 = rep("", length(samples_mt)),
                            ATTR2 = rep("", length(samples_mt)),
                            CATEGORY_CD = rep("Biomarker_Data+Genomics+Methylation+CpG_Methylation_Percentage", length(samples_mt)),
                            SOURCE_CD= rep("STD", length(samples_mt)))
head(SubjSampMap_mt)


### Batches file
batches_mt = data.frame(patient_id = gsub("\\-\\d+\\-\\w+\\-\\d+\\-\\w+$", "", mt_batches$SampleID),
                        sample_id = mt_batches$SampleID,
                        Methylation_lanes_within_batch = paste0("lane",mt_batches$Lane),
                        stringsAsFactors = FALSE)
head(batches_mt)


### Save files
# Parental tumors
write.table(mt_data, file = "./readyForUpload/parentalTumors/mRNA_MT/raw_gex_matrix_mRNA_R.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(annotations_mt, file = "./readyForUpload/parentalTumors/mRNA_MT/MethylationGTdata.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(SubjSampMap_mt, file = "./readyForUpload/parentalTumors/mRNA_MT/Subject_Sample_Mapping.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = ".")
write.table(batches_mt, file = "./readyForUpload/parentalTumors/transmart/data_MethylationBatches.txt", 
            sep = '\t', col.names = T, row.names = F, quote = FALSE, na = "")

