#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
#options(warn = -1)
suppressMessages(require("dplyr"))
suppressMessages(require("doBy"))
suppressMessages(require("data.table"))
args <- commandArgs(TRUE)
strainName <- args[1]

#library(tictoc)
################################################################
# This script averages the depth for species, chromosomes, and windows with the -byBP option
#
#Requirements: dplyr
#Input: Window averaged depth text file
#
################################################################

# docker vars
workingDir <- "/tmp/sppIDer/working/"

spcAvgFile <- paste(workingDir, strainName, "_speciesAvgDepth-d.txt", sep="")
chrAvgFile <- paste(workingDir, strainName, "_chrAvgDepth-d.txt", sep="")
outputFileName <- paste(workingDir, strainName, "_winAvgDepth-d.txt", sep="")
dataFileName <- paste(workingDir, strainName, "-d.bedgraph", sep="")

chrLens <- read.table(paste(workingDir, strainName, "_chrLens.txt", sep=""), header=F, col.names=c("chrom", "length"))
genomeEnd <- sum(as.numeric(chrLens$length))
stepSize <- signif(genomeEnd%/%10000, digits = 2)
print(paste("Step size:", toString(stepSize), sep=" "))
valueCol <- fread(dataFileName, select=3, col.name="value")
totalMean <- mean(valueCol$value)
maxValue <- max(valueCol$value)
medianValue <- median(valueCol$value)
rm(valueCol)

split <- strsplit(as.character(chrLens$chrom), "-")
chrLens$species <- unlist(lapply(split, function(x) x[1]))
chrLens$chr <- as.numeric(unlist(lapply(split, function(x) x[2])))
uniSpecies <- unique(chrLens$species)

speciesSummary <-data.frame(Genome_Pos=numeric(), species=character(), genomeLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
spcAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "species" = "all", "genomeLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(spcAvgOutput, scientific=FALSE), file=spcAvgFile, row.names=F, sep = "\t", quote = FALSE)
spcCumLen <- 0
i=1
chrSummary <-data.frame(Genome_Pos=numeric(), chrom=character(), chrLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
chrAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "chrom" = "all", "chrLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(chrAvgOutput, scientific=FALSE), file=chrAvgFile, row.names=F, sep = "\t", quote = FALSE)
summaryInfo <- data.frame("Genome_Pos"="wholeGenome","species"="all", "chrom" = "all", "winStart"=0, "winEnd"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(summaryInfo, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t", quote = FALSE) 
winSummary <-data.frame(Genome_Pos=numeric(), species=character(), chrom=character(), winStart=numeric(), winEnd=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(), stringsAsFactors = FALSE)
chrCumLen <- 0
tic("Whole loop")
for (species in uniSpecies) {
  spcChrLens <- chrLens[which(chrLens$species==species),]
  speciesData <- fread(paste("grep", species, dataFileName, sep=" "), col.names = c("chrom", "chromPos", "value"))
  print(species)
  dataUniChr <- unique(speciesData$chrom)
  missingChr <- setdiff(spcChrLens$chrom, dataUniChr)
  tic("Add missing data")
  if (length(missingChr)>0){
    for(chrName in missingChr){
      chrLen <- spcChrLens$length[which(spcChrLens$chrom==chrName)]
      chrDF <- data.frame(chrom=rep(chrName, chrLen), chromPos=seq(1,chrLen, by=1), value=rep(0, chrLen))
      speciesData <- rbind(speciesData, chrDF)
    }
  }
  toc()
  tic("Split")
  split <- strsplit(as.character(speciesData$chrom), "-")
  speciesData$species <- unlist(lapply(split, function(x) x[1]))
  speciesData$chr <- as.numeric(unlist(lapply(split, function(x) x[2])))
  toc()
  #tic("Species Summary")
  speciesDataOrdered <- speciesData[order(speciesData$chr),]
  spcLen <- sum(spcChrLens$length)
  speciesSummary[i,1] <- spcCumLen
  speciesSummary[i,2] <- species
  speciesSummary[i,3] <- spcLen
  speciesSummary[i,4] <- mean(speciesData$value)
  speciesSummary[i,6] <- max(speciesData$value, na.rm=T)
  log2 <- log2(mean(speciesData$value/spcLen)/totalMean)
  if (log2<0) {
    log2 <- 0
  } else if (is.infinite(log2)) {
    log2 <- 0
  }
  speciesSummary[i,5] <- log2
  speciesSummary[i,7] <- median(speciesData$value)
  spcCumLen = spcCumLen + spcLen
  i=i+1
  #toc()
  j=1
  m=1
  chrData <- splitBy("chr", speciesDataOrdered)
  for (k in 1:length(chrData)){
    #tic("Chr Summary")
    #print(as.character(chrData[[k]]$chrom[1]))
    chrLen <- length(chrData[[k]]$chrom)
    chrSummary[j,1] <- chrCumLen
    chrSummary[j,2] <- as.character(chrData[[k]]$chrom[1])
    chrSummary[j,3] <- chrLen
    chrSummary[j,4] <- mean(chrData[[k]]$value)
    chrSummary[j,6] <- max(chrData[[k]]$value, na.rm=T)
    log2 <- log2(mean(chrData[[k]]$value)/totalMean)
    if (log2<0) {
      log2 <- 0
    } else if (is.infinite(log2)) {
      log2 <- 0
    }
    chrSummary[j,5] <- log2
    chrSummary[j,7] <- median(chrData[[k]]$value)
    #toc()
    oneChrData <- chrData[[k]]
    oneChrData$Genome_Pos <- oneChrData$chromPos + chrCumLen
    windowBounds <- seq(0, oneChrData$chromPos[nrow(oneChrData)], stepSize)
    windowBounds <- c(windowBounds, max(oneChrData$chromPos))
    oneChrData$bin <- cut(oneChrData$chromPos, breaks=windowBounds, include.lowest = T)
    #tic("Windows")
    winData <- splitBy("bin", oneChrData)
    for (n in 1:length(winData)){
      winLen <- sum(as.numeric(winData[[n]]$regionLen))
      winSummary[m,1] <- winData[[n]]$Genome_Pos[1]
      winSummary[m,2] <- winData[[n]]$species[1]
      winSummary[m,3] <- as.character(winData[[n]]$chrom[1])
      winSummary[m,4] <- winData[[n]]$chromPos[1]
      winSummary[m,5] <- winData[[n]]$chromPos[nrow(winData[[n]])]
      winSummary[m,6] <- mean(winData[[n]]$value)
      log2 <- log2(mean(winData[[n]]$value)/totalMean)
      if (log2<0) {
        log2 <- 0
      } else if (is.infinite(log2)) {
        log2 <- 0
      }
      winSummary[m,7] <- log2
      winSummary[m,8] <- max(winData[[n]]$value, na.rm=T)
      winSummary[m,9] <- median(winData[[n]]$value)
      m <- m+1
    }
    #toc()
    chrCumLen = chrCumLen + chrLen
    j = j+1
  }
  write.table(format(chrSummary, scientific=FALSE), file=chrAvgFile, col.names = F, row.names=F, sep = "\t", quote = FALSE, append=T)
  write.table(format(winSummary, scientific=FALSE), file=outputFileName, col.names = F, row.names=F, sep = "\t", quote = FALSE, append=T) 
}
write.table(format(speciesSummary, scientific=FALSE), file=spcAvgFile, col.names = F, row.names=F, sep = "\t", quote = FALSE, append=T)
toc()

# #Read in data
# strain <- read.table(paste(strainName, "-d.bedgraph", sep=""), header=FALSE, col.names = c("chrom", "chromPos", "value"))
# #Set window size so that the full genome is split into 10000 chunks 
# dataGenomeLen <- length(strain$chromPos)
# chrLens <- read.table(paste(strainName, "chrLens.txt", sep="_"), header=F, col.names=c("chrom", "length"))
# genomeEnd <- sum(as.numeric(chrLens$length))
# if(genomeEnd!=dataGenomeLen){
#   dataUniChr <- unique(strain$chrom)
#   missingChr <- setdiff(chrLens$chrom, dataUniChr)
#   for(chrName in missingChr){
#     chrLen <- chrLens$length[which(chrLens$chrom==chrName)]
#     chrDF <- data.frame(chrom=rep(chrName, chrLen), chromPos=seq(1,chrLen, by=1), value=rep(0, chrLen))
#     strain <- rbind(strain, chrDF)
#   }
# }
# 
# stepSize <- signif(genomeEnd%/%10000, digits = 2)
# print(paste("Step size:", toString(stepSize), sep=" "))
# totalMean <- mean(strain$value)
# maxValue <- max(strain$value)
# medianValue <- median(strain$value)
# 
# spcAvgFile <- paste(strainName, "speciesAvgDepth-d.txt", sep="_")
# chrAvgFile <- paste(strainName, "chrAvgDepth-d.txt", sep="_")
# outputFileName <- paste(strainName, "winAvgDepth-d.txt", sep="_")
# 
# #completeChromList <- unlist(list(unique(strain$chrom)))
# split <- strsplit(as.character(strain$chrom), "-")
# strain$species <- unlist(lapply(split, function(x) x[1]))
# strain$chr <- as.numeric(unlist(lapply(split, function(x) x[2])))
# 
# speciesDataOrdered <- list()
# speciesData <- splitBy("species", strain)
# speciesSummary <-data.frame(Genome_Pos=numeric(), species=character(), genomeLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
# spcCumLen <- 0
# for (i in 1:length(speciesData)){
#   spcLen <- length(speciesData[[i]]$chrom)
#   speciesSummary[i,1] <- spcCumLen
#   speciesSummary[i,2] <- as.character(names(speciesData)[[i]])
#   speciesSummary[i,3] <- spcLen
#   speciesSummary[i,4] <- mean(speciesData[[i]]$value)
#   speciesSummary[i,6] <- max(speciesData[[i]]$value, na.rm=T)
#   speciesSummary[i,5] <- log2(mean(speciesData[[i]]$value/spcLen)/totalMean)
#   speciesSummary[i,7] <- median(speciesData[[i]]$value)
#   spcCumLen = spcCumLen + spcLen
#   speciesDataOrdered[[i]] <- speciesData[[i]][order(speciesData[[i]]$chr),]
# }
# spcAvgOutput <- data.frame("Genome_Pos"=0, "species" = "all", "genomeLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
# spcAvgOutput <- rbind(spcAvgOutput, speciesSummary)
# spcAvgOutput[which(is.infinite(spcAvgOutput$log2mean)==T),"log2mean"]=0
# spcAvgOutput[which(spcAvgOutput$log2mean<0),"log2mean"]=0
# spcAvgOutput[1,1] <- "wholeGenome"
# write.table(format(spcAvgOutput, scientific=FALSE), file=spcAvgFile, row.names=F, sep = "\t", quote = FALSE) 
# 
# chrData <- lapply(speciesDataOrdered, function(x) splitBy("chr", x))
# chrSummary <-data.frame(Genome_Pos=numeric(), chrom=character(), chrLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
# winSummary <-data.frame(Genome_Pos=numeric(), species=character(), chrom=character(), winStart=numeric(), winEnd=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(), stringsAsFactors = FALSE)
# chrCumLen <- 0
# j=1
# m=1
# for (i in 1:length(chrData)){
#   for (k in 1:length(chrData[[i]])){
#     chrLen <- length(chrData[[i]][[k]]$chrom)
#     chrSummary[j,1] <- chrCumLen
#     chrSummary[j,2] <- as.character(chrData[[i]][[k]]$chrom[1])
#     chrSummary[j,3] <- chrLen
#     chrSummary[j,4] <- mean(chrData[[i]][[k]]$value)
#     chrSummary[j,6] <- max(chrData[[i]][[k]]$value, na.rm=T)
#     chrSummary[j,5] <- log2(mean(chrData[[i]][[k]]$value)/totalMean)
#     chrSummary[j,7] <- median(chrData[[i]][[k]]$value)
#     oneChrData <- chrData[[i]][[k]]
#     oneChrData$Genome_Pos <- oneChrData$chromPos + chrCumLen
#     windowBounds <- seq(0, oneChrData$chromPos[nrow(oneChrData)], stepSize)
#     windowBounds <- c(windowBounds, max(oneChrData$chromPos))
#     oneChrData$bin <- cut(oneChrData$chromPos, breaks=windowBounds, include.lowest = T)
#     winData <- splitBy("bin", oneChrData)
#     for (n in 1:length(winData)){
#       winLen <- sum(as.numeric(winData[[n]]$regionLen))
#       winSummary[m,1] <- winData[[n]]$Genome_Pos[1]
#       winSummary[m,2] <- winData[[n]]$species[1]
#       winSummary[m,3] <- as.character(winData[[n]]$chrom[1])
#       winSummary[m,4] <- winData[[n]]$chromPos[1]
#       winSummary[m,5] <- winData[[n]]$chromPos[nrow(winData[[n]])]
#       winSummary[m,6] <- mean(winData[[n]]$value)
#       winSummary[m,7] <- log2(mean(winData[[n]]$value)/totalMean)
#       winSummary[m,8] <- max(winData[[n]]$value, na.rm=T)
#       winSummary[m,9] <- median(winData[[n]]$value)
#       m <- m+1
#     }
#     chrCumLen = chrCumLen + chrLen
#     j = j+1
#   }
# }
# chrAvgOutput <- data.frame("Genome_Pos"=0, "chrom" = "all", "chrLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
# chrAvgOutput <- rbind(chrAvgOutput, chrSummary)
# chrAvgOutput[which(is.infinite(chrAvgOutput$log2mean)==T),"log2mean"]=0
# chrAvgOutput[which(chrAvgOutput$log2mean<0),"log2mean"]=0
# chrAvgOutput[1,1] <- "wholeGenome"
# write.table(format(chrAvgOutput, scientific=FALSE), file=chrAvgFile, row.names=F, sep = "\t", quote = FALSE)
# summaryInfo <- data.frame("Genome_Pos"=0,"species"="all", "chrom" = "all", "winStart"=0, "winEnd"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
# summaryInfo <- rbind(summaryInfo, winSummary)
# summaryInfo[which(is.infinite(summaryInfo$log2mean)==T),"log2mean"]=0
# summaryInfo[which(summaryInfo$log2mean<0),"log2mean"]=0
# summaryInfo[1,1] <- "wholeGenome"
# write.table(format(summaryInfo, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t", quote = FALSE) 
