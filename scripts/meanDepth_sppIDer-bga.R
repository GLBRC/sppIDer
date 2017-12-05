#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
#options(warn = -1)
suppressMessages(require("dplyr"))
suppressMessages(require("doBy"))
args <- commandArgs(TRUE)
strainName <- args[1]

################################################################
# This script averages the depth for species, chromosomes, and windows with the -byGroup option
#
#Requirements: dplyr
#Input: Bedgraph -bga output
#
################################################################

# docker vars
workingDir <- "/tmp/sppIDer/working/"

#Read in data
strain <- read.table(paste(workingDir, strainName, ".bedgraph", sep=""), header=FALSE, col.names = c("chrom", "regionStart", "regionEnd",  "value"))
strain$regionLen <- strain$regionEnd-strain$regionStart
dataGenomeLen <- sum(as.numeric(strain$regionLen))
chrLens <- read.table(paste(workingDir, strainName, "_chrLens.txt", sep=""), header=F, col.names=c("chrom", "length"))
genomeEnd <- sum(as.numeric(chrLens$length))
if(genomeEnd!=dataGenomeLen){
  dataUniChr <- unique(strain$chrom)
  missingChr <- setdiff(chrLens$chrom, dataUniChr)
  for(chrName in missingChr){
    chrLen <- chrLens$length[which(chrLens$chrom==chrName)]
    chrDF <- data.frame(chrom=chrName, regionStart=0, regionEnd=chrLen, value=0, regionLen=chrLen)
    strain <- rbind(strain, chrDF)
  }
}

strain$wtValue <- strain$value*strain$regionLen
#Set window size so that the full genome is split into 10000 chunks 
stepSize <- signif(genomeEnd%/%10000, digits = 2)
windowBounds <- seq(0, genomeEnd, by=stepSize)
print(paste("Step size:", toString(stepSize), sep=" "))
totalMean <- sum(as.numeric(strain$wtValue))/genomeEnd
medianValue<-median(rep(strain$value, strain$regionLen))
maxValue<-max(strain$value)

spcAvgFile <- paste(workingDir, strainName, "_speciesAvgDepth-g.txt", sep="")
chrAvgFile <- paste(workingDir, strainName, "_chrAvgDepth-g.txt", sep="")
outputFileName <- paste(workingDir, strainName, "_winAvgDepth-g.txt", sep="")

#completeChromList <- unlist(list(unique(strain$chrom)))
split <- strsplit(as.character(strain$chrom), "-")
strain$species <- unlist(lapply(split, function(x) x[1]))
strain$chr <- as.numeric(unlist(lapply(split, function(x) x[2])))

speciesDataOrdered <- list()
speciesData <- splitBy("species", strain)
speciesSummary <-data.frame(Genome_Pos=numeric(), species=character(), genomeLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
spcCumLen <- 0
for (i in 1:length(speciesData)){
  spcLen <- sum(as.numeric(speciesData[[i]]$regionLen))
  speciesSummary[i,1] <- spcCumLen
  speciesSummary[i,2] <- as.character(names(speciesData)[[i]])
  speciesSummary[i,3] <- spcLen
  speciesSummary[i,4] <- sum(as.numeric(speciesData[[i]]$wtValue))/spcLen
  speciesSummary[i,6] <- max(speciesData[[i]]$value, na.rm=T)
  speciesSummary[i,5] <- log2(sum(as.numeric(speciesData[[i]]$wtValue)/spcLen)/totalMean)
  speciesSummary[i,7] <- median(rep(speciesData[[i]]$value, speciesData[[i]]$regionLen))
  spcCumLen = spcCumLen + spcLen
  speciesDataOrdered[[i]] <- speciesData[[i]][order(speciesData[[i]]$chr),]
}
spcAvgOutput <- data.frame("Genome_Pos"=0, "species" = "all", "genomeLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
spcAvgOutput <- rbind(spcAvgOutput, speciesSummary)
spcAvgOutput[which(is.infinite(spcAvgOutput$log2mean)==T),"log2mean"]=0
spcAvgOutput[which(spcAvgOutput$log2mean<0),"log2mean"]=0
spcAvgOutput[1,1] <- "wholeGenome"
write.table(format(spcAvgOutput, scientific=FALSE), file=spcAvgFile, row.names=F, sep = "\t", quote = FALSE) 

chrData <- lapply(speciesDataOrdered, function(x) splitBy("chr", x))
chrSummary <-data.frame(Genome_Pos=numeric(), chrom=character(), chrLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
winSummary <-data.frame(Genome_Pos=numeric(), species=character(), chrom=character(), winStart=numeric(), winEnd=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(), stringsAsFactors = FALSE)
chrCumLen <- 0
j=1
m=1
for (i in 1:length(chrData)){
  for (k in 1:length(chrData[[i]])){
    #print(j)
    chrLen <- sum(as.numeric(chrData[[i]][[k]]$regionLen))
    chrSummary[j,1] <- chrCumLen
    chrSummary[j,2] <- as.character(chrData[[i]][[k]]$chrom[1])
    chrSummary[j,3] <- chrLen
    chrSummary[j,4] <- sum(as.numeric(chrData[[i]][[k]]$wtValue))/chrLen
    chrSummary[j,6] <- max(chrData[[i]][[k]]$value, na.rm=T)
    chrSummary[j,5] <- log2(sum(as.numeric(chrData[[i]][[k]]$wtValue)/chrLen)/totalMean)
    chrSummary[j,7] <- median(rep(chrData[[i]][[k]]$value, chrData[[i]][[k]]$regionLen))
    oneChrData <- chrData[[i]][[k]]
    oneChrData$Genome_Pos <- oneChrData$regionStart + chrCumLen
    windowBounds <- seq(0, oneChrData$regionEnd[nrow(oneChrData)], stepSize)
    windowBounds <- c(windowBounds, max(oneChrData$regionEnd))
    oneChrData$bin <- cut(oneChrData$regionStart, breaks=windowBounds, include.lowest = T)
    winData <- splitBy("bin", oneChrData)
    for (n in 1:length(winData)){
      winLen <- sum(as.numeric(winData[[n]]$regionLen))
      winSummary[m,1] <- winData[[n]]$Genome_Pos[1]
      winSummary[m,2] <- winData[[n]]$species[1]
      winSummary[m,3] <- as.character(winData[[n]]$chrom[1])
      winSummary[m,4] <- winData[[n]]$regionStart[1]
      winSummary[m,5] <- winData[[n]]$regionEnd[length(winData[[n]]$regionEnd)]
      winSummary[m,6] <- sum(as.numeric(winData[[n]]$wtValue))/winLen
      winSummary[m,7] <- log2(sum(as.numeric(winData[[n]]$wtValue)/winLen)/totalMean)
      winSummary[m,8] <- max(winData[[n]]$value, na.rm=T)
      winSummary[m,9] <- median(rep(winData[[n]]$value, winData[[n]]$regionLen))
      m <- m+1
    }
    chrCumLen = chrCumLen + chrLen
    j = j+1
  }
}

chrAvgOutput <- data.frame("Genome_Pos"=0, "chrom" = "all", "chrLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
chrAvgOutput <- rbind(chrAvgOutput, chrSummary)
chrAvgOutput[which(is.infinite(chrAvgOutput$log2mean)==T),"log2mean"]=0
chrAvgOutput[which(chrAvgOutput$log2mean<0),"log2mean"]=0
chrAvgOutput[1,1] <- "wholeGenome"
write.table(format(chrAvgOutput, scientific=FALSE), file=chrAvgFile, row.names=F, sep = "\t", quote = FALSE)
summaryInfo <- data.frame("Genome_Pos"=0,"species"="all", "chrom" = "all", "winStart"=0, "winEnd"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
summaryInfo <- rbind(summaryInfo, winSummary)
summaryInfo[which(is.infinite(summaryInfo$log2mean)==T),"log2mean"]=0
summaryInfo[which(summaryInfo$log2mean<0),"log2mean"]=0
summaryInfo[1,1] <- "wholeGenome"
write.table(format(summaryInfo, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t", quote = FALSE) 
