#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
#options(warn = -1)
suppressMessages(require("dplyr"))
suppressMessages(require("doBy"))
suppressMessages(require("data.table"))
args <- commandArgs(TRUE)
strainName <- args[1]

## INLINE TESTS
#print(strainName)
#print("Grep test")

#suppressMessages(require("tictoc"))
################################################################
# This script averages the depth for species, chromosomes, and windows with the -byBP option
#
#Requirements: dplyr
#Input: Window averaged depth text file
#
################################################################

# docker vars
#workingDir <- "/tmp/sppIDer/working/"
workingDir <- ""

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
#i=1
chrAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "chrom" = "all", "chrLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(chrAvgOutput, scientific=FALSE), file=chrAvgFile, row.names=F, sep = "\t", quote = FALSE)
summaryInfo <- data.frame("Genome_Pos"="wholeGenome","species"="all", "chrom" = "all", "winStart"=0, "winEnd"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(summaryInfo, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t", quote = FALSE) 
chrCumLen <- 0
#tic("Whole loop")
for (species in uniSpecies) {
  spcChrLens <- chrLens[which(chrLens$species==species),]

  # pipe the result of "awk {'print $1, $2, $3'} dataFileName | grep species", with supplied column names, into speciesData via read.table
  speciesData <- read.table(pipe(paste("awk {'print $1, $2, $3'}", dataFileName, "| grep", species, sep=" ")), col.names=c("chrom", "chromPos", "value"))

  #print(species)
  dataUniChr <- unique(speciesData$chrom)
  missingChr <- setdiff(spcChrLens$chrom, dataUniChr)
  #tic("Add missing data")
  if (length(missingChr)>0){
    for(chrName in missingChr){
      chrLen <- spcChrLens$length[which(spcChrLens$chrom==chrName)]
      chrDF <- data.frame(chrom=rep(chrName, chrLen), chromPos=seq(1,chrLen, by=1), value=rep(0, chrLen))
      speciesData <- rbind(speciesData, chrDF)
    }
  }
  #toc()
  #tic("Split")
  split <- strsplit(as.character(speciesData$chrom), "-")
  #speciesData$species <- unlist(lapply(split, function(x) x[1]))
  speciesData$chr <- as.numeric(unlist(lapply(split, function(x) x[2])))
  #toc()
  #tic("Species Summary")
  spcLen <- sum(spcChrLens$length)
  speciesSummary[1,1] <- spcCumLen
  speciesSummary[1,2] <- species
  speciesSummary[1,3] <- spcLen
  speciesSummary[1,4] <- mean(speciesData$value)
  speciesSummary[1,6] <- max(speciesData$value, na.rm=T)
  log2 <- log2(mean(speciesData$value)/totalMean)
  if (log2<0) {
    log2 <- 0
  } else if (is.infinite(log2)) {
    log2 <- 0
  }
  speciesSummary[1,5] <- log2
  speciesSummary[1,7] <- median(speciesData$value)
  spcCumLen = spcCumLen + spcLen
  write.table(format(speciesSummary, scientific=FALSE), file=spcAvgFile, col.names = F, row.names=F, sep = "\t", quote = FALSE, append=T)
  #i=i+1
  #toc()
  chrSummary <-data.frame(Genome_Pos=numeric(), chrom=character(), chrLen=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(),stringsAsFactors = FALSE)
  winSummary <-data.frame(Genome_Pos=numeric(), species=character(), chrom=character(), winStart=numeric(), winEnd=numeric(), meanValue=numeric(), log2mean=numeric(), max=numeric(), median=numeric(), stringsAsFactors = FALSE)
  #j=1
  #m=1
  speciesDataOrdered <- speciesData[order(speciesData$chr),]
  rm(speciesData)
  chrData <- splitBy("chr", speciesDataOrdered)
  for (k in 1:length(chrData)){
    #tic("Chr Summary")
    #print(as.character(chrData[[k]]$chrom[1]))
    chrLen <- length(chrData[[k]]$chrom)
    chrSummary[1,1] <- chrCumLen
    chrSummary[1,2] <- as.character(chrData[[k]]$chrom[1])
    chrSummary[1,3] <- chrLen
    chrSummary[1,4] <- mean(chrData[[k]]$value)
    chrSummary[1,6] <- max(chrData[[k]]$value, na.rm=T)
    log2 <- log2(mean(chrData[[k]]$value)/totalMean)
    if (log2<0) {
      log2 <- 0
    } else if (is.infinite(log2)) {
      log2 <- 0
    }
    chrSummary[1,5] <- log2
    chrSummary[1,7] <- median(chrData[[k]]$value)
    write.table(format(chrSummary, scientific=FALSE), file=chrAvgFile, col.names = F, row.names=F, sep = "\t", quote = FALSE, append=T)
    #toc()
    oneChrData <- chrData[[k]]
    oneChrData$Genome_Pos <- oneChrData$chromPos + chrCumLen
    windowBounds <- seq(0, oneChrData$chromPos[nrow(oneChrData)], stepSize)
    if (max(oneChrData$chromPos) %in% windowBounds == F) {
      windowBounds <- c(windowBounds, max(oneChrData$chromPos))
    }
    oneChrData$bin <- cut(oneChrData$chromPos, breaks=windowBounds, include.lowest = T)
    #tic("Windows")
    winData <- splitBy("bin", oneChrData)
    for (n in 1:length(winData)){
      winLen <- sum(as.numeric(winData[[n]]$regionLen))
      winSummary[1,1] <- winData[[n]]$Genome_Pos[1]
      winSummary[1,2] <- species
      winSummary[1,3] <- as.character(winData[[n]]$chrom[1])
      winSummary[1,4] <- winData[[n]]$chromPos[1]
      winSummary[1,5] <- winData[[n]]$chromPos[nrow(winData[[n]])]
      winSummary[1,6] <- mean(winData[[n]]$value)
      log2 <- log2(mean(winData[[n]]$value)/totalMean)
      if (log2<0) {
        log2 <- 0
      } else if (is.infinite(log2)) {
        log2 <- 0
      }
      winSummary[1,7] <- log2
      winSummary[1,8] <- max(winData[[n]]$value, na.rm=T)
      winSummary[1,9] <- median(winData[[n]]$value)
      write.table(format(winSummary, scientific=FALSE), file=outputFileName, col.names = F, row.names=F, sep = "\t", quote = FALSE, append=T)
      #m <- m+1
    }
    #toc()
    chrCumLen = chrCumLen + chrLen
    #j = j+1
  }
  rm(chrData)
  rm(oneChrData)
  rm(speciesDataOrdered)
}
#toc()