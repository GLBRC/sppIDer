#options(warn = -1)
suppressMessages(require("dplyr"))
args <- commandArgs(TRUE)
strainName <- args[1]

################################################################
# This script averages the depth for species, chromosomes, and windows with the -byBP option
#
#Requirements: dplyr
#Input: Window averaged depth text file
#
################################################################


#Read in data
strain <- read.table(paste(strainName, "-d.bedgraph", sep=""), header=FALSE, col.names = c("chrom", "chromPos", "value"))
#Set window size so that the full genome is split into 10000 chunks 
stepSize <- signif(length(strain$chrom)%/%10000, digits = 2)
print(paste("Step size:", toString(stepSize), sep=" "))
spcAvgFile <- paste(strainName, "speciesAvgDepth-d.txt", sep="_")
chrAvgFile <- paste(strainName, "chrAvgDepth-d.txt", sep="_")

completeChromList <- unlist(list(unique(strain$chrom)))
species <- c()
#Split chr name from spcecies name
for (i in 1:length(completeChromList)) {
  completeChrString <- toString(completeChromList[[i]])
  split <- unlist(strsplit(completeChrString, "-"))
  species <- c(species, split[1])
}
uniSpecies <- unique(species)

genomeEnd <- strain$chromPos[length(strain$chromPos)]
totalMean <- mean(strain$value)
maxValue <- max(strain$value)
medianValue <- median(strain$value)
spcAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "species" = "all", "genomeLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(spcAvgOutput, scientific=FALSE), file=spcAvgFile, row.names=F, sep = "\t", quote = FALSE) 
spcCumLen <- 0
options(warn = -1)
#Go though by species and print summary data about depth of coverage to a file
for (k in 1:length(uniSpecies)){
  spcName <- toString(uniSpecies[[k]])
  spc <- filter(strain, grepl(spcName, strain$chrom))
  spcLength <- length(spc$chrom)
  spcMean <- mean(spc$value)
  spcMedian <- median(spc$value)
  spcAvg <- data.frame("Genome_Pos"=spcCumLen, "species" = spcName, "genomeLen"=spcLength, "meanValue" =spcMean, "log2mean"=log2(spcMean/totalMean), "max"=max(spc$value), "median"=spcMedian)
  write.table(format(spcAvg, scientific=FALSE), file=spcAvgFile, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE) 
  spcAvgOutput <- rbind(spcAvgOutput, spcAvg)
  spcCumLen <- spcCumLen + spcLength
}
#Go though by chromosome and print summary data about depth of coverage to a file
chrAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "chrom" = "all", "chrLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(chrAvgOutput, scientific=FALSE), file=chrAvgFile, row.names=F, sep = "\t", quote = FALSE)
chrCumLen <- 0
wGenomeLen <- data.frame()
for (i in 1:length(completeChromList)) {
  chrName <- toString(completeChromList[[i]])
  chr <- subset(strain, chrom==chrName)
  chrLength <- length(chr$chrom)
  chrMean <- mean(chr$value)
  chrMedian <- median(chr$value)
  chrAvg <- data.frame("Genome_Pos"=chrCumLen, "chrom" = chrName, "chrLen"=chrLength, "meanValue" = chrMean, "log2mean"=log2(chrMean/totalMean), "max"=max(chr$value), "median"=chrMedian)
  write.table(format(chrAvg, scientific=FALSE), file=chrAvgFile, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE)
  chrCumLen <- chrCumLen + chrLength
}
#Go though by windows and print summary data about depth of coverage to a file
outputFileName <- paste(strainName, "winAvgDepth-d.txt", sep="_")
output <- data.frame()
totalMean <- mean(strain$value)
summaryInfo <- data.frame("Genome_Pos"="wholeGenome","species"="all", "chrom" = "all", "winStart"=0, "winEnd"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(summaryInfo, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t", quote = FALSE) 
spcCounter = 0
for (i in 1:length(completeChromList)) {
  chr_name <- toString(completeChromList[[i]])
  split <- unlist(strsplit(chr_name, "-"))
  speciesName <- split[1]
  chr <- subset(strain, chrom==chr_name)
  chr.out <- data.frame()
  counter = 0
  start = 0
  if (0 < length(chr[,1])) {
    while (start < length(chr[,1])) {
      window.out <- data.frame(Genome_Pos=spcCounter, species=NA, chrom = NA, chromstart = NA, chromend = NA, meanValue = NA, log2mean=NA, max=NA, median=NA)
      window <- subset(chr, (chr$chromPos>start & chr$chromPos<=counter+stepSize))
      window.out$species <- speciesName
      window.out$chrom <- chr_name  
      window.out$chromstart <- counter
      window.out$chromend <- counter + stepSize
      if (!is.na(window[1,2])) {
        window.out$meanValue <-  mean(window$value)
        window.out$log2mean <- log2((mean(window$value))/totalMean)
        window.out$max <- max(window$value)
        window.out$median <- median(window$value)
        start <- start + stepSize
      } else if (is.na(window[1,2])) {
        window.out$meanValue <-  0
        window.out$log2mean <- 0
        start <- counter + stepSize
      }
      chr.out <- rbind(chr.out, window.out)   
      counter = counter + stepSize  
      write.table(format(window.out, scientific=FALSE), file=outputFileName, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE)
      spcCounter <- spcCounter+stepSize
    }
  } else {
    window.out <- data.frame(chrom = NA, start = NA, end = NA, meanValue = NA)
    window <- subset(chr, (chr$chromPos>start & chr$chromPos<=counter+stepSize))
    window.out$chrom <- chr_name
    window.out$start <- 0
    window.out$end <- 0 + stepSize
    window.out$meanValue <-  0
    window.out$log2mean <- 0
    window.out$max <- 0
    window.out$median <- 0
    chr.out <- rbind(chr.out, window.out) 
    write.table(format(window.out, scientific=FALSE), file=outputFileName, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE)
    spcCounter <- spcCounter+stepSize
  }
  output <- rbind(output, chr.out)
}
