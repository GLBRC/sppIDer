#options(warn = -1)
suppressMessages(require("dplyr"))
args <- commandArgs(TRUE)
strainName <- args[1]

################################################################
# This script averages the depth for species, chromosomes, and windows with the -byGroup option
#
#Requirements: dplyr
#Input: Window averaged depth text file
#
################################################################


#Read in data
strain <- read.table(paste(strainName, ".bedgraph", sep=""), header=FALSE, col.names = c("chrom", "regionStart", "regionEnd",  "value"))
strain$regionLen <- strain$regionEnd-strain$regionStart
strain$wtValue <- strain$value*strain$regionLen
genomeEnd <- sum(as.numeric(strain$regionLen))
#Set window size so that the full genome is split into 10000 chunks 
stepSize <- signif(genomeEnd%/%10000, digits = 2)
print(paste("Step size:", toString(stepSize), sep=" "))
totalMean <- sum(as.numeric(strain$wtValue))/genomeEnd
medianValue=median(rep(strain$value, strain$regionLen))
maxValue=max(strain$value)
spcAvgFile <- paste(strainName, "speciesAvgDepth-g.txt", sep="_")
chrAvgFile <- paste(strainName, "chrAvgDepth-g.txt", sep="_")

completeChromList <- unlist(list(unique(strain$chrom)))
species <- c()
#Split chr name from spcecies name
for (i in 1:length(completeChromList)) {
  completeChrString <- toString(completeChromList[[i]])
  split <- unlist(strsplit(completeChrString, "-"))
  species <- c(species, split[1])
}
uniSpecies <- unique(species)

spcAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "species" = "all", "genomeLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(spcAvgOutput, scientific=FALSE), file=spcAvgFile, row.names=F, sep = "\t", quote = FALSE) 
spcCumLen <- 0
options(warn = -1)
#Go though by species and print summary data about depth of coverage to a file
for (k in 1:length(uniSpecies)){
  spcName <- toString(uniSpecies[[k]])
  spc <- filter(strain, grepl(spcName, strain$chrom))
  spcLength <- sum(as.numeric(spc$regionLen))
  spcMean <- sum(as.numeric(spc$wtValue))/spcLength
  spcMedian <- median(rep(spc$value, spc$regionLen))
  spcAvg <- data.frame("Genome_Pos"=spcCumLen, "species" = spcName, "genomeLen"=spcLength, "meanValue" =spcMean, "log2mean"=log2(spcMean/totalMean), "max"=max(spc$value), "median"=spcMedian)
  write.table(format(spcAvg, scientific=FALSE), file=spcAvgFile, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE) 
  spcAvgOutput <- rbind(spcAvgOutput, spcAvg)
  spcCumLen <- spcCumLen + spcLength
}

chrAvgOutput <- data.frame("Genome_Pos"="wholeGenome", "chrom" = "all", "chrLen"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(chrAvgOutput, scientific=FALSE), file=chrAvgFile, row.names=F, sep = "\t", quote = FALSE)
chrCumLen <- 0
wGenomeLen <- data.frame()
#Go though by chromosome and print summary data about depth of coverage to a file
for (i in 1:length(completeChromList)) {
  chrName <- toString(completeChromList[[i]])
  chr <- subset(strain, chrom==chrName)
  chrLength <- sum(as.numeric(chr$regionLen))
  chrMean <- sum(as.numeric(chr$wtValue))/chrLength
  chrMedian <- median(rep(chr$value, chr$regionLen))
  chrAvg <- data.frame("Genome_Pos"=chrCumLen, "chrom" = chrName, "chrLen"=chrLength, "meanValue" = chrMean, "log2mean"=log2(chrMean/totalMean), "max"=max(chr$value), "median"=chrMedian)
  write.table(format(chrAvg, scientific=FALSE), file=chrAvgFile, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE)
  chrAvgOutput <- rbind(chrAvgOutput, chrAvg)
  chrWgenomeLen <- chr
  chrWgenomeLen$Genome_Pos <- chr$regionStart + chrCumLen
  wGenomeLen <- rbind(wGenomeLen, chrWgenomeLen)
  chrCumLen <- chrCumLen + chrLength
}

outputFileName <- paste(strainName, "winAvgDepth-g.txt", sep="_")
summaryInfo <- data.frame("Genome_Pos"="wholeGenome","species"="all", "chrom" = "all", "winStart"=0, "winEnd"=genomeEnd, "meanValue" = totalMean, "log2mean"=1, "max"=maxValue, "median"=medianValue)
write.table(format(summaryInfo, scientific=FALSE), file=outputFileName, row.names=F, sep = "\t", quote = FALSE) 
win.output <- data.frame()
counter = 0
#Go though by windows and print summary data about depth of coverage to a file
while (counter < genomeEnd) {
    window <- subset(wGenomeLen, (wGenomeLen$Genome_Pos>=counter & wGenomeLen$Genome_Pos<=counter+stepSize))
    if (length(window$chrom)>0) {
      window.out <- data.frame(Genome_Pos=counter, species=NA, chrom = NA, chromstart = NA, chromend = NA, meanValue = NA, log2mean=NA, max=NA, median=NA)
      chrName <- window$chrom[1]
      split <- unlist(strsplit(toString(chrName), "-"))
      speciesName <- split[1]
      window.out$species <- speciesName
      window.out$chrom <- chrName 
      if (!is.na(window[1,2])) {
        window.out$chromstart <- window$regionStart[1]
        window.out$chromend <- window$regionEnd[length(window$regionEnd)]
        window.out$meanValue <-  sum(as.numeric(window$wtValue))/stepSize
        window.out$log2mean <- log2((sum(as.numeric(window$wtValue))/stepSize)/totalMean)
        window.out$median <- median(rep(window$value, window$regionLen))
        window.out$max <- max(window$value)
      } else if (is.na(window[1,2])) {
        window.out$chromstart <- newStart
        window.out$chromend <- counter + stepSize
        window.out$meanValue <-  0
        window.out$log2mean <- 0
        window.out$median <- 0
        window.out$max <- 0
      }
      write.table(format(window.out, scientific=FALSE), file=outputFileName, append = T, col.names=F, row.names=F, sep = "\t", quote = FALSE)
      counter = counter + stepSize 
      win.output <- rbind(win.output, window.out)
    } else {
      counter = counter + stepSize 
    }
}
