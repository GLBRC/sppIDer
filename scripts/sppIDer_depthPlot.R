#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
require(ggplot2)
require(modes)
plot <- ggplot()
args <- commandArgs(TRUE)
outputPrefix <- args[1]

################################################################
# This script makes the windowed depth plots with the -byBP option
#
#Requirements: ggplot2
#Input: Window averaged depth text file
#
################################################################

# docker vars
#workingDir <- "/tmp/sppIDer/working/"
workingDir <- ""

introgressCutoff <- 0.01
#Read in data, get info on window size and spread of mean values. Add log2 and a rescaled mean column to be plotted later.
#bedData <- read.table(paste(workingDir, outputPrefix, "_winAvgDepth-d.txt", sep=""), header=T)
#bedData <- read.table(paste(workingDir, outputPrefix, "_winAvgDepth-d.txt", sep=""), header=F, skip=2, col.names = c("Genome_Pos", "species", "chrom", "start",  "end", "meanValue", "log2mean", "max", "median"))
windowedDataFileName <- paste(workingDir, outputPrefix, "_winAvgDepth-d.txt", sep="")
if (file.exists(windowedDataFileName)){
  bedData <- read.table(windowedDataFileName, header=F, skip=2, col.names = c("Genome_Pos", "species", "chrom", "start",  "end", "meanValue", "log2rawCov", "max", "median"))
} else {
  windowedDataFileName <- paste(workingDir, outputPrefix, "_winAvgDepth-g.txt", sep="")
  bedData <- read.table(windowedDataFileName, header=F, skip=2, col.names = c("Genome_Pos", "species", "chrom", "start",  "end", "meanValue", "log2rawCov", "max", "median"))
}
completeChromList <- unlist(list(unique(bedData$chrom)))
stepSize <- bedData[2,1]
checkMean <- mean(bedData$meanValue)
#meanQuant <- quantile(bedData$meanValue, 0.99, na.rm=T)
bedData$log2 <- log2(bedData$meanValue/checkMean)
bedData$log2[bedData$log2<0] <- 0
logQuant <- quantile(bedData$log2, 0.99, na.rm=T)

#Assign colors based on number of species. Here you can use your own color scheme. 
pdf(NULL)
uniSpecies <- unlist(list(unique(bedData$species)))
uniSpecies <- factor(uniSpecies, levels=uniSpecies)
if (length(uniSpecies)>1) {
  colorList <- palette(rainbow(length(uniSpecies)))
  colorList <- palette(rainbow(length(uniSpecies)))
} else {
  colorList <- c("red")
}
spcLabels <- gsub("_", "\n", uniSpecies)
colors <- c()
for (i in uniSpecies){
  index <- which(uniSpecies %in% i)
  colors[i] <- colorList[index]
}
fillScale <- scale_fill_manual(name="Species", values = colors, breaks=uniSpecies)
colorScale <- scale_color_manual(name="Species", values = colors, breaks=uniSpecies)
theme <- theme_bw()+theme(legend.text=element_text(face="italic"))

data <- bedData$log2
ampOut <- amps(data)
#plot(ampOut$Peaks[,1]~ampOut$Peaks[,2])
antimodes <- ampOut$Antimode[,1]
binInfo <- data.frame("binName"=character(), "log2BinStart"=numeric(), "log2BinEnd"=numeric(), "numWin"=numeric(), "perData"=numeric(),stringsAsFactors=FALSE)
k=1
for (i in 1:(length(antimodes)+1)) {
  if (i==1){
    binStart<-0
  }
  if (i==length(antimodes)+1) {
    binEnd <- max(data)
  } else {
    binEnd <- antimodes[i]
  }
  numWin <- length(which(data>=binStart & data<=binEnd))
  binName <- paste("bin",(k-1),sep="")
  binInfo[k,1] <- binName
  binInfo[k,2] <- binStart
  binInfo[k,3] <- binEnd
  binInfo[k,4] <- numWin
  binInfo[k,5] <- round(100*(numWin/length(data)),2)
  k <- k+1
  binStart <- binEnd
}
#head(binInfo)
numBinsWanted <- max(which(binInfo$numWin>=500))+1
binsWanted <- binInfo[1:numBinsWanted,]
binsWanted$meanCovUpper <- (2^binsWanted$log2BinEnd)*checkMean
binInfo$meanCovUpper <- (2^binInfo$log2BinEnd)*checkMean
write.table(binsWanted, file=paste(outputPrefix, "_globalBinsIncluded.txt", sep=""), row.names = F)

binOut <- bedData
binOut$covBin <- "aboveUpperBin"
for (i in 1:nrow(binsWanted)){
  binName <- binsWanted$binName[i]
  binOut$covBin[which(binOut$log2>=binsWanted[i,2] & binOut$log2<=binsWanted[i,3])] <- binName
}
write.table(binOut, file=paste(outputPrefix, "_winAvgDepth_wCovBins.txt", sep=""), row.names = F)

limitValue <- binInfo$meanCovUpper[nrow(binsWanted)+1]
bedData$meanValueLimited <- bedData$meanValue
bedData$meanValueLimited[bedData$meanValueLimited>limitValue] <- limitValue

lowerCutoff <- binsWanted$log2BinEnd[1]
upperCutoff <- binsWanted$log2BinEnd[nrow(binsWanted)]

spSum <- data.frame(species=character(), perAboveCutoff=character(), contributes=character(), stringsAsFactors = F)
binSum <- c()
k=1
for (spName in uniSpecies){
  spData <- bedData[which(bedData$species==spName),]
  nWinPer <- nrow(spData)*introgressCutoff
  cutoffData <- spData[which(spData$log2>=lowerCutoff),]
  percentage <- paste(round((nrow(cutoffData)/nrow(spData))*100, 2), "%" , sep="")
  spSum[k,1] <- spName
  spSum[k,2] <- percentage
  spSum[k,3] <- nrow(cutoffData)>=nWinPer
  binString <- spName
  for (j in 1:nrow(binsWanted)+1){
    if (j == nrow(binsWanted)+1) {
      binData <- spData[which(spData$log2>=binsWanted$log2BinEnd[j-1]),]
      binPer <- paste(round((nrow(binData)/nrow(spData))*100, 2), "%" , sep="")
      binString <- paste(binString, binPer, sep=",")
    } else {
      binData <- spData[which(spData$log2>=binsWanted$log2BinStart[j] & spData$log2<=binsWanted$log2BinEnd[j]),]
      binPer <- paste(round((nrow(binData)/nrow(spData))*100, 2), "%" , sep="")
      binString <- paste(binString, binPer, sep=",")
    }
  }
  if (k==1) {
    binSum <- binString
  } else {
    binSum <- c(binSum, binString)
  }
  k <- k+1
}
#print(binSum)
binHeader <- "Speices"
for (i in 1:nrow(binsWanted)+1){
  if (i==nrow(binsWanted)+1) {
    binName <- "avboveUpperBin"
    binHeader <- c(binHeader, binName)
  } else {
    binName <- paste("bin", i-1, ":", round(binsWanted$meanCovUpper[i-1], 2), "X-", round(binsWanted$meanCovUpper[i], 2), "X", sep="")
    binHeader <- c(binHeader, binName)
  }
}
binTemp <- data.frame(data=binSum)
binDF <- data.frame(do.call('rbind', strsplit(as.character(binTemp$data),',',fixed=TRUE)))
colnames(binDF) <- binHeader
spSum <- cbind(spSum, binDF[,2:ncol(binDF)])
write.table(spSum, file=paste(outputPrefix, "_covBinsSummary.txt", sep=""), row.names = F)

sigSpecies <- spSum$species[which(spSum$contributes==T)]
sigSpecies <- factor(sigSpecies, levels=sigSpecies)
sigSpBins <- data.frame(species=character(), binName=character(), log2binStart=numeric(), log2binEnd=numeric())
sigData <- data.frame()
sigDataOut <- data.frame()
sigSpColNames <- colnames(binOut)
binOutTemp <- binOut
for (spc in sigSpecies) {
  sigData <- rbind(sigData, bedData[which(bedData$species==spc),])
  spcWin <- binOutTemp[which(binOutTemp$species==spc),]
  sigDataOut <- rbind(sigDataOut, spcWin)
  spLog2 <- spcWin$log2
  ampOut <- amps(spLog2)
  antimodes <- ampOut$Antimode[,1]
  binInfo <- data.frame("species"=character(), "binName"=character(), "log2BinStart"=numeric(), "log2BinEnd"=numeric(), "numWin"=numeric(), "perData"=numeric(),stringsAsFactors=FALSE)
  k=1
  for (i in 1:(length(antimodes)+1)) {
    if (i==1){
      binStart<-0
    }
    if (i==length(antimodes)+1) {
      binEnd <- max(spLog2)
    } else {
      binEnd <- antimodes[i]
    }
    numWin <- length(which(spLog2>=binStart & spLog2<=binEnd))
    binName <- paste("bin",(k-1),sep="")
    binInfo[k,1] <- spc
    binInfo[k,2] <- binName
    binInfo[k,3] <- binStart
    binInfo[k,4] <- binEnd
    binInfo[k,5] <- numWin
    binInfo[k,6] <- round(100*(numWin/length(spLog2)),2)
    k <- k+1
    binStart <- binEnd
  }
  #head(binInfo)
  numBinsWanted <- max(which(binInfo$numWin>=500))+1
  if (numBinsWanted=="-Inf") {
    numBinsWanted <- 4
  }
  spBinsWanted <- binInfo[1:numBinsWanted,]
  sigSpBins <- rbind(sigSpBins, spBinsWanted[,1:4])
  newColNum <- ncol(sigDataOut)+1
  sigDataOut[newColNum] <- "NA"
  binOutTemp[newColNum] <- "NA"
  for (i in 1:(nrow(spBinsWanted)+1)){
    if (i==(nrow(spBinsWanted)+1)){
      binName <- "aboveUpperBin"
      sigDataOut[which(sigDataOut$log2>=spBinsWanted[(i-1),4] & sigDataOut$species==spc),newColNum] <- binName
    } else {
      binName <- spBinsWanted$binName[i]
      sigDataOut[which(sigDataOut$log2>=spBinsWanted[i,3] & sigDataOut$log2<=spBinsWanted[i,4] & sigDataOut$species==spc),newColNum] <- binName
    }
  }
  sigSpColNames <- c(sigSpColNames, paste(spc, "covBin", sep="-"))
}
colnames(sigDataOut) <- sigSpColNames
sigSpBins$meanCovUpper <- (2^sigSpBins$log2BinEnd)*checkMean
## Can add back once I get the bins working for each spcies
#write.table(sigDataOut, file=paste(outputPrefix, "_contributingSpecies_covBins.txt", sep=""), row.names = F)
#write.table(sigSpBins, file=paste(outputPrefix, "_contributingSpecies_covBinsInfo.txt", sep=""), row.names = F)

binAvgLines <- geom_vline(xintercept = binsWanted$meanCovUpper, color="gray30", linetype="dotdash")
binLogLines <- geom_vline(xintercept =binsWanted$log2BinEnd, color="gray30", linetype="dotdash")
sigBinLogLines <- geom_vline(data=sigSpBins, aes(xintercept = log2BinEnd, color=species), linetype="dotdash")
hBinAvgLines <- geom_hline(yintercept = binsWanted$meanCovUpper, color="gray30", linetype="dotdash")
hBinLogLines <- geom_hline(yintercept =binsWanted$log2BinEnd, color="gray30", linetype="dotdash")
hSigBinLogLines <- geom_hline(data=sigSpBins, aes(yintercept =log2BinEnd, color=species), linetype="dotdash")

pdf(paste(workingDir, outputPrefix, "_covDistPlots.pdf", sep=""))
ggplot(bedData, aes(x=meanValue))+geom_density(aes(y=..scaled..))+ggtitle(paste(outputPrefix, "Windowed Mean Coverage Density Distribution", sep=" ")) + coord_flip() + theme + binAvgLines
ggplot(bedData, aes(x=log2))+geom_density(aes(y=..scaled..))+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Density Distribution", sep=" ")) + coord_flip() + theme + binLogLines

ggplot(bedData, aes(x=log2, fill=species))+geom_density(alpha=.75)+fillScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Species Density Distribution", sep=" ")) + coord_flip()+theme + binLogLines
ggplot(bedData, aes(x=species, y=log2, fill=species))+geom_violin(scale="width")+fillScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Violin Plot", sep=" "))+scale_x_discrete(limits = uniSpecies)+theme + hBinLogLines

ggplot(sigData, aes(x=log2, fill=species))+geom_density(alpha=.75)+fillScale+colorScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Species Density Distribution", sep=" ")) + coord_flip()+theme + binLogLines
ggplot(sigData, aes(x=species, y=log2, fill=species))+geom_violin()+fillScale+colorScale+ggtitle(paste(outputPrefix, "Windowed log2 Coverage Violin Plot", sep=" "))+scale_x_discrete(limits = sigSpecies)+theme  + hBinLogLines

## Can add back once I get the bins working for each spcies
# for (spc in sigSpecies) {
#   spcData <- sigData[which(sigData$species==spc),]
#   spBins <- sigSpBins[which(sigSpBins$species==spc),]
#   binLines <- geom_vline(data=spBins, aes(xintercept = log2BinEnd), color="gray30", linetype="dotdash")
#   hBinLines <- geom_hline(data=spBins, aes(yintercept = log2BinEnd), color="gray30", linetype="dotdash")
#   print(ggplot(spcData, aes(x=log2, fill=species))+geom_density(alpha=1)+fillScale+colorScale+ggtitle(paste(outputPrefix, spc, "Windowed log2 Coverage Species Density Distribution", sep=" ")) + coord_flip()+theme + binLines)
#   print(ggplot(spcData, aes(x=species, y=log2, fill=species))+geom_violin()+fillScale+colorScale+ggtitle(paste(outputPrefix, spc, "Windowed log2 Coverage Violin Plot", sep=" "))+theme  + hBinLines)
# }

dev.off()


colors <- c()
#sigSpecies <- c()
meanData <- data.frame(species=character(), meanValue=numeric(), median=numeric(), log2=numeric())
sigTable <- c()
speciesBreaks <- c()
speciesLabelPos <- c()
chrBreaks <- c()
prevChrEnd <- 0
prevSpcEnd <- 0
#Go though by species and determine mean depth and chromosome lengths.
#Which species will be included on the stacked plots is determined by if the species 99% mean depth is greater than the overall mean depth.
for (i in uniSpecies){
  index <- which(uniSpecies %in% i)
  colors[i] <- colorList[index]
  speciesName <- i
  speciesData <- subset(bedData, species==speciesName)
  speciesStart <- speciesData$Genome_Pos[1]
  speciesBreaks <- append(speciesBreaks, speciesStart)
  speciesEnd <- speciesData$Genome_Pos[length(speciesData$Genome_Pos)]
  speciesMidpoint <- speciesStart+(speciesEnd-speciesStart)/2
  speciesLabelPos <- append(speciesLabelPos, speciesMidpoint)
  if (speciesName %in% sigSpecies) {
    #print(speciesName)
    speciesTable <- speciesData
    speciesTable$speciesPos <- speciesTable$Genome_Pos - prevSpcEnd
    sigTable <- rbind(sigTable, speciesTable)
  }
  spMean <- data.frame(species=i, meanValue=mean(speciesData$meanValue),  median=mean(speciesData$median), log2=mean(speciesData$log2))
  meanData <- rbind(meanData, spMean)
  speciesChrs <- unlist(list(unique(speciesData$chrom)))
  for (chr in speciesChrs) {
    chrData <- subset(speciesData, chrom==chr)
    chrEnd <- chrData$Genome_Pos[length(chrData$Genome_Pos)]
    chrLength <- chrEnd-prevChrEnd
    if (chrLength >= 20*stepSize) {
      chrBreaks <-  append(chrBreaks, chrEnd)  
    }
    prevChrEnd <- chrData$Genome_Pos[length(chrData$Genome_Pos)]
  }
  prevSpcEnd <- speciesEnd
}
#Add extras to plot including vertical lines to delinate species and if labels will be added to axes
vertLines <- geom_vline(xintercept = speciesBreaks)
chrLinesMean <- geom_segment(aes(x = chrBreaks, xend = chrBreaks, y =0, yend=(max(meanQuant)/4)), alpha=0.5, linetype="dotted")
chrLinesLog <- geom_segment(aes(x = chrBreaks, xend = chrBreaks, y =0, yend=(max(logQuant)/4)), alpha=0.5, linetype="dotted")
line <- geom_abline(intercept=0, slope=0)
if (length(uniSpecies)<11){
  xaxis <- scale_x_continuous(breaks=speciesLabelPos, labels=spcLabels, name="Genome Position", limits=c(0,NA))
} else {
  xaxis <- scale_x_continuous(breaks=NULL, labels=NULL, name="Genome Position", limits=c(0,NA))
}
#Assign colors for species 
fillLegend <- scale_fill_manual(name="Species", values = colors, breaks=uniSpecies, labels=gsub("_", " ", uniSpecies))
pointLegend <- scale_color_manual(name="Species", values=colors, breaks=uniSpecies, labels=gsub("_", " ", uniSpecies))

#Plot the data
pdf(paste(workingDir, outputPrefix, "_sppIDerDepthPlot-covBins.pdf", sep=""), width=14)
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))
totalPoint <- geom_point(data=bedData, aes(x=Genome_Pos, y=meanValue, colour=species))
totalPointLim <- geom_point(data=bedData, aes(x=Genome_Pos, y=meanValueLimited, colour=species))
totalFill <- geom_ribbon(data=bedData, aes(x=Genome_Pos, ymin=0, ymax=meanValue, fill=species))
totalFillLim <- geom_ribbon(data=bedData, aes(x=Genome_Pos, ymin=0, ymax=meanValueLimited, fill=species))
yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))
plot(plot+plotTitle+xaxis+yaxis + hBinAvgLines+totalPoint+vertLines+line+ theme_classic()+pointLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic"))) 
plot(plot+plotTitle+xaxis+yaxis + hBinAvgLines+totalFill+vertLines+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
yaxis <- scale_y_continuous(name="log2(Average Depth)", limits = c(0,NA))
plotTitle <- ggtitle(paste(outputPrefix, "Log2(Average Depth)", sep=" "))
totalPointLog <- geom_point(data=bedData, aes(x=Genome_Pos, y=log2, colour=species))
totalFillLog <- geom_ribbon(data=bedData, aes(x=Genome_Pos, ymin=0, ymax=log2, fill=species))
plot(plot+plotTitle+xaxis+yaxis + hBinLogLines+totalPointLog+vertLines+line+ theme_classic()+pointLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
plot(plot+plotTitle+xaxis+yaxis + hBinLogLines+totalFillLog+vertLines+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
yaxis <- scale_y_continuous(name="Average Depth (limited)", limits = c(0,limitValue))
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))
plot(plot+plotTitle+xaxis+yaxis + hBinAvgLines+totalPointLim+vertLines+line+ theme_classic()+pointLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
plot(plot+plotTitle+xaxis+yaxis + hBinAvgLines+totalFillLim+vertLines+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))

#Make stacked plots based on which species have mean 99% is greater than the overall mean
#sigMeanQuant <- quantile(sigTable$meanValue, 0.99, na.rm=T)
#sigSpecies <- unlist(list(unique(sigTable$species)))
#write(paste(outputPrefix, " highCoveSpecies ", sep=""), file=paste(outputPrefix, "_highCovSpecies.txt", sep=""))
#for (spName in sigSpecies) {
#  write(spName, file=paste(outputPrefix, "_highCovSpecies.txt", sep=""), append=T)
#}
endPos <- 0
facetDF <- data.frame()
endsDF <- data.frame()
prevChrEnd <- 0
#Restructure the data for just the species of interest
for (speciesName in sigSpecies) {
  speciesData <- subset(sigTable, species==speciesName)
  sigFinalEnd <- speciesData$speciesPos[length(speciesData$speciesPos)]
  if (sigFinalEnd>endPos)
    endPos <- sigFinalEnd
  speciesChrs <- unlist(list(unique(speciesData$chrom)))
  for (chr in speciesChrs) {
    chrData <- subset(speciesData, chrom==chr)
    chrEnd <- chrData[1,"speciesPos"]
    chrLength <- chrEnd - prevChrEnd
    if (chrLength >= 20*stepSize){
      spcFacetBreaks <- data.frame("Breaks"= chrData[1,"speciesPos"], "species"= gsub("_", " ", speciesName)) 
      facetDF <- rbind(facetDF, spcFacetBreaks) 
    }
    prevChrEnd <- chrEnd
  }
  endsDF <- rbind(endsDF, data.frame("end"=sigFinalEnd, "species"=gsub("_", " ", speciesName)))
}
breakInterval <- endPos/10
breakPos <- seq(0,endPos, by = breakInterval)
breakLabels <- c(toString(0))
#Reformat xaxis so that it shows information about basepair position
for (i in 2:length(breakPos)) {
  level <- breakPos[i]
  if (level > 1000000) {
    label <- paste(toString(round((level/1000000), digits=3)), "Mb", sep="")
    breakLabels <- c(breakLabels, label)
  } else if (1000 < level && level < 1000000){
    label <- paste(toString(round((level/1000000), digits=3)), "Kb", sep="")
    breakLabels <- c(breakLabels, label)
  } else {
    label <- paste(toString(round((level/1000000), digits=3)), "Gb", sep="")
    breakLabels <- c(breakLabels, label)
  }    
}
xaxis <- scale_x_continuous(breaks=breakPos, labels=breakLabels, name="Genome Position", limits=c(0,endPos)) 
plotVert <- theme_classic()
plotVertEnd <- geom_vline(xintercept = endPos)
if(nrow(facetDF) != 0) {
  facetVertLines <- geom_vline(aes(xintercept = Breaks), facetDF, linetype=2)
  plotVertEnd <- geom_vline(aes(xintercept = end), endsDF)
  plotVert <- facetVertLines
}
yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))
sigSpeciesColors <- c()
for (i in sigSpecies){
  index <- which(uniSpecies %in% i)
  sigSpeciesColors[gsub("_", " ", i)] <- colorList[index]
}
sigTable$species <- gsub("_", " ", sigTable$species)
sigSpecies <- gsub("_", " ", sigSpecies)
sigSpBins$species <- gsub("_", " ", sigSpBins$species)

sigSpBins$species <- factor(sigSpBins$species, levels=sigSpecies)
hSigBinAvgLines <- geom_hline(data=sigSpBins, aes(yintercept = meanCovUpper), color="gray30", linetype="dotdash")
hSigBinLogLines <- geom_hline(data=sigSpBins, aes(yintercept =log2BinEnd), color="gray30", linetype="dotdash")

#Plot data
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos, meanValue, colour = species)) + hBinAvgLines+geom_point()+facet_grid(species ~ .)+line+theme_classic()+scale_colour_manual(values=sigSpeciesColors, name="Species", breaks=sigSpecies)+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos)) + hBinAvgLines+geom_ribbon(aes(ymin=0, ymax=meanValue, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd

yaxis <- scale_y_continuous(name="log2(Average Depth)", limits = c(0,NA))
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos, log2, colour = species)) + hBinLogLines+geom_point()+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Log2(Avg depth of coverage)", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos)) + hBinLogLines+geom_ribbon(aes(ymin=0, ymax=log2, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Log2(Avg depth of coverage)", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd

yaxis <- scale_y_continuous(name="Average Depth (limited)", limits = c(0,limitValue))
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos, meanValueLimited, colour = species)) + hBinAvgLines+geom_point()+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos)) + hBinAvgLines+geom_ribbon(aes(ymin=0, ymax=meanValueLimited, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd

dev.off()

