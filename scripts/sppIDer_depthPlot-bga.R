#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
require(ggplot2)
plot <- ggplot()
args <- commandArgs(TRUE)
outputPrefix <- args[1]

################################################################
# This script makes the windowed depth plots with the -byGroup option
#
#Requirements: ggplot2
#Input: Window averaged depth text file
#
################################################################

#Read in data, get info on window size and spread of mean values. Add log2 and a rescaled mean column to be plotted later.
bedData <- read.table(paste(outputPrefix, "winAvgDepth-g.txt", sep="_"), header=T)
bedData <- read.table(paste(outputPrefix, "winAvgDepth-g.txt", sep="_"), header=F, skip=2, col.names = c("Genome_Pos", "species", "chrom", "start",  "end", "meanValue", "log2mean", "max", "median"))
completeChromList <- unlist(list(unique(bedData$chrom)))
stepSize <- bedData[2,1]
checkMean <- mean(bedData$meanValue)
meanQuant <- quantile(bedData$meanValue, 0.99, na.rm=T)
bedData$meanValueLimited <- bedData$meanValue
bedData$meanValueLimited[bedData$meanValueLimited>meanQuant] <- meanQuant
bedData$log2 <- log2(bedData$meanValue/checkMean)
bedData$log2[bedData$log2<0] <- 0
logQuant <- quantile(bedData$log2, 0.99, na.rm=T)

#Assign colors based on number of species. Here you can use your own color scheme. 
pdf(NULL)
uniSpecies <- unlist(list(unique(bedData$species)))
if (length(uniSpecies)>1) {
  colorList <- palette(rainbow(length(uniSpecies)))
  colorList <- palette(rainbow(length(uniSpecies)))
} else {
  colorList <- c("red")
}
spcLabels <- gsub("_", "\n", uniSpecies)

colors <- c()
sigSpecies <- c()
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
  spc99Mean <- quantile(speciesData$meanValue, 0.985)[[1]]
  if (is.na(spc99Mean)){
    spc99Mean <- 0
  }
  if (spc99Mean>checkMean) {
    speciesTable <- speciesData
    spcStart <- speciesTable$Genome_Pos[1]
    speciesTable$speciesPos <- speciesTable$Genome_Pos - prevSpcEnd
    sigTable <- rbind(sigTable, speciesTable)
  }
  speciesChrs <- unlist(list(unique(speciesData$chrom)))
  for (chr in speciesChrs) {
    chrData <- subset(speciesData, chrom==chr)
    chrEnd <- chrData$Genome_Pos[length(chrData$Genome_Pos)]
    chrLength <- chrEnd-prevChrEnd
    if (chrLength >= 10*stepSize) {
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
pdf(paste(outputPrefix, "sppIDerDepthPlot-g.pdf", sep='_'), width=14)
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))
totalPoint <- geom_point(data=bedData, aes(x=Genome_Pos, y=meanValue, colour=species))
totalFill <- geom_ribbon(data=bedData, aes(x=Genome_Pos, ymin=0, ymax=meanValue, fill=species))
totalPointLim <- geom_point(data=bedData, aes(x=Genome_Pos, y=meanValueLimited, colour=species))
totalFillLim <- geom_ribbon(data=bedData, aes(x=Genome_Pos, ymin=0, ymax=meanValueLimited, fill=species))
yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))
plot(plot+plotTitle+xaxis+yaxis+totalPoint+vertLines+chrLinesMean+line+ theme_classic()+pointLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
plot(plot+plotTitle+xaxis+yaxis+totalFill+vertLines+chrLinesMean+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
yaxis <- scale_y_continuous(name="log2(Average Depth)", limits = c(0,NA))
plotTitle <- ggtitle(paste(outputPrefix, "Log2(Average Depth)", sep=" "))
totalPointLog <- geom_point(data=bedData, aes(x=Genome_Pos, y=log2, colour=species))
totalFillLog <- geom_ribbon(data=bedData, aes(x=Genome_Pos, ymin=0, ymax=log2, fill=species))
plot(plot+plotTitle+xaxis+yaxis+totalPointLog+vertLines+chrLinesLog+line+ theme_classic()+pointLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
plot(plot+plotTitle+xaxis+yaxis+totalFillLog+vertLines+chrLinesLog+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))
yaxis <- scale_y_continuous(name="Average Depth (limited)", limits = c(0,meanQuant))
plot(plot+plotTitle+xaxis+yaxis+totalPointLim+vertLines+chrLinesMean+line+ theme_classic()+pointLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
plot(plot+plotTitle+xaxis+yaxis+totalFillLim+vertLines+chrLinesMean+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))

#Make stacked plots based on which species have mean 99% is greater than the overall mean
sigMeanQuant <- quantile(sigTable$meanValue, 0.99, na.rm=T)
sigSpecies <- unlist(list(unique(sigTable$species)))
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
    if (chrLength >= 10*stepSize){
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

#Plot data
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos, meanValue, colour = species))+geom_point()+facet_grid(species ~ .)+line+theme_classic()+scale_colour_manual(values=sigSpeciesColors, name="Species", breaks=sigSpecies)+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos))+geom_ribbon(aes(ymin=0, ymax=meanValue, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd

yaxis <- scale_y_continuous(name="log2(Average Depth)", limits = c(0,NA))
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos, log2, colour = species))+geom_point()+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Log2(Avg depth of coverage)", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos))+geom_ribbon(aes(ymin=0, ymax=log2, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Log2(Avg depth of coverage)", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd

yaxis <- scale_y_continuous(name="Average Depth (limited)", limits = c(0,meanQuant))
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos, meanValueLimited, colour = species))+geom_point()+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd
ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(speciesPos))+geom_ribbon(aes(ymin=0, ymax=meanValueLimited, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=sigSpeciesColors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+theme(legend.text=element_text(face="italic"), strip.text=element_text(face="italic"))+plotVert+geom_vline(xintercept = 0)+plotVertEnd

dev.off()