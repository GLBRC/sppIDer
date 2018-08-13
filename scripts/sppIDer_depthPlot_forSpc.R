#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
require(ggplot2)
args <- commandArgs(TRUE)
outputPrefix <- args[1]
#depthOpt <- args[2]

################################################################
# This script makes the depth plots by species
#
#Requirements: ggplot2
#Input: Species averaged depth text file
#
################################################################

# docker vars
#workingDir <- "/tmp/sppIDer/working/"
workingDir <- ""

#Read in data, get info on window size and spread of mean values. Add log2 and a rescaled mean column to be plotted later.
dataFileName <- paste(workingDir, outputPrefix, "_speciesAvgDepth-d.txt", sep="")
if (file.exists(dataFileName)){
  bedData <- read.table(dataFileName, header=F, skip=2, col.names = c("Genome_Pos", "species",  "end", "meanValue", "relativeMean", "max", "median"))
} else {
  dataFileName <- paste(workingDir, outputPrefix, "_speciesAvgDepth-g.txt", sep="")
  bedData <- read.table(dataFileName, header=F, skip=2, col.names = c("Genome_Pos", "species",  "end", "meanValue", "relativeMean", "max", "median"))
}
#bedData <- read.table(paste(workingDir, outputPrefix, "_speciesAvgDepth-", depthOpt, ".txt", sep=""), header=F, skip=2, col.names = c("Genome_Pos", "species",  "end", "meanValue", "relativeMean", "max", "median"))
completeChromList <- unlist(list(unique(bedData$species)))
## INLINE TESTS
#print ("============== length(completeChromList)")
#print (length(completeChromList))
output <- data.frame()

species <- c()
#Split chr name from spcecies name
for (i in 1:length(completeChromList)) {
  completeChrString <- toString(completeChromList[[i]])
  ## INLINE TESTS
  #print ('============== completeChrString')
  #print (completeChrString)
  split <- unlist(strsplit(completeChrString, "-"))
  species <- c(species, split[1])
}
uniSpecies <- unique(species)

plot <- ggplot()
checkMean <- mean(bedData$meanValue)
maxMean <- max(bedData$meanValue)
meanQuant <- quantile(bedData$meanValue, 0.99, na.rm=T)
bedData$log2 <- log2(bedData$meanValue/checkMean)
bedData$log2[bedData$log2<0] <- 0
logQuant <- quantile(bedData$log2, 0.99, na.rm=T)

#Assign colors based on number of species. Here you can use your own color scheme. 
pdf(NULL)
if (length(uniSpecies)>1) {
  colorList <- palette(rainbow(length(uniSpecies)))
  colorList <- palette(rainbow(length(uniSpecies)))
} else {
  colorList <- c("red")
}
spcLabels <- gsub("_", "\n", uniSpecies)

#Assign colors for each species
colors <- c()
for (k in 1:length(uniSpecies)){
  speciesName <- uniSpecies[[k]]
  color <- colorList[k]
  colors[speciesName] <- color
}
speciesBreaks <- c()
speciesLabel <- c()
labelPos <- c()
sigSpecies <- c()
legendValues <- c()
sigTable <- c()
#Cycle through Species
#pull out important species by if the spcies mean cov is greater than the overall mean
spcLabeled <- data.frame()
for (k in 1:length(uniSpecies)) {
  speciesName <- toString(uniSpecies[[k]])
  spc <- subset(bedData, grepl(speciesName, species))
  subLen <- spc$end[length(spc$end)] 
  labelPos <- append(labelPos, (spc[1,"Genome_Pos"])+subLen/2)
  speciesBreaks <- append(speciesBreaks, spc[1,"Genome_Pos"])
  spc <- cbind(species_name = speciesName, spc)
  spcLabeled <- rbind(spcLabeled, spc)
  spc99Mean <- quantile(spc$meanValue, 0.99)[[1]]
  if (is.na(spc99Mean)){
    spc99Mean <- 0
  }
  if (spc99Mean>checkMean) {
    sigSpecies <- c(sigSpecies, speciesName)
    #speciesTable <- cbind(spc_Chr_Pos= seq(0, (length(spc$chrom)-1)*stepSize, stepSize), spc)
    sigTable <- rbind(sigTable, spc)
    meanQuant <- quantile(spc$meanValue, 0.99, na.rm=T)
    logQuant <- quantile(spc$log2, 0.99, na.rm=T)
  }
  speciesLabel <- append(speciesLabel, speciesName)
  legendValues = c(legendValues, paste(speciesName, colorList[k]), sep="=")
}
#Add extras to plot including vertical lines to delinate species and if labels will be added to axes
GenomeEnds <- spcLabeled$Genome_Pos-1
GenomeEnds <- GenomeEnds[-1]
finalEnd <- spcLabeled$Genome_Pos[length(spcLabeled$Genome_Pos)]+spcLabeled$end[length(spcLabeled$end)]
GenomeEnds <- c(GenomeEnds, finalEnd)
spcLabeledDup <- spcLabeled
spcLabeledDup$Genome_Pos <- GenomeEnds
spcLabeledDup <- rbind(spcLabeled, spcLabeledDup)
  
totalFill <- geom_ribbon(data=spcLabeledDup, aes(x=Genome_Pos, ymin=0, ymax=meanValue, fill=species_name))
totalFillMean <- geom_ribbon(data=spcLabeledDup, aes(x=Genome_Pos, ymin=0, ymax=log2, fill=species_name))
totalPoint <- geom_point(data=spcLabeled, aes(x=Genome_Pos, y=meanValue, colour=species_name))
totalPointMean <- geom_point(data=spcLabeled, aes(x=Genome_Pos, y=log2, colour=species_name))
line <- geom_abline(intercept=0, slope=0)
vertLines <- geom_vline(xintercept = speciesBreaks)
#Assign colors for species 
fillLegend <- scale_fill_manual(name="Species", values = colors, breaks=uniSpecies, labels=gsub("_", " ", uniSpecies))
pointLegend <- scale_color_manual(name="Species", values=colors, breaks=uniSpecies, labels=gsub("_", " ", uniSpecies))

pdf(paste(workingDir, outputPrefix, "_speciesDepth.pdf", sep=""), width=14)

if (length(uniSpecies)<11){
  xaxis <- scale_x_continuous(breaks=labelPos, labels=spcLabels, name="Genome Position", limits=c(0,NA))
} else {
  xaxis <- scale_x_continuous(breaks=NULL, labels=NULL, name="Genome Position", limits=c(0,NA))
}
#Plot the data
yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,maxMean+(.1*maxMean)))
plotTitle <- ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+totalFill+vertLines+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))
yaxis <- scale_y_continuous(name="log2(avg/whole genome avg)", limits = c(0,max(spcLabeled$log2)+(.1*max(spcLabeled$log2))))
plotTitle <- ggtitle(paste(outputPrefix, "log2 Mean Avg depth of coverage", sep=" "))
plot(plot+plotTitle+xaxis+yaxis+totalFillMean+vertLines+line+ theme_classic()+fillLegend+theme(axis.text.x=element_text(face="italic"), legend.text=element_text(face="italic")))

dev.off()
