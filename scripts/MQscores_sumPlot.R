#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
require(ggplot2)
args <- commandArgs(TRUE)
strainName <- args[1]

################################################################
# This script makes plots of percent mapped reads and MQ
#
#Requirements: ggplot2
#Input: text file of the MQ scores parsed from the sam file.
#
################################################################

# docker vars
workingDir <- "/tmp/sppIDer/working/"
#workingDir <- ""

#Read in data and determine species
chiSqOutName <- paste(workingDir, strainName, "_MQ_chiSquared.txt", sep="")
summaryOutputName <- paste(workingDir, strainName, "_MQsummary.txt", sep="")
plotOutputName <- paste(workingDir, strainName, "_plotMQ.pdf", sep="")
MQdf <- read.table(paste(workingDir, strainName, "_MQ.txt", sep=""), header=T)
levels(MQdf$Species) <- c(levels(MQdf$Species), "Unmapped")
MQdf$Species[MQdf$Species=="*"] <- "Unmapped"
MQdfForChiSq <- MQdf
uniSpecies <- unlist(list(unique(MQdf$Species)))
uniSpecies <- factor(uniSpecies, levels=uniSpecies)
spRe <- " "
if (length(uniSpecies)<=11) {
  spRe <- "\n"
}

#Count reads and subset the data based on Mapping Quality
numReads <- sum(MQdf$count)
maps <- subset(MQdf, Species!="Unmapped")
mapsNonZero <- subset(maps, MQscore>0)
MQscoresZero <- subset(MQdf, MQscore==0)
MQscoresNonZero <- subset(MQdf, MQscore>0)
propZero <- (sum(MQscoresZero$count)/sum(MQdf$count))*100
#Print summary stats about mapped and unmapped reads
write(paste(strainName, " Num reads = ", toString(numReads), sep=""), file=summaryOutputName)
write(paste(strainName, " Num mapped reads = ", toString(sum(maps$count)), sep=""), file=summaryOutputName, append=T)
write(paste(strainName, " Unmapped reads = ", toString(round(propZero, digits=2)), "%", sep=""), file=summaryOutputName, append=T)
write(paste(strainName, " Average MQ = ", toString(weighted.mean(MQdf$MQscore, MQdf$count)), sep=""), file=summaryOutputName, append=T)
write(paste(strainName, " Median MQ = ", toString(median(rep(MQdf$MQscore, times=MQdf$count))), sep=""), file=summaryOutputName, append=T)
write("Species\tTotal mapped reads\t% of all reads\t% Nonzero mapped reads\tAll average MQ\tNonZero average MQ\tAll Median MQ\tNonzero Median MQ", file=summaryOutputName, append=T)

MQprop <- data.frame("Species"=character(), "Percentage"=numeric())
MQpropPos <- data.frame()
MQcounts <- data.frame(Species=character(), totalCount = numeric(), nonZeroCount = numeric())
speciesCountDF <- data.frame()
violin_df <- data.frame()
#Go though data by species and quantify how many and how well reads mapped
for (species in uniSpecies){
    speciesNonZeroMQscores <- subset(MQscoresNonZero, Species==species)
    speciesCount <- sum(speciesNonZeroMQscores$count)
    speciesCountDF <- rbind(speciesCountDF, data.frame("species"=species, "sum"=speciesCount))
    propSpeciesNonZero <-  (speciesCount/sum(MQscoresNonZero$count))*100
    speciesMQscores <- subset(MQdf, Species==species)
    speciesTotalCount <- sum(speciesMQscores$count)
    propSpecies <- (sum(speciesMQscores$count)/sum(MQdf$count))*100
    speciesMean <- weighted.mean(speciesMQscores$MQscore, speciesMQscores$count)
    if (is.na(speciesMean)){
      speciesMean <- 0
    }
    speciesMedian <- median(rep(speciesMQscores$MQscore, times=speciesMQscores$count))
    if (is.na(speciesMedian)){
      speciesMedian <- 0
    }
    speciesMeanNonZero <- weighted.mean(speciesNonZeroMQscores$MQscore, speciesNonZeroMQscores$count)
    if (is.na(speciesMeanNonZero)){
      speciesMeanNonZero <- 0
    }
    speciesMedianNonZero <- median(rep(speciesNonZeroMQscores$MQscore, times=speciesNonZeroMQscores$count))
    if (is.na(speciesMedianNonZero)){
      speciesMedianNonZero <- 0
    }
    write(paste(species, "\t", toString(speciesTotalCount), "\t", toString(round(propSpecies, digits=2)), "%\t", toString(round(propSpeciesNonZero, digits=2)), "%\t", toString(round(speciesMean, digits=2)), "\t", toString(round(speciesMeanNonZero, digits=2)), "\t", toString(round(speciesMedian, digits=2)), "\t", toString(round(speciesMedianNonZero, digits=2)),  sep=""), file=summaryOutputName, append=T)
    MQprop <- rbind(MQprop, data.frame("Species"=species, "Percentage"=propSpecies))
    MQpropPos <- rbind(MQpropPos, data.frame("Species"=species, "Percentage"=propSpeciesNonZero))
    MQcounts <- rbind(MQcounts, data.frame(Species=species, totalCount=speciesTotalCount, nonZeroCount =speciesCount))
    if (speciesCount>0){
      violin_df <- rbind(violin_df, MQdf[which(MQdf$Species==species),])
    }
}
MQpropPos <- MQpropPos[-1,]
noDataSpecies <- MQcounts$Species[which(MQcounts$totalCount==0)]
mapsExisting <- maps
for(sp in noDataSpecies){
  mapsExisting <- mapsExisting[-which(mapsExisting$Species==sp),]
}

###ChiSquared test###
countSpecies <- matrix(0,nrow=61, ncol=1)
countSpecies <- data.frame(countSpecies)
for (species in uniSpecies){
  spData <- MQdfForChiSq[which(MQdfForChiSq$Species==species),]
  countData <- spData$count
  if (species=="Unmapped"){
    extraData <- rep(0,60)
    countData <- c(countData, extraData)
  }
  countSpecies <- cbind(countSpecies, data.frame(species=countData))
}
header <- c("temp", as.character(uniSpecies))
colnames(countSpecies) <- header
rownames(countSpecies) <- c(0:60)
countSpecies <- countSpecies[,-1]

#Removing Zeros and Unmapped
unmappedCount <- countSpecies[1,1][[1]]
posCount <- countSpecies[-1,-1]
posSum <- colSums(posCount)
posSum <- c(posSum, unmappedCount)
names(posSum)[length(posSum)] <- "Unmapped"
countChisq <- chisq.test(posSum)
#print(countChisq$p.value)
#print(countChisq$residuals)
countRes <- t(data.frame(countChisq$residuals))
namesRes <- t(data.frame(names(countChisq$residuals)))
countRes <- rbind(namesRes, countRes)

MQ60 <- posCount[60,]
mq60chisq <- chisq.test(x=MQ60)
mq60Res <- t(data.frame(mq60chisq$residuals))
namesRes60 <- namesRes[-length(namesRes)]
mq60Res <- rbind(namesRes60, mq60Res)

write(paste(strainName, " count of reads mapped - Chi Squared P-value ", toString(countChisq$p.value), sep=""), file=chiSqOutName)
write(countRes, file=chiSqOutName, append=T)
write(paste(strainName, " MQ60 - Chi Squared P-value ", toString(mq60chisq$p.value), sep=""), file=chiSqOutName, append=T)
write(mq60Res, file=chiSqOutName, append=T)

#Set color for number of species
pdf(NULL)
if (length(uniSpecies)>2) {
  colorList <- palette(rainbow(length(uniSpecies)-1))
  colorList <- palette(rainbow(length(uniSpecies)-1))
} else {
  colorList <- c("red")
}
#Assign colors for each species
colors <- c()
for (i in uniSpecies[2:length(uniSpecies)]){
  index <- which(uniSpecies %in% i)-1
  colors[gsub("_", spRe, i)] <- colorList[index]
}
colors["Unmapped"] <- "darkgrey"
#Replace "_" with spaces for dataframes
MQdf$Species <- gsub("_", spRe, MQdf$Species)
MQpropPos$Species <- gsub("_", spRe, MQpropPos$Species)
MQprop$Species <- gsub("_", spRe, MQprop$Species)
maps$Species <- gsub("_", spRe, maps$Species)
mapsExisting$Species <- gsub("_", spRe, mapsExisting$Species)
mapsNonZero$Species <- gsub("_", spRe, mapsNonZero$Species)
spcLabels <- gsub("_", spRe, uniSpecies)
violin_df$Species <- gsub("_", spRe, violin_df$Species)
#Set color legends
fillLegend <- scale_fill_manual(name="Species", values = colors, breaks=spcLabels)
colorLegend <- scale_colour_manual(name="Species", values = colors, breaks=spcLabels)
theme <- theme_bw()+theme(legend.text=element_text(face="italic"), axis.text.x=element_text(face="italic"))
#Remove labels if too many species will make it cluttered.
if (length(uniSpecies)>11){
  themeUnmap <- theme_bw()+theme(axis.text.x=element_blank(), legend.text=element_text(face="italic"))
  theme <- theme_bw()+theme(axis.text.x=element_blank(), legend.text=element_text(face="italic"))
} else if (length(uniSpecies)==11) { 
  themeUnmap <- theme_bw()+theme(axis.text.x=element_blank(), legend.text=element_text(face="italic"))
} else {
  themeUnmap <- theme 
}

#Plot the data
pdf(plotOutputName, compress=T, width=10)
ggplot(MQprop, aes(x=factor(Species, levels=spcLabels)))+geom_bar(aes(fill=Species, weight=Percentage))+ggtitle(paste(strainName, "Mapping bar plot w/ unmapped", sep=" "))+labs(x = "Species", y="Percentage")+fillLegend+themeUnmap+scale_y_continuous(limits=c(0,100))
ggplot(MQpropPos, aes(x=factor(Species, levels=spcLabels)))+geom_bar(aes(fill=Species, weight=Percentage))+ggtitle(paste(strainName, "Mapping bar plot w/out unmapped", sep=" "))+labs(x = "Species", y = "Percentage")+fillLegend+theme(axis.text.x=element_text(face="italic"))+theme+theme(legend.text=element_text(face="italic"))+scale_y_continuous(limits=c(0,100))
options(warn = -1)
ggplot(mapsExisting, aes(factor(Species, levels=spcLabels), MQscore, weight=count))+ geom_violin(bw=1, scale = "count", draw_quantiles=c(.25,.5,.75), aes(fill=Species))+fillLegend+colorLegend+ggtitle(paste(strainName, "Mapping quality of species with mapped reads", sep=" "))+labs(x = "Species")+theme(axis.text.x=element_text(face="italic"))+theme+theme(legend.text=element_text(face="italic"))
#ggplot(mapsNonZero, aes(factor(Species, levels=spcLabels), MQscore, weight=count))+ geom_violin(bw=1, scale = "count", draw_quantiles=c(.25,.5,.75), aes(fill=Species))+fillLegend+colorLegend+ggtitle(paste(strainName, "Positive Mapping violin (count scaled)", sep=" "))+ labs(x = "Species")+theme(axis.text.x=element_text(face="italic"))+theme+theme(legend.text=element_text(face="italic"))

options(warn = 1)

dev.off()
