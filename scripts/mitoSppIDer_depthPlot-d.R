require(ggplot2)
plot <- ggplot()
args <- commandArgs(TRUE)
outputPrefix <- args[1]

################################################################
# This script makes the depth plots for mito sppIDer
#
#Requirements: ggplot2
#Input: Window averaged depth text file, and gff of just CDS if available. 
#
################################################################


#Read in data, get info on window size and spread of mean values. Add log2 and a rescaled mean column to be plotted later.
bedData <- read.table(paste(outputPrefix, "winAvgDepth-d.txt", sep="_"), header=T)
bedData <- read.table(paste(outputPrefix, "winAvgDepth-d.txt", sep="_"), header=F, skip=2, col.names = c("Genome_Pos", "species", "chrom", "start",  "end", "meanValue", "relativeMean", "max", "median"))
#Read in gff dat if it exists
if (length(args)==2){
  gffKey <- read.table(args[2], header=T)
  names(gffKey)[names(gffKey)=="Species"] <- "species"
}
covMeanFile <- paste(outputPrefix, "mitoAvgDepthSummary-d.txt", sep="_")
stepSize <- bedData[2,1]
checkMean <- mean(bedData$meanValue)
meanQuant <- quantile(bedData$meanValue, 0.97, na.rm=T)
maxMean <- max(bedData$meanValue)
mean99 <- max(quantile(bedData$meanValue, 0.99, na.rm=T))
mean992 <- max(quantile(bedData$meanValue, 0.992, na.rm=T))
bedData$meanValueLimited <- bedData$meanValue
bedData$meanValueLimited[bedData$meanValueLimited>meanQuant] <- mean99
bedData$log2 <- log2(bedData$meanValue/checkMean)
bedData$log2[bedData$log2<0] <- 0
logQuant <- quantile(bedData$log2, 0.99, na.rm=T)
maxLog2 <- max(bedData$log2)

#Assign colors based on number of species. Here you can use your own color scheme. 
pdf(NULL)
uniSpecies <- unlist(list(unique(bedData$species)))
if (length(uniSpecies)>1) {
  colorList <- palette(rainbow(length(uniSpecies)))
  colorList <- palette(rainbow(length(uniSpecies)))
} else {
  colorList <- c("red")
}
dev.off()

colors <- c()
sigSpecies <- c()
sigTable <- data.frame()
#Go though by species and determine mean depth and chromosome lengths.
#Which species will be included on the stacked plots is determined by if the species 99% mean depth is greater than the overall mean depth.
for (i in uniSpecies){
  index <- which(uniSpecies %in% i)
  colors[i] <- colorList[index]
  speciesData <- subset(bedData, species==i)
  spc99Mean <- quantile(speciesData$meanValue, 0.99)[[1]]
  if (is.na(spc99Mean)){
    spc99Mean <- 0
  }
  if (spc99Mean>meanQuant) {
    sigSpecies <- append(sigSpecies, i)
    sigTable <- rbind(sigTable, speciesData)
  }
}
pdf(paste(outputPrefix, "mitoSppIDerDepthPlot-d.pdf", sep='_'), width=14, onefile = TRUE)
line <- geom_abline(intercept=0, slope=0)
xaxis <- scale_x_continuous(name="Genome Position", limits=c(0,NA)) 
yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))

#If there is a gff key then assign where the coding regions are on the plot and plot the data
if (exists("gffKey")){
  print("wGFF")
  yLabPosFull <- c()
  yLabPosLimited <- c()
  yLabPosLog2 <- c()
  for (i in 1:((length(gffKey$End))/3)){
    yLabPosFull <- append(yLabPosFull, 0)
    yLabPosFull <- append(yLabPosFull, maxMean/2)
    yLabPosFull <- append(yLabPosFull, maxMean)
    yLabPosLimited <- append(yLabPosLimited, 0)
    yLabPosLimited <- append(yLabPosLimited, mean99/2)
    yLabPosLimited <- append(yLabPosLimited, mean99)
    yLabPosLog2 <- append(yLabPosLog2, 0)
    yLabPosLog2 <- append(yLabPosLog2, logQuant/2)
    yLabPosLog2 <- append(yLabPosLog2, logQuant)
  }
  yLabPosFull <- yLabPosFull[1:length(gffKey$End)]
  yLabPosLimited <- yLabPosLimited[1:length(gffKey$End)]
  yLabPosLog2 <- yLabPosLog2[1:length(gffKey$End)]
  gffKey$yPosFull <- yLabPosFull
  gffKey$yPosLimited <- yLabPosLimited
  gffKey$yPosLog2 <- yLabPosLog2
  bottomLabs <- subset(gffKey, yPosLimited==0)
  middleLabs <- subset(gffKey, yPosLimited==mean99/2)
  topLabs <- subset(gffKey, yPosLimited==mean99)
  bottomLabsLog2 <- subset(gffKey, yPosLog2==0)
  middleLabsLog2 <- subset(gffKey, yPosLog2==logQuant/2)
  topLabsLog2 <- subset(gffKey, yPosLog2==logQuant)
  ##Full all species
  CDSstartLines <- geom_segment(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(x = Start, xend = Start, y =0, yend=maxMean), colour="black", alpha=0.5)
  CDSendLines <- geom_segment(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(x = End, xend = End, y =0, yend=maxMean), colour="black", alpha=0.5)
  shade <- geom_rect(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(xmin = Start, xmax = End, ymin =0, ymax=maxMean), alpha=0.5)
  bottomLabels <- geom_text(data = transform(bottomLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'bottom', nudge_y=0.5, aes(x=Midpoint, label =Name, y=yPosFull))
  middleLabels <- geom_text(data = transform(middleLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'center', aes(x=Midpoint, label =Name, y=yPosFull))
  topLabels <- geom_text(data = transform(topLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 1, aes(x=Midpoint, label =Name, y=yPosFull))
  xaxis <- scale_x_continuous(name="Genome Position", limits=c(0,NA)) 
  yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))
  print(ggplot()+shade+geom_point(data =transform(bedData, species=factor(species, levels=uniSpecies)), aes(start, meanValue, colour = species), size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species", breaks=uniSpecies)+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_ribbon(data =transform(bedData, species=factor(species, levels=uniSpecies)), aes(x=start, ymin=0, ymax=meanValue, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  #Log2 all species
  CDSstartLines <- geom_segment(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(x = Start, xend = Start, y =0, yend=maxLog2), colour="black", alpha=0.5)
  CDSendLines <- geom_segment(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(x = End, xend = End, y =0, yend=maxLog2), colour="black", alpha=0.5)
  shade <- geom_rect(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(xmin = Start, xmax = End, ymin =0, ymax=maxLog2), alpha=0.5)
  bottomLabels <- geom_text(data = transform(bottomLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'bottom', nudge_y=0.5, aes(x=Midpoint, label =Name, y=yPosLog2))
  middleLabels <- geom_text(data = transform(middleLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'center', aes(x=Midpoint, label =Name, y=yPosLog2))
  topLabels <- geom_text(data = transform(topLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 1, aes(x=Midpoint, label =Name, y=yPosLog2))
  yaxis <- scale_y_continuous(name="log2(Average Depth)", limits = c(0,NA))
  print(ggplot()+shade+geom_point(data=transform(bedData, species=factor(species, levels=uniSpecies)), aes(start, log2, colour = species),size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Log2 Avg depth of coverage", sep=" "))+CDSstartLines+CDSendLines+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_point(data=transform(bedData, species=factor(species, levels=uniSpecies)), aes(start, log2, colour = species),size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Log2 Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_ribbon(data=transform(bedData, species=factor(species, levels=uniSpecies)), aes(x=start, ymin=0, ymax=log2, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  #Limited all species
  CDSstartLines <- geom_segment(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(x = Start, xend = Start, y =0, yend=mean99), colour="black", alpha=0.5)
  CDSendLines <- geom_segment(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(x = End, xend = End, y =0, yend=mean99), colour="black", alpha=0.5)
  shade <- geom_rect(data=transform(gffKey, species=factor(species, levels=uniSpecies)), aes(xmin = Start, xmax = End, ymin =0, ymax=mean99), alpha=0.5)
  bottomLabels <- geom_text(data = transform(bottomLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'bottom', nudge_y=0.5, aes(x=Midpoint, label =Name, y=yPosLimited))
  middleLabels <- geom_text(data = transform(middleLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'center', aes(x=Midpoint, label =Name, y=yPosLimited))
  topLabels <- geom_text(data = transform(topLabs, species=factor(species, levels=uniSpecies)), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 1, aes(x=Midpoint, label =Name, y=yPosLimited))
  yaxis <- scale_y_continuous(name="Average Depth (limited)", limits = c(0,mean99))
  print(ggplot()+shade+geom_point(data=transform(bedData, species=factor(species, levels=uniSpecies)), aes(start, meanValueLimited, colour = species),size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_ribbon(data=transform(bedData, species=factor(species, levels=uniSpecies)),aes(x=start, ymin=0, ymax=meanValueLimited, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  ##For just sig species
  sigGFF <- data.frame()
  for (speciesName in sigSpecies){
    spcGFF <- subset(gffKey, species==speciesName)
    sigGFF <- rbind(sigGFF, spcGFF)
  }
  bottomLabs <- subset(sigGFF, yPosLimited==0)
  middleLabs <- subset(sigGFF, yPosLimited==mean99/2)
  topLabs <- subset(sigGFF, yPosLimited==mean99)
  CDSstartLines <- geom_segment(data=transform(sigGFF, species=factor(species, levels=sigSpecies), rm.na=T), aes(x = Start, xend = Start, y =0, yend=maxMean), colour="black", alpha=0.5)
  CDSendLines <- geom_segment(data=transform(sigGFF, species=factor(species, levels=sigSpecies), rm.na=T), aes(x = End, xend = End, y =0, yend=maxMean), colour="black", alpha=0.5)
  shade <- geom_rect(data=transform(sigGFF, species=factor(species, levels=uniSpecies)), aes(xmin = Start, xmax = End, ymin =0, ymax=maxMean), alpha=0.5)
  bottomLabels <- geom_text(data = transform(bottomLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'bottom', nudge_y=0.5, aes(x=Midpoint, label =Name, y=yPosFull))
  middleLabels <- geom_text(data = transform(middleLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'center', aes(x=Midpoint, label =Name, y=yPosFull))
  topLabels <- geom_text(data = transform(topLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 1, aes(x=Midpoint, label =Name, y=yPosFull))
  yaxis <- scale_y_continuous(name="Average Depth", limits = c(0,NA))
  print(ggplot()+shade+geom_point(data=transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start, meanValue, colour = species), size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_ribbon(data=transform(sigTable, species=factor(species, levels=sigSpecies)), aes(x=start, ymin=0, ymax=meanValue, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  #Log2 for sig species
  CDSstartLines <- geom_segment(data=transform(sigGFF, species=factor(species, levels=sigSpecies), rm.na=T), aes(x = Start, xend = Start, y =0, yend=maxLog2), colour="black", alpha=0.5)
  CDSendLines <- geom_segment(data=transform(sigGFF, species=factor(species, levels=sigSpecies), rm.na=T), aes(x = End, xend = End, y =0, yend=maxLog2), colour="black", alpha=0.5)
  shade <- geom_rect(data=transform(sigGFF, species=factor(species, levels=uniSpecies)), aes(xmin = Start, xmax = End, ymin =0, ymax=maxLog2), alpha=0.5)
  bottomLabels <- geom_text(data = transform(bottomLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'bottom', nudge_y=0.5, aes(x=Midpoint, label =Name, y=yPosLog2))
  middleLabels <- geom_text(data = transform(middleLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'center', aes(x=Midpoint, label =Name, y=yPosLog2))
  topLabels <- geom_text(data = transform(topLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 1, aes(x=Midpoint, label =Name, y=yPosLog2))
  yaxis <- scale_y_continuous(name="log2(Average Depth)", limits = c(0,NA))
  print(ggplot()+shade+geom_point(data=transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start, log2, colour = species), size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_ribbon(data=transform(sigTable, species=factor(species, levels=sigSpecies)), aes(x=start, ymin=0, ymax=log2, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  #Limited for sig species
  CDSstartLines <- geom_segment(data=transform(sigGFF, species=factor(species, levels=sigSpecies), rm.na=T), aes(x = Start, xend = Start, y =0, yend=mean99), colour="black", alpha=0.5)
  CDSendLines <- geom_segment(data=transform(sigGFF, species=factor(species, levels=sigSpecies), rm.na=T), aes(x = End, xend = End, y =0, yend=mean99), colour="black", alpha=0.5)
  shade <- geom_rect(data=transform(sigGFF, species=factor(species, levels=uniSpecies)), aes(xmin = Start, xmax = End, ymin =0, ymax=mean99), alpha=0.5)
  bottomLabels <- geom_text(data = transform(bottomLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'bottom', nudge_y=0.5, aes(x=Midpoint, label =Name, y=yPosLimited))
  middleLabels <- geom_text(data = transform(middleLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 'center', aes(x=Midpoint, label =Name, y=yPosLimited))
  topLabels <- geom_text(data = transform(topLabs, species=factor(species, levels=sigSpecies), rm.na=T), colour="black", size = 2.5, show.legend=F, check_overlap = T, angle = 0, vjust = 1, aes(x=Midpoint, label =Name, y=yPosLimited))
  yaxis <- scale_y_continuous(name="Average Depth (limited)", limits = c(0,mean99))
  print(ggplot()+shade+geom_point(data=transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start, meanValueLimited, colour = species), size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
  print(ggplot()+shade+geom_ribbon(data=transform(sigTable, species=factor(species, levels=sigSpecies)), aes(x=start, ymin=0, ymax=meanValueLimited, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" "))+bottomLabels+middleLabels+topLabels)
} else {
  ###Without GFF
  plot(ggplot(transform(bedData, species=factor(species, levels=uniSpecies)))+geom_point(aes(start, meanValue, colour = species), size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species", breaks=uniSpecies)+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" ")))
  plot(ggplot(transform(bedData, species=factor(species, levels=uniSpecies)), aes(start))+geom_ribbon(aes(ymin=0, ymax=meanValue, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" ")))
  plot(ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start, meanValue, colour = species))+geom_point(size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" ")))
  plot(ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start))+geom_ribbon(aes(ymin=0, ymax=meanValue, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" ")))
  yaxis <- scale_y_continuous(name="log2 Average Depth", limits = c(0,NA))
  plot(ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start, log2, colour = species))+geom_point(size=0.5)+facet_grid(species ~ .)+theme_classic()+line+scale_colour_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" ")))
  plot(ggplot(transform(sigTable, species=factor(species, levels=sigSpecies)), aes(start))+geom_ribbon(aes(ymin=0, ymax=log2, fill=species))+facet_grid(species ~ .)+theme_classic()+line+scale_fill_manual(values=colors, name="Species")+xaxis+yaxis+ggtitle(paste(outputPrefix, "Avg depth of coverage", sep=" ")))
}

