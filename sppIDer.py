__author__ = 'Quinn'

import sys, re, subprocess, time, argparse

################################################################
# This script runs the full sppIDer pipeline.
# Before running this you must make your combination reference genome with the script combineRefGenomes.py
# This script will map short-read data to a combination reference genome and parse the outputs to create a summary of
# where and how well the reads map to each species in the combinaiton reference genome.
#
#Program Requirements: bwa, samtools, bedtools, R, Rpackage ggplot2, Rpackage dplyr
#Input: Output name, Combination reference genome, fastq short-read files
#
################################################################

parser = argparse.ArgumentParser(description="Run full sppIDer")
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--ref', help="Reference Genome, required", required=True)
parser.add_argument('--r1', help="Read1, required", required=True)
parser.add_argument('--r2', help="Read2, optional")
parser.add_argument('--byBP', help="Calculate coverage by basepair, optional, DEFAULT, can't be used with -byGroup", dest='bed', action='store_true')
parser.add_argument('--byGroup', help="Calculate coverage by chunks of same coverage, optional, can't be used with -byBP", dest='bed', action='store_false')
parser.set_defaults(bed=True)
args = parser.parse_args()

#Replace paths

outputPrefix = args.out
refGen=args.ref
read1Name = args.r1
if args.r2: read2Name = args.r2
start = time.time()
def calcElapsedTime( endTime ):
    trackedTime = str()
    if 60 < endTime < 3600:
        min = int(endTime) / 60
        sec = int(endTime - (min * 60))
        trackedTime = "%s mins %s secs" % (min, sec)
    elif 3600 < endTime < 86400:
        hr = int(endTime) / 3600
        min = int((endTime - (hr * 3600)) / 60)
        sec = int(endTime - ((hr * 60) * 60 + (min * 60)))
        trackedTime = "%s hrs %s mins %s secs" % (hr, min, sec)
    elif 86400 < endTime < 604800:
        day = int(endTime) / 86400
        hr = int((endTime - (day * 86400)) / 3600)
        min = int((endTime - (hr * 3600 + day * 86400)) / 60)
        sec = int(endTime - ((day * 86400) + (hr * 3600) + (min * 60)))
        trackedTime = "%s days %s hrs %s mins %s secs" % (day, hr, min, sec)
    elif 604800 < endTime:
        week = int(endTime) / 604800
        day = int((endTime)-(week * 604800) / 86400)
        hr = int((endTime - (day * 86400 + week * 604800)) / 3600)
        min = int((endTime - (hr * 3600 + day * 86400 + week * 604800)) / 60)
        sec = int(endTime - ((week * 604800) + (day * 86400) + (hr * 3600) + (min * 60)))
        trackedTime = "%s weeks %s days %s hrs %s mins %s secs" % (week, day, hr, min, sec)
    else:
        trackedTime = str(int(endTime)) + " secs"
    return trackedTime

trackerOut = open(outputPrefix+"_sppIDerRun.info", 'w')
trackerOut.write("outputPrefix="+args.out+"\n")
trackerOut.write("ref="+refGen+"\n")
trackerOut.write("read1=" + read1Name + "\n")
if args.r2: trackerOut.write("read2=" + read2Name + "\n")
if args.bed == False:
    trackerOut.write("coverage analysis option =  by coverage groups, bedgraph format -bga\n")
else: trackerOut.write("coverage analysis option = by each base pair -d\n")
trackerOut.close()

########################## BWA ###########################
bwaOutName = outputPrefix+"_aln-pe.sam"
bwaOutFile = open(bwaOutName, 'w')
if args.r2:
    print("Read1=" + read1Name + "\nRead2=" + read2Name)
    #trackerOut.write("read1=" + read1Name + "\n")
    #trackerOut.write("read2=" + read2Name + "\n")
    #trackerOut.close()
    subprocess.call(["bwa", "mem", refGen, read1Name, read2Name], stdout=bwaOutFile)
else:
    #read1Name = args.r1
    print("Read1="+read1Name)
    #trackerOut.write("read1=" + read1Name + "\n")
    #trackerOut.close()
    subprocess.call(["bwa", "mem", refGen, read1Name], stdout=bwaOutFile)
print("BWA complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("BWA complete\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## samtools ###########################
#samViewOut = outputPrefix+"_aln-pe.view.bam"
samViewOutQual = outputPrefix+"_aln-pe.view.bam"
#samSortOut = outputPrefix+".sorted.sam"
bamSortOut = outputPrefix+"_aln-pe.sort.bam"
#samViewFile = open(samViewOut, 'w')
samViewQualFile = open(samViewOutQual, 'w')
subprocess.call(["samtools", "view", "-q", "3", "-bhSu",  bwaOutName], stdout=samViewQualFile)
subprocess.call(["samtools", "sort", samViewOutQual, "-o", bamSortOut])
print("SAMTOOLS complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("\nSAMTOOLS complete\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## parse SAM file ###########################
#parseInput = outputPrefix+"_aln-pe"
subprocess.call(["python2.7", "scripts/parseSamFile.py", outputPrefix])
print("Parsed SAM file")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("\nParsed SAM\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## plot MQ scores ###########################
subprocess.call(["Rscript", "scripts/MQscores_sumPlot.R", outputPrefix])
print("Plotted MQ scores")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("\nMQ scores plotted\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## bedgraph Coverage ###########################
sortOut = bamSortOut
if args.bed == True:
    bedOutD = outputPrefix + "_coverage-d.bedgraph"
    bedFileD = open(bedOutD, 'w')
    subprocess.call(["genomeCoverageBed", "-d", "-ibam", sortOut], stdout=bedFileD)
else:
    bedOut = outputPrefix+"_coverage.bedgraph"
    bedFile = open(bedOut, 'w')
    subprocess.call(["genomeCoverageBed", "-bga", "-ibam", sortOut], stdout=bedFile)
print("bedgraph complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("\nbedgraph complete\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## average Bed ###########################
if args.bed == True:
    subprocess.call(["Rscript", "scripts/meanDepth_sppIDer-d.R", outputPrefix])
else:
    subprocess.call(["Rscript", "scripts/meanDepth_sppIDer-bga.R", outputPrefix])
print("Found mean depth")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("\nFound mean depth\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## make plot ###########################
if args.bed == True:
    subprocess.call(["Rscript", "scripts/sppIDer_depthPlot_forSpc.R", outputPrefix, "d"])
    subprocess.call(["Rscript", "scripts/sppIDer_depthPlot-d.R", outputPrefix])
else:
    subprocess.call(["Rscript", "scripts/sppIDer_depthPlot_forSpc.R", outputPrefix, "g"])
    subprocess.call(["Rscript", "scripts/sppIDer_depthPlot-bga.R", outputPrefix])
print("Plot complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(outputPrefix+"_sppIDerRun.info", 'a')
trackerOut.write("\nPlot complete\nElapsed time: " + elapsedTime + "\n")
trackerOut.close()

