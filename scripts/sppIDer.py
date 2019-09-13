__author__ = 'Quinn'

import argparse, multiprocessing, sys, re, subprocess, time, os

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
parser.add_argument('--delete', help="Delete BAM/SAM and Bedgraph files after completion, optional")
parser.add_argument('--byBP', help="Calculate coverage by basepair, optional, DEFAULT, can't be used with -byGroup", dest='bed', action='store_true')
parser.add_argument('--byGroup', help="Calculate coverage by chunks of same coverage, optional, can't be used with -byBP", dest='bed', action='store_false')
parser.set_defaults(bed=True)
args = parser.parse_args()

# docker vars
scriptDir = "/tmp/sppIDer/"
workingDir = "/tmp/sppIDer/working/"
numCores = str(multiprocessing.cpu_count())

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

trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'w')
trackerOut.write("outputPrefix="+args.out+"\n")
trackerOut.write("ref="+refGen+"\n")
trackerOut.write("read1=" + read1Name + "\n")
if args.r2: trackerOut.write("read2=" + read2Name + "\n")
if args.bed == False:
    trackerOut.write("coverage analysis option =  by coverage groups, bedgraph format -bga\n")
else: trackerOut.write("coverage analysis option = by each base pair -d\n")
trackerOut.close()

########################## BWA ###########################
bwaOutName = outputPrefix + ".sam"
bwaOutFile = open(workingDir + bwaOutName, 'w')
if args.r2:
    print("Read1=" + read1Name + "\nRead2=" + read2Name)
    subprocess.call(["bwa", "mem", "-t", numCores, refGen, read1Name, read2Name], stdout=bwaOutFile, cwd=workingDir)
else:
    print("Read1=" + read1Name)
    subprocess.call(["bwa", "mem", "-t", numCores, refGen, read1Name], stdout=bwaOutFile, cwd=workingDir)
bwaOutFile.close()
print("BWA complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("BWA complete\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## samtools ###########################
samViewOutQual = outputPrefix + ".view.bam"
bamSortOut = outputPrefix + ".sort.bam"
samViewQualFile = open(workingDir + samViewOutQual, 'w')
subprocess.call(["samtools", "view", "-@", numCores, "-q", "3", "-bhSu", bwaOutName], stdout=samViewQualFile, cwd=workingDir)
samViewQualFile.close()
subprocess.call(["samtools", "sort", "-@", numCores, samViewOutQual, "-o", bamSortOut], cwd=workingDir)
print("SAMTOOLS complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("\nSAMTOOLS complete\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## parse SAM file ###########################
subprocess.call(["python2.7", scriptDir + "parseSamFile.py", outputPrefix], cwd=workingDir)
print("Parsed SAM file")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("\nParsed SAM\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## plot MQ scores ###########################
subprocess.call(["Rscript", scriptDir + "MQscores_sumPlot.R", outputPrefix], cwd=workingDir)
print("Plotted MQ scores")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("\nMQ scores plotted\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## bedgraph Coverage ###########################
sortOut = bamSortOut
if args.bed == True:
    bedOutD = outputPrefix + "-d.bedgraph"
    bedFileD = open(workingDir + bedOutD, 'w')
    subprocess.call(["genomeCoverageBed", "-d", "-ibam", sortOut], stdout=bedFileD, cwd=workingDir)
    bedFileD.close()
else:
    bedOut = outputPrefix + ".bedgraph"
    bedFile = open(workingDir + bedOut, 'w')
    subprocess.call(["genomeCoverageBed", "-bga", "-ibam", sortOut], stdout=bedFile, cwd=workingDir)
    bedFile.close()
print("bedgraph complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("\nbedgraph complete\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## average Bed ###########################
if args.bed == True:
    subprocess.call(["Rscript", scriptDir + "meanDepth_sppIDer-d.R", outputPrefix], cwd=workingDir)
else:
    subprocess.call(["Rscript", scriptDir + "meanDepth_sppIDer-bga.R", outputPrefix], cwd=workingDir)
print("Found mean depth")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("\nFound mean depth\nElapsed time: " + elapsedTime)
trackerOut.close()

########################## make plot ###########################
subprocess.call(["Rscript", scriptDir + "sppIDer_depthPlot_forSpc.R", outputPrefix], cwd=workingDir)
subprocess.call(["Rscript", scriptDir + "sppIDer_depthPlot.R", outputPrefix], cwd=workingDir)
# if args.bed == True:
#     subprocess.call(["Rscript", scriptDir + "sppIDer_depthPlot_forSpc.R", outputPrefix, "d"], cwd=workingDir)
#     subprocess.call(["Rscript", scriptDir + "sppIDer_depthPlot-d.R", outputPrefix], cwd=workingDir)
# else:
#     subprocess.call(["Rscript", scriptDir + "sppIDer_depthPlot_forSpc.R", outputPrefix, "g"], cwd=workingDir)
#     subprocess.call(["Rscript", scriptDir + "sppIDer_depthPlot-bga.R", outputPrefix], cwd=workingDir)
print("Plot complete")
currentTime = time.time()-start
elapsedTime = calcElapsedTime(currentTime)
print("Elapsed time: " + elapsedTime)
trackerOut = open(workingDir + outputPrefix + "_sppIDerRun.info", 'a')
trackerOut.write("\nPlot complete\nElapsed time: " + elapsedTime + "\n")
trackerOut.close()

if args.delete
   print("Cleaning up working files")
   os.remove(outputPrefix + ".sam")
   os.remove(outputPrefix + ".sort.bam")
   os.remove(outputPrefix + ".view.bam")
   os.remove(outputPrefix + "-d.bedgraph")

