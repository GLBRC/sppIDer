__author__ = 'Quinn'

import sys, re, time

################################################################
# This script will take the sam formatted output of bwa and parse for mapping quality and which genomes reads map to.
#
# Input: sam file of reads mapped to a combination reference genome
################################################################


inputName = sys.argv[1]
samName = inputName+".sam"
outputName = inputName+"_MQ.txt"
start = time.time()

speciesDict = {}
MQscoreDict = {}
speciesDict["*"] = {}
speciesDict["*"][0] = 0
speciesList = ['*']
sam = open(samName, 'r')
samLines = sam.read().splitlines()
for line in samLines:
    if re.match('^(@SQ)', line):
        headerInfo = line.split('\t')
        chrInfo = headerInfo[1].split(":")[1]
        chrName = chrInfo.split("-")
        speciesName = chrName[0]
        chrNum = int(chrName[1])
        if chrNum==1:
            speciesList.append(speciesName)
            speciesDict[speciesName] = {}
            for i in range(0, 61):
                speciesDict[speciesName][i] = 0
    elif re.match('^(?!@)', line):
        lineSplit = line.split('\t')
        chr = lineSplit[2].split("-")
        species = chr[0]
        pos = lineSplit[3]
        MQscore = int(lineSplit[4])
        count = speciesDict[species][MQscore]
        speciesDict[species][MQscore] = count+1

output = open(outputName, 'w')
output.write("Species\tMQscore\tcount\n")
for species in speciesList:
    for score in speciesDict[species].keys():
        output.write(species+"\t"+str(score)+"\t"+str(speciesDict[species][score])+"\n")

currentTime = time.time()-start
print(str(currentTime)+" secs\n")
