#!/usr/bin/env python
import sys, subprocess, argparse
from Bio import SeqIO

################################################################
# This script will create the combination reference genome and all additional files needed to run the full sppIDer script.
#
#Program Requirements: bwa, samtools
#Input: output name, text file key of genome names, and desired minimum length of contigs included
#
################################################################

parser = argparse.ArgumentParser(description="Combine desired reference genomes")
parser.add_argument('--out', help="Output prefix, required", required=True)
parser.add_argument('--key', help="Key to reference genomes, required", required=True)
parser.add_argument('--trim', help="Maximum contig length to trim", default=0)
args = parser.parse_args()
comboGenomeName = args.out
listName = args.key
trimLength = args.trim

#Replace paths

lengthFile = open("comboLength_"+comboGenomeName+".txt", 'w')
lengthFile.write(comboGenomeName+"\tTrimmed contigs >"+str(trimLength)+"\n")
outGenome = open(comboGenomeName, 'w')
list = open(listName, 'r')
lines = list.readlines()
for line in lines:
    line = line.strip().split('\t')
    uniID = line[0]
    genomeName = line[1]
    sumGenomeLen = 0
    counter = 0
    fasta = open(genomeName, 'r')
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if len(seq_record.seq)>=trimLength:
            outGenome.write(">"+uniID+"-"+str(counter+1)+"\n")
            outGenome.write(str(seq_record.seq)+"\n")
            lengthFile.write(uniID+"-"+str(counter)+"\t"+str(len(seq_record.seq))+"\n")
            sumGenomeLen += len(seq_record.seq)
            counter += 1
    strGenomeLen = ''
    if sumGenomeLen > 1000000000:
        genomeLen = sumGenomeLen / 1000000000
        strGenomeLen = str(genomeLen) + " Gb"
    elif sumGenomeLen>1000000:
        genomeLen = sumGenomeLen/1000000
        strGenomeLen = str(genomeLen) + " Mb"
    elif sumGenomeLen>1000:
        genomeLen = sumGenomeLen/1000
        strGenomeLen = str(genomeLen)+" Kb"
    else:
        strGenomeLen = str(sumGenomeLen) + "bp"
    lengthFile.write(uniID + "-totalGenome\t" + strGenomeLen + "\n")

outGenome.close()
lengthFile.close()

subprocess.call(["/opt/bifxapps/bin/bwa", "index", comboGenomeName])
subprocess.call(["/opt/bifxapps/bin/samtools", "faidx", comboGenomeName])
