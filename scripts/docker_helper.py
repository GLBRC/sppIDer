#!/usr/bin/env python

class colors:
    pink = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    end = '\033[0m'
    bold = '\033[1m'
    underline = '\033[4m'
    reverse = '\033[07m'

def main():
    print ''
    print colors.blue + colors.bold + 'sppIDer bioinformatics pipeline' + colors.end
    print colors.blue + '_______________________________' + colors.end
    print ''
    print 'sppIDer is a pipeline for looking at genome composition in hybrid genomes and checking for chromosomal copy variants in single species strains.'
    print ''
    print 'sppIDer.py is the main wrapper that calls established bioinformatic tools and custom scripts. this pipeline needs a combination reference genome and one or more short read (fastq) files.'
    print ''
    print colors.blue + 'pipeline help/syntax:' + colors.end
    print ''
    print 'docker run --rm -it glbrc/sppider [pipeline_script] --help'
    print ''
    print 'available pipeline scripts:'
    print '  combineGFF.py'
    print '  combineRefGenomes.py'
    print '  sppIDer.py'
    print '  mitoSppIDer.py'
    print ''
    print colors.blue + 'sppIDer.py usage example:' + colors.end
    print ''
    print 'this command will take a script\'s input files from your current working directory, outputting resultant pipeline files there upon completion.'
    print ''
    print 'docker run \ '
    print '--rm -it \ '
    print '--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \ '
    print '--user "$UID:$(id -g $USERNAME)" \ '
    print 'glbrc/sppider \ '
    print '  sppIDer.py \ '
    print '  --out OUT \ '
    print '  --ref REF.fasta \ '
    print '  --r1 R1.fastq \ '
    print '  --r2 R2.fastq'
    print ''


if __name__ == '__main__':
    main()

