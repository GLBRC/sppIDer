# sppIDer

sppIDer is a pipeline for looking at genome composition in hybrid genomes and check for chromosomal copy variants in single species strains.  

sppIDer.py is the main wrapper that calls established bioinformatic tools and custom scripts. This pipeline needs a combination reference genome and one or more short read (fastq) files. 

_The sppIDer docker image is a self-contained platform capable of executing its pipeline without requiring cumbersome managment and installation of prerequisite tools._   

_Changes to this source repo are automatically built into an updated docker image, available from docker hub at [glbrc/sppider](https://hub.docker.com/r/glbrc/sppider/)._


### Getting Started

##### pipeline help/syntax:  
```
docker run --rm -it glbrc/sppider [pipeline_script] --help

    pipeline scripts:
      sppIDer.py
      mitoSppIDer.py
      combineRefGenomes.py
```

##### example: sppIDer.py help  
```
docker run --rm -it glbrc/sppider sppIDer.py -h

    usage: sppIDer.py [-h] --out OUT --ref REF --r1 R1 [--r2 R2] [--byBP]
                      [--byGroup]
    
    Run full sppIDer
    
    optional arguments:
      -h, --help  show this help message and exit
      --out OUT   Output prefix, required
      --ref REF   Reference Genome, required
      --r1 R1     Read1, required
      --r2 R2     Read2, optional
      --byBP      Calculate coverage by basepair, optional, DEFAULT, can't be used
                  with -byGroup
      --byGroup   Calculate coverage by chunks of same coverage, optional, can't
                  be used with -byBP
```

### Pipeline Usage

Workflow:
- The combination reference genome must be built first using combineRefGenomes.py. The outputs can be used many times with sppIDer.py with different data sets.
- The main pipeline, sppIDer.py, takes fastq input(s) and maps the reads to the combined reference genome made with combineRefGenomes.py.
- The pipeline then uses bioinfromatic tools and custom scripts to pares this output for where, how well, and how deeply the reads map to combined reference genome by species, chromosomes, and windows. 
- The output is several pdfs with plots of precentage and quality of reads mapped and plots for coverage by species and in windows. Addionally several summary text files are created.
- All files are kept from intermediate steps and could be used in other anlyses. 

Notes:  
- Execute the container with a host volume mount, as shown below, to retrieve pipeline output files into the host machine's current working directory  
- Providing the example "--user" switch will write to output files using permissions of the host user  


##### example: executing a combineRefGenome.py  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  combineRefGenomes.py
  --out REF.fasta \ 
  --key KEY.txt
```
An optional --trim can be used to trim short uninformative contigs for reference genomes with many short contigs. All contigs shorter than the supplied interger will be ignored.
The KEY.txt file must be tab delimited and the reference genome unique name cannot contain hyphens. See example.

##### example: executing a sppIDer.py pipeline  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  sppIDer.py \
  --out OUT \
  --ref REF.fasta \
  --r1 R1.fastq \
  --r2 R2.fastq
```
An optional --byGroup flag can be used for very large combination genomes. This produce a bedfile that doesn't have coverage information for each basepair but by groups. Which speeds up the run.

###### For mitoSppIDer
##### example: executing combineGFF.py  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  combineRefGenomes.py
  --out REF.gff \ 
  --key GFF_KEY.txt
```

##### example: executing a mitoSppIDer.py pipeline  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  mitoSppIDer.py \
  --out OUT \
  --ref MITO_REF.fasta \
  --r1 R1.fastq \
  --r2 R2.fastq
```
An optional --gff can be used if you are providing a combined gff of the regions that should be marked on the final plots.


### System Requirements 

This pipeline has been tested on a CentOS 7.4 (1708) host running [Docker Community Edition (CE) Stable](https://docs.docker.com/engine/installation/) (17.09.0.ce).

