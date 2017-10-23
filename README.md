# sppIDer

sppIDer is a pipeline for looking at genome composition in hybrid genomes and check for chromosomal copy variants in single species strains.  

sppIDer.py is the main wrapper that calls established bioinformatic tools and custom scripts. This pipeline needs a combination reference genome and one or more short read (fastq) files. 

_The sppIDer docker image is a self-contained platform capable of executing its pipeline without requiring cumbersome managment and installation of prerequisite tools._  


### Usage

sppIDer help/syntax:  
```
docker run --rm -it glbrc/sppider -h

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

### Pipeline Examples

Notes:  
- Execute the container with a host volume mount, as shown below, to retrieve pipeline output files into the host machine's current working directory  
- Providing the example "--user" switch will write to output files using permissions of the host user  

Execute a sppIDer pipeline using full paths for clarity:  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  --out OUT \
  --ref /tmp/sppIDer/working/inputs/REF.fasta \
  --r1 /tmp/sppIDer/working/inputs/R1.fastq \
  --r2 /tmp/sppIDer/working/inputs/R2.fastq
```

Run using relative paths of input files, for brevity:  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  --out OUT \
  --ref REF.fasta \
  --r1 R1.fastq \
  --r2 R2.fastq
```


### Generate Pipeline Input Files

_Use of the included combineRefGenomes.py script produces an output that can be used multiple times with different data. You may need to first execute this script to prepare files for your pipeline. We can override the default sppIDer image behavior and do this as follows._  

Usage info for combining desired reference genomes:
```
docker run --rm -it --entrypoint=/usr/bin/python2.7 glbrc/sppider /tmp/sppIDer/scripts/combineRefGenomes.py -h

	usage: combineRefGenomes.py [-h] --out OUT --key KEY [--trim TRIM]

	Combine desired reference genomes

	optional arguments:
	  -h, --help   show this help message and exit
	  --out OUT    Output prefix, required
	  --key KEY    Key to reference genomes, required
	  --trim TRIM  Maximum contig lenght to trim
```

Generate input files for subsequent sppIDer runs:
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer/working \
--user "$UID:$(id -g $USERNAME)" \
--entrypoint=/usr/bin/python2.7 \
glbrc/sppider \
/tmp/sppIDer/scripts/combineRefGenomes.py \
  --out OUT \
  --key KEY \
  --trim TRIM
```


### System Requirements 

This pipeline has been tested on a CentOS 7.4 host running Docker CE Stable.


