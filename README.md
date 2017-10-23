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

### Examples

Notes:  
- Execute the container with a host volume mount, as shown below, to retrieve pipeline output files into the host machine's current working directory  
- Providing the example "--user" switch will write to output files using permissions of the host user  

Execute pipeline using full paths for clarity:  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  --out GLBRC \
  --ref /tmp/sppIDer/inputs/SaccharomycesCombo.fasta \
  --r1 /tmp/sppIDer/inputs/R1_1k.fastq \
  --r2 /tmp/sppIDer/inputs/R2_1k.fastq
```

Run in the container's working directory, for input brevity:  
```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/sppIDer \
--user "$UID:$(id -g $USERNAME)" \
glbrc/sppider \
  --out GLBRC \
  --ref SaccharomycesCombo.fasta \
  --r1 R1_1k.fastq \
  --r2 R2_1k.fastq
```

### System Requirements 

This pipeline has been tested on a CentOS 7.4 host running Docker CE Stable.


