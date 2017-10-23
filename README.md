# sppIDer

sppIDer is a pipeline for looking at genome composition in hybrid genomes and check for chromosomal copy variants in single species strains.  

sppIDer.py is the main wrapper that calls established bioinformatic tools and custom scripts. This pipeline needs a combination reference genome and one or more short read (fastq) files. 


### Usage

The sppIDer docker image is a self-contained platform capable of executing its pipeline without requiring cumbersome managment and installation of prerequisite tools.

Execute the container with a host volume mount, as shown below, to retrieve pipeline output files into the host machine's current working directory. Providing the example "--user" switch will 


### Examples

sppIDer help/syntax:  
```
docker run --rm -it glbrc/sppider -h
```

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


