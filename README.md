# sppIDer

sppIDer is a pipeline for looking at genome composition in hybrid genomes and check for chromosomal copy variants in single species strains.  

sppIDer.py is the main wrapper that calls established bioinformatic tools and custom scripts. This pipeline needs a combination reference genome and one or more short read (fastq) files. 

### usage example

```
docker run \
--rm -it \
--mount type=bind,src=$(pwd),target=/tmp/working \
glbrc/sppider \
--out GLBRC --ref /tmp/working/comboGenome/SaccharomycesCombo.fasta --r1 /tmp/working/R1_1k.fastq --r2 /tmp/working/R2_1k.fastq

```

