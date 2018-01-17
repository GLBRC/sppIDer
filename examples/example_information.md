# Example Information for sppIDer  
The end of the [sppIDer manual](../sppIDerManual.md) had how to run the sppIDer pipeline with these examples.  

## Example reference genome input  
Files in [exampleCombineRef](exampleCombineRef.tar.gz)  
SaccharomycesGenomesKey.txt - Tab separated text file to match reference genome with desired unique name.  
Reference Genomes:  
File Name | Species | Source | Publication  
--------- | ------- | ------ | -----------  
Scer.fasta | *Saccharomyces cerevisiae* | [Saccharomyces Sensu Stricto](http://www.saccharomycessensustricto.org/cgi-bin/s3.cgi?data=Assemblies&version=current) | [SGD](https://www.yeastgenome.org/)  
Spar.fasta | *Saccharomyces paradoxus* | [Saccharomyces Sensu Stricto](http://www.saccharomycessensustricto.org/cgi-bin/s3.cgi?data=Assemblies&version=current) | Liti and Carter et al. 2009 Nature  
Smik.fasta | *Saccharomyces mikatae* | [Saccharomyces Sensu Stricto](http://www.saccharomycessensustricto.org/cgi-bin/s3.cgi?data=Assemblies&version=current) | Scannell and Zill et al. 2011 G3  
SkudZP.fasta | *Saccharomyces kudriavzeii* | [Saccharomyces Sensu Stricto](http://www.saccharomycessensustricto.org/cgi-bin/s3.cgi?data=Assemblies&version=current) | Scannell and Zill et al. 2011 G3  
Sarb.fasta | *Saccharomyces arboricola* | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA_000292725.1) | Liti et al. 2013 BMC Genomics  
Suva.fasta | *Saccharomyces uvarum* | [Saccharomyces Sensu Stricto](http://www.saccharomycessensustricto.org/cgi-bin/s3.cgi?data=Assemblies&version=current) | Scannell and Zill et al. 2011 G3  
Seub_wMito.fasta | *Saccharomyces eubayanus* | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA\_001298625.1) | Baker et al. 2015 MBE  

The published *S. uvarum* genome (Scannell and Zill et al 2011) had chromosome X swapped with chromosome XII, which is fixed in this version.  

####References:  
Baker E, Wang B, Bellora N, Peris D, Hulfachor AB, Koshalek JA, Adams M, Libkind D, Hittinger CT. 2015. The genome sequence of *Saccharomyces eubayanus* and the domestication of lager-brewing yeasts. Mol. Biol. Evol. 32:2818–2831.  
Liti G, Carter DM, Moses AM, Warringer J, Parts L, James SA, Davey RP, Roberts IN, Burt A, Koufopanou V, et al. 2009. Population genomics of domestic and wild yeasts. Nature 458:337–341.  
Liti G, Nguyen Ba AN, Blythe M, Müller CA, Bergström A, Cubillos FA, Dafhnis-Calas F, Khoshraftar S, Malla S, Mehta N, et al. 2013. High quality de novo sequencing and assembly of the *Saccharomyces arboricolus* genome. BMC Genomics 14:1-14.  
Scannell DR, Zill OA, Rokas A, Payen C, Dunham MJ, Eisen MB, Rine J, Johnston M, Hittinger CT. 2011. The Awesome Power of Yeast Evolutionary Genetics: New Genome Sequences and Strain Resources for the *Saccharomyces sensu stricto* Genus. G3 1:11–25.  

## Example test data  
Files in [exampleFastqs](exampleFastq.tar.gz)  
SRA | Species | Strain | Publication  
--- | ------- | ------ | -----------  
[SRR2586160](https://www.ncbi.nlm.nih.gov/sra/SRR2586160/) | *Saccharomyces eubayanus* | yHRVM108 | Peris and Langdon et al. 2016 PLoS Genet.  
[SRR2586169](https://www.ncbi.nlm.nih.gov/sra/SRR2586169/) | *Saccharomyces cerevisiae X S. eubayanus* | Weihenstephan 34/70 (syn. yHAB47) | Peris and Langdon et al. 2016 PLoS Genet.  
[SRR1119201](https://www.ncbi.nlm.nih.gov/sra/SRR1119201/) | *Saccharomyces cerevisiae X S. kudriavzevii X S. uvarum X S. eubayanus* | CBS2834 | Almeida et al. 2014 Nature communications  

####References: 
Almeida P, Gonçalves C, Teixeira S, Libkind D, Bontrager M, Masneuf-Pomarède I, Albertin W, Durrens P, Sherman DJ, Marullo P, et al. 2014. A Gondwanan imprint on global diversity and domestication of wine and cider yeast *Saccharomyces uvarum*. Nat. Commun. 5:4044.  
Peris D, Langdon QK, Moriarty R V, Sylvester K, Bontrager M, Charron G, Leducq J, Landry CR, Libkind D, Hittinger CT. 2016. Complex Ancestries of Lager-Brewing Hybrids Were Shaped by Standing Variation in the Wild Yeast *Saccharomyces eubayanus*. PLoS Genet. 12: e1006155.  

## Example outputs  
Files in [exampleOutputs](exampleOutputs.tar.gz). Each set of input fastqs have corresponding output files.   
File Suffix | Description  
----------- | -----------  
SRR\*\_sppIDerRun.info | A text file that contains the options and inputs for that run and the time to run.  
SRR\*\_MQsummary.txt | Text file with summary of how many and how well reads map to each genome.  
SRR\*\_plotMQ.pdf | Plot of reads mapped per genome and Mapping Quality per genome.  
SRR\*\_speciesAvgDepth-d.txt | Text file summary of coverage for each species including: mean, relativeMean (speciesMean/globalMean), max, and median coverage.  
SRR\*\_speciesDepth.pdf | Plot of coverage by species.  
SRR\*\_sppIDerDepthPlot-d.pdf | Plot of coverage by genome split into 10,000 windows.  

## Example mitochondrial inputs  
Files in [exampleMito](exampleMito.tar.gz)  
mitoRefKey.txt - Tab separated text file to match mitochondrial (mito) reference genome with desired unique name.  
mitoGFFKey.txt - Tab separated text file to match mitochondrial (mito) reference GFF with desired unique name.  
Species Files:  
Species | Reference | GFF | Accession | Publication | Genes reannotated  
------- | --------- | --- | --------- | ----------- | -----------------  
*S. cerevisiae* | S288c\_mtDNA.fasta | S288c\_mtDNA.gff | [NC_001224](https://www.ncbi.nlm.nih.gov/nuccore/NC_001224) | Foury et al. 1998 FEBS Lett. | None  
*S. paradoxus* | CBS432\_mtDNA.fasta | CBS432\_mtDNA.gff | [NC_018044](https://www.ncbi.nlm.nih.gov/nuccore/NC_018044) | Prochazka et al. 2012 FEMS Yeast Res. | RF2  
*S. mikatae* | KX707788\_mtDNA.fasta | KX707788\_mtDNA.gff | [KX707788](https://www.ncbi.nlm.nih.gov/nuccore/KX707788) | NA | RF2, F-SmikIII,VAR1,COX1  
*S. kudriavzeii* | KX707787-Skud\_mtDNA.fasta | KX707787-Skud\_mtDNA.gff | [KX707787](https://www.ncbi.nlm.nih.gov/nuccore/KX707787) | NA | F-SkudIII  
*S. arboricola* | CBS10644\_mtDNA.fasta | CBS10644\_mtDNA.gff | [KX657740](https://www.ncbi.nlm.nih.gov/nuccore/KX657740) | Sulo et al 2017 | COX1,RF3,F-SarbIII,COB,VAR1,COX3,ATP6,COX2,ATP9,ATP8  
*S. uvarum* | CBS395\_mtDNA.fasta | CBS395\_mtDNA.gff | [KX657742](https://www.ncbi.nlm.nih.gov/nuccore/KX657742) | Sulo et al 2017 | COX1,COB,F-SuvaIII,VAR1,COX3,ATP6,COX2,ATP9,ATP8  
*S. eubayanus* | FM1318\_mtDNA.fasta | FM1318\_mtDNA.gff | [NW_017264706.1](https://www.ncbi.nlm.nih.gov/nuccore/NW_017264706.1) | Baker et al. 2015 MBE | None  

Genes reannotated: Genes that are not in the published gff but were added here.  

####References:  
Baker E, Wang B, Bellora N, Peris D, Hulfachor AB, Koshalek JA, Adams M, Libkind D, Hittinger CT. 2015. The genome sequence of *Saccharomyces eubayanus* and the domestication of lager-brewing yeasts. Mol. Biol. Evol. 32:2818–2831.  
Foury F, Roganti T, Lecrenier N, Purnelle B. 1998. The complete sequence of the mitochondrial genome of *Saccharomyces cerevisiae*. FEBS Lett. 440:325–331.  
Procházka E, Franko F, Poláková S, Sulo P. 2012. A complete sequence of *Saccharomyces paradoxus* mitochondrial genome that restores the respiration in *S. cerevisiae*. FEMS Yeast Res. 12:819–830.  
Sulo P, Szabóová D, Bielik P, Polákova S, Šoltys K, Jatzová K, Szemes T. 2017. The evolutionary history of *Saccharomyces* species inferred from completed mitochondrial genomes and revision in the 'yeast mitochondrial genetic code'. DNA Res. 24:571-583.  