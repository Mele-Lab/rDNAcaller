# rDNAcaller
Repository to call rDNA variants from short read WGS

# Software dependencies

| Package | Version | 
| -------- | ------- | 
| MUMMER | 4.0.0rc1 | 
| ANACONDA | 2023.07 | 
| seqkit | 2.8.0 |
| Samtools | 1.19.2 |
| bwa | 0.7.17 |
| gatk | 4.5.0.0 |
| java-openjdk | 17.0.11+9|
| bcftools | 1.19 |
| R | 4.3.2 |

We created a singularity image with all necessary software
```
singularity pull library://jmiguelramirez/jmiguelramirez/rdna_variant_caller:latest
```


# Running the code
Adding reference data for the code to work
```
#Reference fasta file:
curl -L -o Reference/Human_hs1-rDNA_genome_v1.0.tar.gz https://github.com/vikramparalkar/rDNA-Mapping-Genomes/raw/main/Human_hs1-rDNA_genome_v1.0.tar.gz
#Untar
tar -xvzf Reference/Human_hs1-rDNA_genome_v1.0.tar.gz
#Create .dict file
gatk-launch CreateSequenceDictionary -R Reference/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa

#Create bwa index
bwa index -p Reference/new_index Reference/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
```

To run the code sequentially
```
./Pipeline.sh -h
This code runs a variant caller for rDNA reads

Syntax: Pipeline.sh [-h|-t|-n|-f|-i|-o]
parameters:
h     Print the help message.
t     Type of analysis: Either WGS to call all variants or vcf to call a subset of known variants.
n     Basename of the input data: (e.g., for a bam file named sample.bam: -n sample, for fastqs named sample_1.fastq.gz and sample_2.fasq.gz: -n sample. This means that paired fastqs should end by _1.fastq.gz and _2.fastq.gz)
f     Type of input data format: Either bam (for bam/cram/sam) or fastq.
i     Input folder.
o     Output folder.
```

Example using a small fraction of reads from a CHM13 sample, for which input fastq files and expected output vcf are in Data/
```
./Pipeline.sh -t WGS -n SRR3189741 -f fastq -i Data/ -o Results/
```

If using the singularity image:
```
singularity exec rdna_variant_caller_latest.sif ./Pipeline.sh -t WGS -n SRR14773542 -f fastq -i Data/ -o Results/
```

To run a parallel version of the code using SLURM, check inside Scripts_parallel or using greasy, Scripts_greasy/

# Running complementary code

Other code used to complement the manuscript: "rDNAcaller: a fast and robust workflow to call ribosomal DNA variants" is available inside Downstream_analysis/ such as the code to compute mean nucleotide diversity score

Code to run the simulator from "rDNAcaller: a fast and robust workflow to call ribosomal DNA variants" is available here: https://github.com/Rbbt-Workflows/SyntheticRDNA and the singularity image here:
```
singularity pull --arch amd64 library://jmiguelramirez/jmiguelramirez/simulator:latest
```



