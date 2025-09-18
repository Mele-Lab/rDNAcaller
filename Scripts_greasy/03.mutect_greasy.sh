#!/bin/bash

# load modules
module load java-openjdk gatk bcftools
#module load java-openjdk/17.0.11+9
#module load gatk/4.5.0.0
       
sample=$1
input_folder=$2
output_folder=$3
species=$4
threads=112

bam=${input_folder}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam

if [[ $species == "mouse" ]]; then
	reference=/gpfs/projects/bsc83/Data/assemblies/Mm39/Mouse_mm39-rDNA_genome_v1.0/mm39-rDNA_v1.0.fa
	intervals=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/11_Mice/data/pre-rRNA_47S.regions.bed
elif [[ $species == "human" ]]; then
	reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
	#The reference comes from this paper: https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext 
	#intervals=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.included_intervals.bed
	intervals=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.regions.bed
else
	echo "Species not available. Try mouse or human"
fi

echo "$bam has started at $(date)"

gatk --java-options "-Xmx115g -Xms100g" Mutect2 \
	-I ${bam} \
	-R ${reference} \
	-O ${output_folder}/${sample}.mutect_1.g.vcf.gz \
	--native-pair-hmm-threads $(nproc) \
	-L $intervals

bcftools index -t -f ${output_folder}/${sample}.mutect_1.g.vcf.gz	


#Removing PL from the vcf, this is necessary for bcftools norm
bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${output_folder}/${sample}.mutect_1.g.vcf.gz > ${output_folder}/${sample}_filtering.g_wo_PL.vcf

rm ${output_folder}/${sample}.mutect_1.g.vcf.gz ${output_folder}/${sample}.mutect_1.g.vcf.gz.tbi ${output_folder}/${sample}.mutect_1.g.vcf.gz.stats 

#Remove alternative alleles that have not been called from g.vcf. Sometime GATK adds a variant in the ALT column when it has not been called
#-O z is to get the output as .gz
bcftools view --trim-alt-alleles ${output_folder}/${sample}_filtering.g_wo_PL.vcf -O z -o ${output_folder}/${sample}_trim.vcf.gz

rm ${output_folder}/${sample}_filtering.g_wo_PL.vcf

#Split multiallelic positions into different rows of "biallelic" positions. If we add the reference fasta, reference alleles are adjusted so all donors call them in the same way
#-O z is to get the output as .gz.
bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref $reference \
        ${output_folder}/${sample}_trim.vcf.gz -O z --threads ${threads} > ${output_folder}/${sample}.vcf.gz

bcftools index -t ${output_folder}/${sample}.vcf.gz #Creating index in .tbi format (-t), we will need index if we plan to use bcftools merge
        	
rm ${output_folder}/${sample}_trim.vcf.gz

echo "$bam has finished at $(date)"
