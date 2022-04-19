#!/usr/bin/perl -w 
# Day 20 hOM single cell RNA-seq alignment 
# 
# Merge fastq files
system "cat hOM_part_1_R1.fastq.gz hOM_part_2_R1.fastq.gz hOM_part_3_R1.fastq.gz hOM_part_4_R1.fastq.gz > hOM_part_1234_R1.fastq.gz";
system "cat hOM_part_1_R2.fastq.gz hOM_part_2_R2.fastq.gz hOM_part_3_R2.fastq.gz hOM_part_4_R2.fastq.gz > hOM_part_1234_R2.fastq.gz";

# Alignment: Cellranger (v.6.0.1)
system "cellranger count --id=MidOrganoid --fastqs=./fastqs --sample=hOM_part_1234 --transcriptome=./refdata-gex-GRCh38-2020-A"; 
# /refdata-gex-GRCh38-2020-A is downloaded from 10x official website

# Prepare LOOM file by velocyto (v.0.17.17)
system "velocyto run -o ./velocyto -@ 16 ./MidOrganoid/outs/possorted_genome_bam.bam ./refdata-gex-GRCh38-2020-A/genes/genes.gtf"

