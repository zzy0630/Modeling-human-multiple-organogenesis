#!/usr/bin/perl -w 
# RNA-seq analysis: alignment 
# 
my $h_genome = 'hg38'; # human genome hg38 version reference
my $h_gtf = 'hg38.gtf'; # human genome hg38 annotation

my $m_genome = 'mm10'; # mouse mm10 version reference
my $m_gtf = 'mm10.gtf'; # mouse mm10 annotation

# hEA samples 
my @file qw(
Large-B-1
Large-B-2
Large-B-3
Large-BB-1
Large-BB-2
Large-BB-3
Large-C-1
Large-C-2
Large-C-3
Large-CB-1
Large-CB-2
Large-CB-3
Large-CC-1
Large-CC-2
Large-CC-3
Large-E8-1
Large-E8-2
Large-E8-3
Small-B-1
Small-B-2
Small-B-3
Small-BB-1
Small-BB-2
Small-BB-3
Small-C-1
Small-C-2
Small-C-3
Small-CB-1
Small-CB-2
Small-CB-3
Small-CC-1
Small-CC-2
Small-CC-3
Small-E8-1
Small-E8-2
Small-E8-3
);

foreach my $one (@file) {
	# HISAT2 (v.2.0.5) alignment 
	# paired-end
	system "hisat2 -p 8 --dta -3 2 -5 2 -x $h_genome -1 $one\-R1.fastq.gz -2 $one\-R2.fastq.gz -S $one.sam"; 
	
	# samtools (v.1.2) and Stringtie (v.1.3.3b) assembly
	system "samtools view -bS $one.sam > $one.bam";
	system "samtools sort -@ 8 $one.bam -o $one.srt.bam";
	system "stringtie -p 8 -G $h_gtf -e- o ./Stringtie_hEA/$one.gtf -l $one $one.srt.bam";
}


# Mouse temporal gastruloids 
# GSE106227
my @file = qw(
H120-1
H120-2
H144-1
H144-2
H168-1
H168-2
H24-1
H24-2
H48-1
H48-2
H72-1
H72-2
H96-1
H96-2
);

foreach my $one (@file) {
	# HISAT2 alignment
	# Single end
	system "hisat2 -p 8 --dta -3 2 -5 2 -x $m_genome -U $one.fastq.gz -S $one.sam"; 
	
	# samtools and Stringtie assembly
	system "samtools view -bS $one.sam > $one.bam";
	system "samtools sort -@ 8 $one.bam -o $one.srt.bam";
	system "stringtie -p 8 -G $m_gtf -e- o ./Stringtie_gastruloid/$one.gtf -l $one $one.srt.bam";
}

# Integrate into gene expression matrix for DEseq2
# use Stringtie provided tool: prep.py2.py
# Integrate hEA matrix
system "prep.py2.py -i ./Stringtie_hEA/ -g ./DESeq2_hEA/gene_count_matrix.csv -t ./DESeq2/transcript_count_matrix.csv";

# Integrate gastruloid matrix
system "prep.py2.py -i ./Stringtie_gastruloid/ -g ./DESeq2_gastruloid/gene_count_matrix.csv -t ./DESeq2/transcript_count_matrix.csv";





























