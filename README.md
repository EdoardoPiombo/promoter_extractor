# Requirements
To run the program you need python3, with the modules pandas, argparse and Bio.

# What does it do?
This is a program to add promoter coordinates to genes in a gff file. It will produce a new gff file with genes and promoter coordinates.
To generate the promoter coordinates an arbitrary number is subtracted to the start of the gene (for genes on the plus strand) or added to the end (for genes on the minus strand).
There is no default added length to generate promoters, this number should be chosen according to the existing literature on the organism you are working with.
Added promoters will not overlap existing genes. If two genes are already overlapping in the input gff file, then a token promoter of 1 bp will be added in front of the gene.
When two promoters would be generated in a way that makes them overlap (example: a gene on the plus strand immediately following one on the minus strand) the distance between the two genes is equally divided between their promoters, without overlapping. This will result in promoters shorter than the requested length.

# Input files
Promoter extractor requires a gff file with gene coordinates and a fasta file with the contigs indicated in the gff file. This fasta file is used to count the bases in each contig and make sure that no added promoter will go beyond the length of the contig.
