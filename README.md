# MEeR-seq
Pipeline for mapping the 5' ends and 3' ends of pre-rRNA precursors as a proxy for assessing ribosome biogenesis.

## Follow the detailed script called merr.seq.pipeline.sh, which is summarized below

1.	Preprocess raw sequencing reads using Cutadapt to remove adapter sequences.
2.	Extract unique molecular identifiers (UMIs).
3.	Prepare a custom rRNA reference using the RNA45SN5 reference sequence (accession no. NT_167214.1). Save the resulting reference sequences in FASTA format as rRNA_processing.fa and generate an index using Bowtie2.
4.	Align the preprocessed reads to the indexed custom rRNA reference using Bowtie2. Sort and index the resulting alignment files using SAMtools.
5.	Remove PCR duplicates using UMICollapse.
6.	Split the deduplicated alignments into first-in-pair (R1) and second-in-pair (R2) using SAMtools. To retain the exact 5'-end of rRNA precursors, retain R1 from the 5'-end libraries. Conversely, to retain the exact 3'-end of rRNA precursors, retain R2 from the 3'-end libraries. 
7.	Convert the retained R1 and R2 alignments into BigWig coverage files using deepTools. Normalize the coverage to reads per kilobase per million mapped reads (RPKM) to facilitate comparisons across samples.

## Analyze the bigwig files using the file bigwig_analysis.Rmd as reference	

