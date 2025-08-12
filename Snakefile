rule all:
	input:
		"output/DESEQ2_results.csv",
		"output/DESEQ2_results_sig.csv"

rule differential_expression:
	input:
		zero_counts="data/ZERO_RNAseq/ZERO_GeneExpression_rawCounts_CNS_NB_Sarcoma_mRNA_24052024.txt",
		normal_counts="data/PsychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",
		metadata="data/various/ZERO_psychencode_metadata_average.xlsx"
	output:
		all="output/DESEQ2_results.csv",
		sig="output/DESEQ2_results_sig.csv"
	script:
		"src/differential_expression.R"

