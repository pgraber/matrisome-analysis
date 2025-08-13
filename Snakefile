rule all:
	input:
		"output/DESEQ2_results.csv",
		"output/DESEQ2_results_sig.csv",
		expand("output/enrichment/GOplot_chord_{ontology}.pdf", ontology=["MF", "BP", "CC"])

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

rule enrichment_analysis:
	input:
		deseq_results="output/DESEQ2_results.csv",
		gprofiler_data="data/gProfiler/gProfiler_hsapiens_12-8-2025_9-23-12 pm__intersections_{ontology}.csv"
	output:
		chord_plot="output/enrichment/GOplot_chord_{ontology}.pdf",
		circ_data="output/enrichment/circ_data_{ontology}.csv"
	script:
		"src/enrichment_analysis.R"

