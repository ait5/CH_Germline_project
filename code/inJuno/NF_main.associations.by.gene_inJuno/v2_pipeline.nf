#!/usr/bin/env nextflow


// Script parameters
wd = "/work/dmp/franches/Germline_project/CH_Germline_project"

df_genes = file("${wd}/code/inJuno/NF_main.associations.by.gene_inJuno/all_germ_cols.txt")
allNames  = df_genes.readLines()



process rscript {
	errorStrategy 'ignore'
	
	input:
  	val NAME from allNames

	"""
	Rscript ${wd}/code/inJuno/NF_main.associations.by.gene_inJuno/vNF_script_BY.germ.GENE_by_ch.gene.R ${NAME}
	"""
} 



// Run in terminal:
//% nextflow run v2_pipeline.nf