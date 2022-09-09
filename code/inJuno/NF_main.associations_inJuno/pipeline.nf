#!/usr/bin/env nextflow


// Script parameters
wd = "/work/dmp/franches/Germline_project/CH_Germline_project"

df_genes = file("${wd}/code/inJuno/NF_main.associations_inJuno/all_germ_conditions.txt")
allNames  = df_genes.readLines()



process rscript {
	errorStrategy 'ignore'
	
	input:
  	val NAME from allNames

	"""
	Rscript ${wd}/code/inJuno/NF_main.associations_inJuno/script_main.associations.R ${NAME}
	"""
} 



// Run in terminal:
//% nextflow run pipeline.nf