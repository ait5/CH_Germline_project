#! /usr/bin/env Rscript 



setwd('/work/dmp/franches/Germline_project/CH_Germline_project')

library(data.table)
library(dplyr)




## LOADING : (data frame with anonymized features from ~50k cancer patients)
	in_dir <- 'data_mined/'
	NAME <- paste(in_dir, 'CH.GERM_association.matrix___2022-09-08.tsv', sep='')
ano <- read.delim(NAME, header=T, sep='\t')



# df_ano <- ano %>% filter(fully_annotated==1) ## ATTENTION!!
df_ano <- ano

##
germ.genes.cols <- colnames(df_ano)[grep('gMutORCNV.', colnames(df_ano), fixed=T)]
	germ.genes.cols = unique(gsub('gMutORCNV.', '', germ.genes.cols))
	
	gMut.cols = paste0('g.', germ.genes.cols)
	gCNV.cols = paste0('cnv.', germ.genes.cols)
g_mut.OR.cnv_genes = paste0('gMutORCNV.', germ.genes.cols)
all.germ.cols = c(g_mut.OR.cnv_genes, gMut.cols, gCNV.cols)

ch.genes.cols <- colnames(df_ano)[grep('ch.DNMT3A', colnames(df_ano))[1]:grep('ch.HIST1H3E', colnames(df_ano))[1]]

mca_cols = paste0('mCA.chr', seq(1:23))
mca_genes = colnames(df_ano)[grep('mCA.gene.', colnames(df_ano), fixed=T)]

ch.OR.mCA_genes = colnames(df_ano)[grep('ch.OR.mCA.DNMT3A', colnames(df_ano))[1]:grep('ch.AND.mCA.HIST1H3E', colnames(df_ano))[1]]
ch.AND.mCA_genes = ch.OR.mCA_genes[grep('ch.AND.mCA.', ch.OR.mCA_genes)]
ch.OR.mCA_genes = ch.OR.mCA_genes[grep('ch.OR.mCA.', ch.OR.mCA_genes)]


two.Hits_genes <- colnames(df_ano)[grep('two.Hits.', colnames(df_ano))]

MajorCancerTypes_cols <- colnames(df_ano)[grep('Non.Small.Cell.Lung.Cancer', colnames(df_ano)):grep('Other.cancer.type', colnames(df_ano))]
##




## Output directory

out_dir <- 'Results/'
dir.create(out_dir)

out_dir_1 <- paste0(out_dir, 'Main_associations/')
dir.create(out_dir_1)

out_dir_2 <- paste0(out_dir_1, 'tables_main/')
	dir.create(out_dir_2)





set.seed(123)
########################################
######## Main associations:
data <- df_ano

# name.dir <- 'all'


v_conditions <- c('germ_event', 'germ_mutated', 'germ_cnv',
 'ch.pd_germ.gene', 'heme_genes_asco', 'DDR_germ.gene', 
	# 'Splicing_germ.gene', 'EpigModif_germ.gene', 
	'adult_MDS.AML', 'bone_marrow_failure', 'all', 'lymphoma', 'heme.malignancy')

v_dependents <- c(
	# 'ch_status', 'n_ch_muts', 'one_ch_muts', 'more_than_one_ch_muts',
  'functional_ch', 'one_ch_muts_functional', 'more_than_one_ch_muts_functional',
   'functional_ch_pd',
   # 'one_ch_PD_muts', 'more_than_one_ch_PD_muts',
    'ch_cnv', 'ch_amp', 'ch_del', 'ch_loh',
    'DDR_CH', 'Splicing_CH', 'EpigModif_CH', 'DNMT3A_CH', 'Other_CH', 'HighVAF_CH', 'functionalCH.OR.mCA', 'functionalCH.AND.mCA',
    paste0('functional_ch__', MajorCancerTypes_cols))

i=6
j=1

# l_conditions <- list()
# for (i in 1:length(v_conditions)){
	# CONDITION <- v_conditions[i]
	CONDITION <- as.character(commandArgs(6)[1])
	print(paste('##--', CONDITION, '--##', sep=''))

	data2 <- data 


	l_dependents <- list()
	for (j in 1:length(v_dependents)){
		DEPENDENT <- v_dependents[j]
		print(paste('--', DEPENDENT, '--', sep=''))


		m <- data2 %>% select(all_of(DEPENDENT), all_of(CONDITION)) %>% table

		Glm.test <- glm(as.factor(get(DEPENDENT)) ~ get(CONDITION)+age_bin+Sex+ancestry_label+smoking_bin+Major.Cancer.Type+binary_therapy , data = data2, family = "binomial")
		# b <- broom::tidy(Glm.test, exponentiate = T)[2,c('estimate', 'p.value')] %>%
		# 	as.data.frame()
		# 	colnames(b) <- c('glm.odds-ratio', 'glm.p.value')
		b <- broom::tidy(Glm.test, exponentiate = T, conf.int=T)[2,c('estimate', 'p.value', 'conf.low', 'conf.high')] %>%
			as.data.frame()
			colnames(b) <- c('glm.odds-ratio', 'glm.p.value', 'glm.conf.low', 'glm.conf.high')

		m2 <- c(CONDITION, DEPENDENT,
				 sum(m),
				  sum(m[,2]),
				  	sum(m[2,]),
					   m[2,2],
					    sum(m[2,])/sum(m),
					     m[2,2]/sum(m[,2])) %>%
			as.data.frame %>% t()
			colnames(m2) <- c('Condition', 'Dependent', 'Total', 'n_Condition', 'n_Dependent', 'n_Dependent_in_Condition', 'f_Dependent_in_Total', 'f_Dependent_in_Condition')

		ans <- as.data.frame( cbind(m2,b) )

		l_dependents[[paste(CONDITION, DEPENDENT, sep='-->')]] <- ans
	} # for (DEPENDENTS)
	ANS <- as.data.frame(do.call(rbind, l_dependents))
	# l_conditions[[CONDITION]] <- ANS
# } # for (CONDITIONS)


## Export:
	NAME <- paste(out_dir_2, 'MainAssociations_', CONDITION, '_', Sys.Date(), '.tsv', sep='')
	write.table(ANS, NAME, col.names=T, row.names=F, quote=F, sep='\t')

NAME 









