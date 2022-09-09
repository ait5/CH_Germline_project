#! /usr/bin/env Rscript 

setwd('/work/dmp/franches/Germline_project/CH_Germline_project')

library(data.table)
library(dplyr)
library(doBy)


	

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





###############################################################################
##############################################################################

TERM <- 'CH'
term_col <- paste(TERM, '.Rate_when_Gene.germ.mut', sep='')

out_dir <- 'Results/'
dir.create(out_dir)

out_dir_1 <- paste(out_dir, 'tables_gene_by_gene/', sep='')
dir.create(out_dir_1)

data <- df_ano


somatic_genes = c(ch.genes.cols, mca_cols,
 ch.OR.mCA_genes, ch.AND.mCA_genes,
    paste0('functional_ch__', MajorCancerTypes_cols), mca_genes)


##
	germ.gene <- commandArgs(6)[1] #parsing $NAME from bash#
	print(paste('Germ gene is', germ.gene))


		stats_by_gene <- lapply(somatic_genes, function(x){
			gene <- x
			print(gene)


			m <- as.matrix(table(data[,c(gene, germ.gene)]))
				m.0.0 <- data %>% filter(get(gene)==0 & get(germ.gene)==0) %>% nrow()
				m.1.0 <- data %>% filter(get(gene)==1 & get(germ.gene)==0) %>% nrow()
				m.0.1 <- data %>% filter(get(gene)==0 & get(germ.gene)==1) %>% nrow()
				m.1.1 <- data %>% filter(get(gene)==1 & get(germ.gene)==1) %>% nrow()

				m1 <- matrix(c(m.0.0, m.1.0, m.0.1, m.1.1), nrow = 2,
				dimnames = list(gene = c(0,1),
                            germ.gene = c(0,1)))
				x.Rate <- sum(m.0.1, m.1.1)/sum(m.0.0, m.1.0, m.0.1, m.1.1)
				y.Rate_when_x <- m.1.1/sum(m.0.1, m.1.1)

			f.test <- fisher.test(m1)
			ans_fisher <- c(as.character(f.test$estimate),
				as.character(f.test$p.value),
				as.character(f.test$conf.int[1]),
				as.character(f.test$conf.int[2]))
			names(ans_fisher) <- c(
					'fisher.odds.ratio', 'fisher.pvalue',
					'fisher.odds.ratio.MIN',
					'fisher.odds.ratio.MAX')



				ans <- c(gene, m.0.0, m.1.0, m.0.1, m.1.1, rep(NA, 2), x.Rate, y.Rate_when_x)
				names(ans) <- c('Gene',
					paste(germ.gene, '.NO__Gene.NO', sep=''),
					paste(germ.gene, '.NO__Gene.YES', sep=''),
					paste(germ.gene, '.YES__Gene.NO', sep=''),
					paste(germ.gene, '.YES__Gene.YES', sep=''),
					'glm.odds.ratio', 'glm.pvalue',
					'Gene.germ.mut.Rate', term_col)


			if ( gene %in% paste0('functional_ch__', MajorCancerTypes_cols)) {

				### In case any confounder is totally biased (e.g with GENDER being M=0)
				v_biases <- c()
				for (i in c('Sex','smoking_bin', 'binary_therapy')){
					table_factor <- data %>% filter(get(germ.gene)==1) %>% select(all_of(i)) %>% table()
					if (length(grep(TRUE, as.vector(table_factor==0))>0)){
						v_biases <- c(v_biases, TRUE)
							print(paste(i, 'is biased!'))
					} else {
						v_biases <- c(v_biases, FALSE)
					}
				}
				names(v_biases) <- c('Sex','smoking_bin', 'binary_therapy')


				if (m.0.0>0 & m.1.0>0 & m.0.1>0 & sum(m.1.0, m.1.1)>3 & sum(m.0.1, m.1.1)>3 & sum(m.0.1, m.1.0, m.1.1)>30 & length(grep(TRUE, v_biases))==0){ # There has to be mutations in four groups (contingency table)

					Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+Sex+ancestry_label+smoking_bin+binary_therapy , data = data, family = "binomial")

					tb.tidy <- broom::tidy(Glm.test, exponentiate = T)[2,c('estimate', 'p.value')] %>%
						as.data.frame()
						colnames(tb.tidy) <- c('estimate', 'p.value')
				# cancer.type.gene.tidy

				ans <- c(gene, m.0.0, m.1.0, m.0.1, m.1.1,
					as.character(tb.tidy[,'estimate']),
					as.character(tb.tidy[,'p.value']),
					x.Rate, y.Rate_when_x)

					names(ans) <- c('Gene',
						paste(germ.gene, '.NO__Gene.NO', sep=''),
						paste(germ.gene, '.NO__Gene.YES', sep=''),
						paste(germ.gene, '.YES__Gene.NO', sep=''),
						paste(germ.gene, '.YES__Gene.YES', sep=''),
						'glm.odds.ratio', 'glm.pvalue',
						'Gene.germ.mut.Rate', term_col)

				} else if (m.0.0>0 & m.1.0>0 & m.0.1>0 & sum(m.1.0, m.1.1)>3 & sum(m.0.1, m.1.1)>3 & length(grep(TRUE, v_biases))>0) {
					co.term <- names(v_biases[grep(TRUE, v_biases)])
					if (co.term=='Sex'){
						Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+ancestry_label+smoking_bin+binary_therapy , data = data, family = "binomial")
						} else if (co.term=='smoking_bin') {
							Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+Sex+ancestry_label+binary_therapy , data = data, family = "binomial")
						}  else if (co.term=='binary_therapy') {
							Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+Sex+ancestry_label+smoking_bin , data = data, family = "binomial")
							}
					tb.tidy <- broom::tidy(Glm.test, exponentiate = T)[2,c('estimate', 'p.value')] %>%
						as.data.frame()
						colnames(tb.tidy) <- c('estimate', 'p.value')

					ans <- c(gene, m.0.0, m.1.0, m.0.1, m.1.1,
					as.character(tb.tidy[,'estimate']),
					as.character(tb.tidy[,'p.value']),
					x.Rate, y.Rate_when_x)

					names(ans) <- c('Gene',
						paste(germ.gene, '.NO__Gene.NO', sep=''),
						paste(germ.gene, '.NO__Gene.YES', sep=''),
						paste(germ.gene, '.YES__Gene.NO', sep=''),
						paste(germ.gene, '.YES__Gene.YES', sep=''),
						'glm.odds.ratio', 'glm.pvalue',
						'Gene.germ.mut.Rate', term_col)
				}

			} else {

					### In case any confounder is totally biased (e.g with GENDER being M=0)
				v_biases <- c()
				for (i in c('Sex','smoking_bin', 'binary_therapy')){
					table_factor <- data %>% filter(get(germ.gene)==1) %>% select(all_of(i)) %>% table()
					if (length(grep(TRUE, as.vector(table_factor==0))>0)){
						v_biases <- c(v_biases, TRUE)
							print(paste(i, 'is biased!'))
					} else {
						v_biases <- c(v_biases, FALSE)
					}
				}
				names(v_biases) <- c('Sex','smoking_bin', 'binary_therapy')


				if (m.0.0>0 & m.1.0>0 & m.0.1>0 & sum(m.1.0, m.1.1)>3 & sum(m.0.1, m.1.1)>3 & sum(m.0.1, m.1.0, m.1.1)>30 & length(grep(TRUE, v_biases))==0){ # There has to be mutations in four groups (contingency table)

					Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+Sex+ancestry_label+smoking_bin+Major.Cancer.Type+binary_therapy , data = data, family = "binomial")

					tb.tidy <- broom::tidy(Glm.test, exponentiate = T)[2,c('estimate', 'p.value')] %>%
						as.data.frame()
						colnames(tb.tidy) <- c('estimate', 'p.value')
				# cancer.type.gene.tidy

				ans <- c(gene, m.0.0, m.1.0, m.0.1, m.1.1,
					as.character(tb.tidy[,'estimate']),
					as.character(tb.tidy[,'p.value']),
					x.Rate, y.Rate_when_x)

					names(ans) <- c('Gene',
						paste(germ.gene, '.NO__Gene.NO', sep=''),
						paste(germ.gene, '.NO__Gene.YES', sep=''),
						paste(germ.gene, '.YES__Gene.NO', sep=''),
						paste(germ.gene, '.YES__Gene.YES', sep=''),
						'glm.odds.ratio', 'glm.pvalue',
						'Gene.germ.mut.Rate', term_col)

				} else if (m.0.0>0 & m.1.0>0 & m.0.1>0 & sum(m.1.0, m.1.1)>3 & sum(m.0.1, m.1.1)>3 & length(grep(TRUE, v_biases))>0) {
					co.term <- names(v_biases[grep(TRUE, v_biases)])
					if (co.term=='Sex'){
						Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+ancestry_label+smoking_bin+Major.Cancer.Type+binary_therapy , data = data, family = "binomial")
						} else if (co.term=='smoking_bin') {
							Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+Sex+ancestry_label+Major.Cancer.Type+binary_therapy , data = data, family = "binomial")
						}  else if (co.term=='binary_therapy') {
							Glm.test <- glm(get(gene) ~ get(germ.gene)+age_bin+Sex+ancestry_label+smoking_bin+Major.Cancer.Type , data = data, family = "binomial")
							}
					tb.tidy <- broom::tidy(Glm.test, exponentiate = T)[2,c('estimate', 'p.value')] %>%
						as.data.frame()
						colnames(tb.tidy) <- c('estimate', 'p.value')

					ans <- c(gene, m.0.0, m.1.0, m.0.1, m.1.1,
					as.character(tb.tidy[,'estimate']),
					as.character(tb.tidy[,'p.value']),
					x.Rate, y.Rate_when_x)

					names(ans) <- c('Gene',
						paste(germ.gene, '.NO__Gene.NO', sep=''),
						paste(germ.gene, '.NO__Gene.YES', sep=''),
						paste(germ.gene, '.YES__Gene.NO', sep=''),
						paste(germ.gene, '.YES__Gene.YES', sep=''),
						'glm.odds.ratio', 'glm.pvalue',
						'Gene.germ.mut.Rate', term_col)
				}

			}


			

			c(ans, ans_fisher)

		}) #end lapply
		names(stats_by_gene) <- somatic_genes

	xx <- as.data.frame(do.call(rbind,stats_by_gene))
	xx$glm.pvalue <- as.numeric(as.character(xx$glm.pvalue))
	
	xx2 <- as.data.frame(orderBy(~glm.pvalue+glm.odds.ratio, xx))
	rownames(xx2) <- 1:nrow(xx2)
	xx2$glm.odds.ratio <- as.numeric(as.character(xx2$glm.odds.ratio))
	xx2$glm.pvalue <- as.numeric(as.character(xx2$glm.pvalue))

	## Adjusting Pvalue for Multiple Testing:
	xx2$glm.FDR.adj.Pvalue <- p.adjust(xx2$glm.pvalue, method="fdr")
	xx2$fisher.FDR.adj.Pvalue <- p.adjust(xx2$fisher.pvalue, method="fdr")

	NAME <- paste(out_dir_1, germ.gene, '_ass.with_CH.genes.tsv', sep='')
	write.table(xx2, NAME, sep='\t', col.names=T, row.names=F, quote=F)


