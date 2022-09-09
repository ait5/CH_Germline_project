

library(data.table)
library(doBy)
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
  'DDR_CH', 'Splicing_CH', 'EpigModif_CH', 'DNMT3A_CH', 'Other_CH', 'HighVAF_CH',
  'functionalCH.OR.mCA', 'functionalCH.AND.mCA',
  paste0('ch.OR.mCA.', MajorCancerTypes_cols))

i=6
j=1
in_dir = paste0(out_dir_1, 'tables_main/')
DATE = '2022-09-08' #@@@@ ATTENTION!!

l_conditions <- list()
for (i in 1:length(v_conditions)){
  CONDITION <- v_conditions[i]
  print(paste('##--', CONDITION, '--##', sep=''))
  
  if( length(grep(CONDITION, list.files(in_dir)))==1 ){
    name_file = paste0(in_dir, 'MainAssociations_', CONDITION, '_', DATE, '.tsv', sep='')
    x_file = read.delim(name_file, sep='\t', header=T)
  } else {
    paste0(CONDITION, " not found!") 
    next
  }
  
  ANS <- x_file
  l_conditions[[CONDITION]] <- ANS
} # for (CONDITIONS)

df <- as.data.frame(do.call(rbind, l_conditions))

MAIN <- df %>% mutate(glm.signif=fifelse(round(glm.p.value, 2)<=0.001, '***',
                                         fifelse(round(glm.p.value, 2)<=0.01, '**', 
                                                 fifelse(round(glm.p.value, 2)<=0.05, '*', 
                                                         fifelse(round(glm.p.value, 2)<0.1, 't', ''))))) 


##### Adjusting Pvalue for Multiple Testing (between all genes by Dependent):
dd = MAIN
dd <- as.data.frame(do.call(rbind, lapply(
  as.character(unique(dd$Dependent)), function(dependent_term){
    tt = dd %>% filter(Dependent == dependent_term)
    tt$Adj.glm.p.value <- p.adjust(tt$glm.p.value, method="fdr")
    tt = tt %>% 
      mutate(Adj.glm.signif = fifelse(
        round(Adj.glm.p.value,2)<=0.001, '***',
        fifelse(round(Adj.glm.p.value, 2)<=0.01, '**',
                fifelse(round(Adj.glm.p.value, 2)<=0.05, '*',
                        fifelse(round(Adj.glm.p.value, 2)<0.1, 't', '')))))
    tt
  })))
MAIN = dd
######

## Export:
NAME <- paste(out_dir_1, 'MainAssociations_', Sys.Date(), '.tsv', sep='')
write.table(MAIN, NAME, col.names=T, row.names=F, quote=F, sep='\t')

NAME 









