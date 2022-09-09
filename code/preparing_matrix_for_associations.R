#####
#####


library(data.table)
library(dplyr)
library(tidyverse)


out_dir <- 'data_mined/'
dir.create(out_dir)

in_dir = './data/'


## LOADING DATA:

## Load CH Impact calls
## CH variants for 50k  MSK-IMPACT (anonymized)
name <- 'FILTERED_50k.CH.calls_2022-02-16_anonymized.tsv_2022-08-23.txt'
NAME <- paste(in_dir, name, sep='')
impact_ch <- read.delim(NAME, sep='\t', header=T) %>% 
  arrange(DMP_PATIENT_ID) %>% unique()
dim(impact_ch)


## Clinical annotation by sample
name <- 'ch.50k_clinical_annots_by_SAMPLE_2021-12-06_anonymized.tsv'
NAME <- paste(in_dir, name, sep='')
clinical_by_sample <- read.delim(NAME) %>%
  unique %>% arrange(DMP_PATIENT_ID, Age_dpn.dob) %>%
  select(-PartC.Constent.Status) %>%
  rename(Teng.ch_cnv = ch_cnv) %>%
  rename(Teng.ch_amp = ch_amp) %>%
  rename(Teng.ch_del = ch_del) %>%
  rename(Teng.ch_loh = ch_loh) 
dim(clinical_by_sample)


## Oncotree conversion codes:
NAME <- paste(in_dir, 'oncotree_conversion_codes.txt', sep='')
df_oncotree <- read.delim(NAME, header=T, sep='\t')
dim(df_oncotree)


## Germline pathogenic calls reviewed by Miika for 90 genes:
name <- 'DUPLI_filtered_Cohort_50k_Pathogenic_90Genes_v2.txt'
NAME <- paste(in_dir, name, sep='')
germ_calls <- read.delim(NAME, header=T, sep='\t') %>% 
  separate('Normal_Sample',
           into=c('P', 'Number', 'NT', 'ImpactVersion'), sep='-', remove=F) %>%
  mutate(DMP_PATIENT_ID = paste(P, Number, sep='-')) %>%
  select(-c('P', 'Number', 'NT')) %>% arrange(DMP_PATIENT_ID) %>% unique()
dim(germ_calls)
# Just created DMP_PATIENT_ID column 

## List of patients assessed for germline variants (only anonymized):
name <- 'List_of_anonymized_Samples.txt'
NAME <- paste(in_dir, name, sep='')
germ_patients_d <- read.delim(NAME, header=T, sep='\t') 
dim(germ_patients_d)
# Just created DMP_PATIENT_ID column 

## (updated) mCA calling by CH-FACETS:
name <- 'mCA.50k.annotations_2022-03-07_anonymized.tsv'
NAME <- paste(in_dir, name, sep='')
mCA_df <- read.delim(NAME, header=T, sep='\t') %>%
  unique
dim(mCA_df)

## (updated) mCA CALLS by CH-FACETS:
name <- 'calls_mCA.50k.annotations_2022-03-07_anonymized.tsv'
NAME <- paste(in_dir, name, sep='')
mCA_calls <- read.delim(NAME, header=T, sep='\t') %>%
  unique
dim(mCA_calls)



## Pathogenic CNVs [only deletions in TSG]:
name <- 'Cohort_50k_CNVs_90Genes_deletions_only_pathogenic_only_v1.txt'
NAME <- paste(in_dir, name, sep='')
germ_cnv_df <- read.delim(NAME, header=T, sep='\t') %>%
  rename(DMP_ASSAY_ID_Normal = sample) %>%
  separate('DMP_ASSAY_ID_Normal',
           into=c('P', 'Number', 'N-T', 'ImpactVersion'), sep='-', remove=F) %>%
  mutate(DMP_PATIENT_ID = paste(P, Number, sep='-')) %>%
  select(-c('P', 'Number')) %>% arrange(DMP_PATIENT_ID) %>% unique
dim(germ_cnv_df)

## Genes in hg19 genome:
name <- 'refGene_hg19.txt'
NAME <- paste(in_dir, name, sep='')
hg19_genes <- read.delim(NAME, sep='\t', header=T) %>% as.data.frame()


#############################################################################
#########################                    ################################
#########################     BY PATIENT     ################################
#########################                    ################################
#############################################################################

## Selecting earliest Dateof_Procedure in clinical_by_sample to select one Oncotree_Code, Age
dim(clinical_by_sample)
data_1 = clinical_by_sample %>%
  arrange(DMP_PATIENT_ID, Age_dpn.dob) %>%
  mutate(dupli = duplicated(DMP_PATIENT_ID)) %>%
  filter(dupli==FALSE) %>% select(-dupli) # Checking patients with multiple age annotations... And keeping the youngest!!
dim(data_1)

data_2 = data_1 %>%
  rename(age=Age_dpn.dob) %>%
  rename(age_bin=Age_bin) %>% 
  select(-MRN) %>%
  mutate(smoking_bin=fifelse(
    grepl('Current / Prior', Smoking_Hx), 1, ifelse(
      grepl('Never', Smoking_Hx), 0, NA
    )
  ))
dim(data_2)


##################################################################
###################### CH MUTATIONS ######################

calls_1 = impact_ch %>% 
  mutate(functional_ch = ifelse(Variant_Classification %in%
                                  c("Frame_Shift_Del", "Frame_Shift_Ins",
                                    "In_Frame_Del", "In_Frame_Ins",
                                    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site"), 1, 0))

## Annotating CH STATUS:
dummy_ch_status = calls_1 %>%
  mutate(ch_status = 1) %>%
  mutate(ch_pd = ifelse(cbp_driver_annotation == 'CH-PD', 1, 0)) %>%
  mutate(functional_ch = ifelse(Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site"), 1, 0)) %>%
  mutate(functional_ch_pd = ifelse(functional_ch == 1 & ch_pd == 1, 1, 0)) %>%
  select(DMP_PATIENT_ID, ch_status, ch_pd, functional_ch, functional_ch_pd) 

ch.status_pts <- dummy_ch_status %>% 
  filter(ch_status == 1) %>%
  select(DMP_PATIENT_ID) %>%
  unique %>%  pull %>% as.character
ch.pd_pts <- dummy_ch_status %>% 
  filter(ch_pd == 1) %>%
  select(DMP_PATIENT_ID) %>% pull %>% as.character

functional.ch_pts <- dummy_ch_status %>% 
  filter(functional_ch == 1) %>%
  select(DMP_PATIENT_ID) %>% 
  unique %>% pull %>% as.character
functional.ch.pd_pts <- dummy_ch_status %>% 
  filter(functional_ch_pd == 1) %>%
  select(DMP_PATIENT_ID) %>% 
  unique %>%pull %>% as.character

data_3 = data_2

data_4 <- data_3 %>%
  # mutate(ch_status = ifelse(DMP_PATIENT_ID %in% ch.status_pts, 1, 0)) %>%
  # mutate(ch_pd = ifelse(DMP_PATIENT_ID %in% ch.pd_pts, 1, 0)) %>%
  mutate(functional_ch = ifelse(DMP_PATIENT_ID %in% functional.ch_pts, 1, 0)) %>%
  mutate(functional_ch_pd = ifelse(DMP_PATIENT_ID %in% functional.ch.pd_pts, 1, 0))
dim(data_4)

## Adding number of mutations:
total_ch_muts__df = dummy_ch_status %>%
  filter(ch_status == 1) %>%
  count(DMP_PATIENT_ID) %>% rename(total_ch_muts = n)
total_ch.pd_muts__df = dummy_ch_status %>%
  filter(ch_pd == 1) %>%
  count(DMP_PATIENT_ID) %>% rename(total_ch.pd_muts = n)		
total_functional.ch_muts__df = dummy_ch_status %>%
  filter(functional_ch == 1) %>%
  count(DMP_PATIENT_ID) %>% rename(total_ch.functional_muts = n)
total_functional.ch.pd_muts__df = dummy_ch_status %>%
  filter(functional_ch_pd == 1) %>%
  count(DMP_PATIENT_ID) %>% rename(total_functional.ch.pd_muts = n)

data_5 <- data_4 %>%
  left_join(total_ch_muts__df, by='DMP_PATIENT_ID') %>%
  mutate(total_ch_muts = ifelse(is.na(total_ch_muts), 0, total_ch_muts)) %>%
  left_join(total_ch.pd_muts__df, by='DMP_PATIENT_ID') %>%
  mutate(total_ch.pd_muts = ifelse(is.na(total_ch.pd_muts), 0, total_ch.pd_muts)) %>%
  left_join(total_functional.ch_muts__df, by='DMP_PATIENT_ID') %>%
  mutate(total_ch.functional_muts = ifelse(is.na(total_ch.functional_muts), 0, total_ch.functional_muts)) %>%
  left_join(total_functional.ch.pd_muts__df, by='DMP_PATIENT_ID') %>%
  mutate(total_functional.ch.pd_muts = ifelse(is.na(total_functional.ch.pd_muts), 0, total_functional.ch.pd_muts))
dim(data_5)

# Checking...
#
# data_5 %>% select(ch_status, ch_pd, functional_ch, functional_ch_pd, total_ch_muts, total_ch.pd_muts, total_ch.functional_muts, total_functional.ch.pd_muts) %>% apply(2, table, useNA='always')


## Annotating binary tags for CH mutated genes:
ch_mutated_genes <- calls_1 %>% select(Hugo_Symbol) %>% unique %>% pull %>% as.character

for (i in 1:length(ch_mutated_genes)) {
  gene <- ch_mutated_genes[i]
  print(gene)
  
  patients_with_mut.gene <- calls_1 %>%
    filter(Hugo_Symbol==gene & functional_ch==1) %>% select(DMP_PATIENT_ID) %>% unique() %>% pull() %>% as.character()
  data_5 <- data_5 %>% mutate(gene=fifelse(DMP_PATIENT_ID%in%patients_with_mut.gene, 1, 0))
  colnames(data_5)[ncol(data_5)] <- paste0('ch.', gene)
}
data_6 <- data_5
dim(data_6)


##################################################################
###################### Extra columns CH features ######################
##
data_7 = data_6 %>%
  mutate(n_ch_muts=ifelse(total_ch_muts>=6, '>=6',total_ch_muts)) %>%
  mutate(one_ch_muts=fifelse(n_ch_muts==1, 1, 0)) %>%
  mutate(more_than_one_ch_muts=fifelse(n_ch_muts!=0 &
                                         n_ch_muts!=1, 1, 0)) %>%
  mutate(one_ch_muts_functional=fifelse(total_ch.functional_muts==1, 1, 0)) %>%
  mutate(more_than_one_ch_muts_functional=fifelse(total_ch.functional_muts>1, 1, 0)) %>%
  mutate(one_ch_PD_muts=fifelse(total_functional.ch.pd_muts==1, 1, 0)) %>%
  mutate(more_than_one_ch_PD_muts=fifelse(total_functional.ch.pd_muts>1, 1, 0))




## Annotationg 'mainType' of cancer
y <- df_oncotree %>% mutate(OncoTree_Code=code) %>% select(OncoTree_Code, mainType)

data_ano <- left_join(data_7, y)
data_ano$mainType <- as.character(data_ano$mainType)
data_ano[which(is.na(data_ano$mainType)), 'mainType'] <- 'Cancer of Unknown Primary'

data_7 = data_ano %>% 
  mutate(mainType = ifelse(mainType == 'Non-Small Cell Lung Cancer',
                           'Non Small Cell Lung Cancer', mainType))
dim(data_7)





## ## ## ## 
## Annotating mCA by patient:
## ## ## ## 

## Removing mCA calls affecting germline CNV deletions locations:
## Annotating affected genes on mCA_calls
mCA_calls$genes_affected = NA
mCA_calls$n_genes = NA
mCA_calls$cnv.conflict = FALSE

mCA_df.2 = mCA_df %>% select(1:4)
mCA_df.2$mCAs = as.character(mCA_df.2$mCAs)

germ_cnv_df$cnv.conflict = FALSE

for (i in 1:nrow(mCA_calls)){
  # i = 50
  print(paste0(i, ' from ', nrow(mCA_calls) ))
  CHROM = paste0('chr', mCA_calls[i, 'chrom'])
  if(CHROM == 'chr23') { CHROM = 'chrX'}
  START = mCA_calls[i, 'start']
  END = mCA_calls[i, 'end']
  PATIENT = mCA_calls[i, 'DMP_PATIENT_ID']
  CALL = mCA_calls[i, 'event_name']
  
  genes = hg19_genes %>% filter(chrom==CHROM & 
                                  ( 
                                    (start<START & end>START) |
                                      (start<END & end>END) |
                                      (start>=START & end<=END)	
                                  )
  ) %>% select(symbol_name) %>%
    unique %>% pull %>% as.character
  n_genes = length(genes)
  mCA_calls[i, 'n_genes'] = n_genes
  
  if (length(genes)>0){
    mCA_calls[i, 'genes_affected'] = as.character(
      paste(genes, collapse=';'))
  }
  
  call_in_cnv = germ_cnv_df %>%
    filter(DMP_PATIENT_ID == as.character(PATIENT) & Gene %in% genes)
  
  if (nrow(call_in_cnv)>0){
    mCA_calls[i, 'cnv.conflict'] = TRUE
    germ_cnv_df[which(germ_cnv_df$DMP_PATIENT_ID == as.character(PATIENT) & germ_cnv_df$Gene %in% genes), 'cnv.conflict'] = TRUE
  } else { mCA_calls[i, 'cnv.conflict'] = FALSE }
  
  # if ( mCA_calls[i, 'cnv.conflict']==TRUE ) {
  # 	mCA_term = as.character(mCA_df.2[which(mCA_df.2$DMP_PATIENT_ID==as.character(PATIENT)),'mCAs'])
  
  # 	mCA_terms = strsplit(
  # 		mCA_df.2 %>%
  # 		filter(DMP_PATIENT_ID == as.character(PATIENT)) %>%
  # 		 select(mCAs) %>% pull %>% as.character, split=';')[[1]]
  # 	mCA_terms_NEW = as.character(paste(mCA_terms[!mCA_terms%in%CALL], collapse=';'))
  # 	print(mCA_terms_NEW)
  
  # 	# mCA_df.2 = mCA_df.2 %>%
  # 	# 	mutate(mCAs_new= ifelse(DMP_PATIENT_ID == as.character(PATIENT), mCA_terms_NEW, mCAs))
  # 	mCA_df.2[which(mCA_df.2$DMP_PATIENT_ID==as.character(PATIENT)), 'mCAs'] = as.character(mCA_terms_NEW)
  
  
  # }
}
# mCA_df.2 %>% filter(DMP_PATIENT_ID=='A-894c080433eb')
# mCA_df.2 %>% filter(DMP_PATIENT_ID=='P-0056891')
# mCA_df.2 %>% filter(DMP_PATIENT_ID=='A-3f50c5a730ca') 

# germ_cnv_df %>% filter(DMP_PATIENT_ID=='A-894c080433eb')

############## EXPORTING ############## 
## Exporting 'mCA_calls' with germ.conflict annotations:
NAME = paste0(out_dir, 'mCA_calls_and_germ.conflict_', Sys.Date(),
              '.txt')
write.table(mCA_calls, NAME, col.names=T, row.names=F, quote=F, sep='\t')

## Exporting 'germ_cnv_df' with germ.conflict annotations:
NAME = paste0(out_dir, 'germ_cnvs_and_mCA.conflict_', Sys.Date(),
              '.txt')
write.table(germ_cnv_df, NAME, col.names=T, row.names=F, quote=F, sep='\t')
##########################################


## Creating big matrix by sample and annotating general categories ('ch_cnv', 'ch_del', 'ch_amp', 'ch_loh')
print('Creating big matrix by sample...')

mm = mCA_df.2 %>%
  mutate(mCAs=ifelse(mCAs=='', 'none', mCAs)) %>%
  mutate(ch_cnv = ifelse(mCAs!='none', 1, 0)) %>%
  mutate(ch_del = ifelse(grepl('del', mCAs), 1, 0)) %>%
  mutate(ch_amp = ifelse(grepl('amp', mCAs), 1, 0)) %>%
  mutate(ch_loh = ifelse(grepl('loh', mCAs), 1, 0))

## By mCA event (different chromosomes)
print('By mCA event (different chromosomes)...')

all_events = c(
  paste0('mCA.chr', seq(1, 23, by=1), '_'),
  paste0('mCA.chr', seq(1, 23, by=1), '_amp'),
  paste0('mCA.chr', seq(1, 23, by=1), '_del'),
  paste0('mCA.chr', seq(1, 23, by=1), '_loh')
)

mm2 = mm
for (j in 1:length(all_events)){
  COLUMN_event = all_events[j]
  mCA_event = gsub('mCA.', '', COLUMN_event)
  
  mm2 = mm2 %>%
    mutate(new = ifelse(grepl(mCA_event, mCAs), 1, 0))
  colnames(mm2)[ncol(mm2)] = COLUMN_event
  if (COLUMN_event %in% paste0('mCA.chr', seq(1, 23, by=1), '_')){
    colnames(mm2)[ncol(mm2)] = gsub('\\_', '', COLUMN_event)
  }
}

mm2$mCA.total = mm2 %>%
  select(all_of(c(
    paste0('mCA.chr', seq(1, 23, by=1), '_amp'),
    paste0('mCA.chr', seq(1, 23, by=1), '_del'),
    paste0('mCA.chr', seq(1, 23, by=1), '_loh')))) %>%
  apply(1, sum)
mm2$mCA.total.chroms = mm2 %>%
  select(all_of(c(
    paste0('mCA.chr', seq(1, 23, by=1))))) %>%
  apply(1, sum)
mm2$mCA.total.amp = mm2 %>%
  select(all_of(c(
    paste0('mCA.chr', seq(1, 23, by=1), '_amp')))) %>%
  apply(1, sum)
mm2$mCA.total.del = mm2 %>%
  select(all_of(c(
    paste0('mCA.chr', seq(1, 23, by=1), '_del')))) %>%
  apply(1, sum)
mm2$mCA.total.loh = mm2 %>%
  select(all_of(c(
    paste0('mCA.chr', seq(1, 23, by=1), '_loh')))) %>%
  apply(1, sum)

#### Annotate patients with mCA affecting TOP20 ch.genes
ch.genes.cols <- colnames(data_7)[grep('ch.DNMT3A', colnames(data_7))[1]:grep('ch.HIST1H3E', colnames(data_7))[1]]

top_ch_genes = names(data_7 %>% select(ch.genes.cols) %>% apply(2,sum) %>% sort(decreasing=T))[1:30]
top_ch_genes = gsub('ch.', '', top_ch_genes)

dd_mat = mm2
for (i in 1:length(top_ch_genes)){
  gene = top_ch_genes[i]
  term_gene = paste0('mCA.gene.', gene)
  print(gene)
  
  # patients_with_mca_in_gene = mCA_calls %>%
  #   filter(grepl(paste0(';', gene, ';'), fixed=TRUE, genes_affected) &
  #    cnv.conflict==FALSE) %>% select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character()
  patients_with_mca_in_gene = mCA_calls %>%
    filter(grepl(paste0(';', gene, ';'), fixed=TRUE, genes_affected)) %>% select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character()
  
  dd_mat = dd_mat %>%
    mutate(new = ifelse( DMP_PATIENT_ID %in% patients_with_mca_in_gene, 1, 0))
  colnames(dd_mat)[ncol(dd_mat)] = term_gene
}

mm3 = dd_mat
dim(mm3)

## Annotationg mCA finally:
y <- mm3 %>% select(-DMP_PATIENT_ID, -PartC.Consent.Status)

data_7.2 <- data_7 %>%
  left_join(y, by='DMP_ASSAY_ID_Normal')
dim(data_7.2)



################# Annotate functional_ch OR mCA -affected genes ################# 
data_7.3 = data_7.2
ch.genes.cols <- colnames(data_7.3)[grep('ch.DNMT3A', colnames(data_7.3))[1]:grep('ch.HIST1H3E', colnames(data_7.3))[1]]
mca_genes = colnames(data_7.3)[grep('mCA.gene.', colnames(data_7.3), fixed=T)]

v_genes = unique( c( gsub('ch.', '', ch.genes.cols), gsub('mCA.gene.', '', mca_genes)) )

for (i in 1:length(v_genes)){
  GENE = v_genes[i]
    print(GENE)
  CH_col = paste0('ch.', GENE)
  mCA_col = paste0('mCA.gene.', GENE)
  
  if ( mCA_col %in% colnames(data_7.3)){
    data_7.3 = data_7.3 %>%
      mutate(new = ifelse((get(CH_col) == 1 | get(mCA_col) == 1) & !is.na(ch_cnv), 1, ifelse(is.na(ch_cnv), NA, 0))) %>%
      mutate(new_2 = ifelse((get(CH_col) == 1 & get(mCA_col) == 1) & !is.na(ch_cnv), 1, ifelse(is.na(ch_cnv), NA, 0)))
    colnames(data_7.3)[(ncol(data_7.3)-1):ncol(data_7.3)] = c( paste0('ch.OR.mCA.', GENE), paste0('ch.AND.mCA.', GENE))
  } else {
    data_7.3 = data_7.3 %>%
      mutate(new = ifelse(get(CH_col) == 1 & !is.na(ch_cnv), 1, ifelse(is.na(ch_cnv), NA, 0))) %>%
      mutate(new_2 = ifelse(get(CH_col) == 1 & !is.na(ch_cnv), 1, ifelse(is.na(ch_cnv), NA, 0)))
    colnames(data_7.3)[(ncol(data_7.3)-1):ncol(data_7.3)] = c( paste0('ch.OR.mCA.', GENE), paste0('ch.AND.mCA.', GENE))
  }
  

}
data_7.4 = data_7.3
dim(data_7.4)


# ########### FILTER for patients that were run through CH-FACETS:
# data_7.5 = data_7.4 %>% filter(!is.na(ch_cnv))
# dim(data_7.5)

data_7.5 = data_7.4
dim(data_7.5)


##################################################################
###################### GERMLINE VARIANTS ######################
table(germ_cnv_df$cnv.conflict)
germ_cnv_df_2 = germ_cnv_df %>% filter(cnv.conflict==FALSE)

## Data must be reduced to those with germline profile done (SNVs and CNVs)!
patients_germline_profiled <- unique(c(
  data_7.5 %>% filter(PartC.Consent.Status==1) %>% select(DMP_PATIENT_ID) %>%
    unique %>% pull %>% as.character ,
  germ_patients_d %>% select(DMP_PATIENT_ID) %>% 
    unique %>% pull %>% as.character ,
  germ_cnv_df_2 %>% select(DMP_PATIENT_ID) %>% 
    unique %>% pull %>% as.character ))

data_8 <- data_7.5 %>% filter(DMP_PATIENT_ID %in% patients_germline_profiled)
dim(data_8)

## Substetting patients with germline variants
# germ_calls2 <- germ_calls %>% filter(FinalPathScore != '-')
germ_calls2 <- germ_calls


## Annotating Germ STATUS:
# Avoiding duplicated CH variants from multiple time-points in patients

data_9 <- data_8 %>%
  mutate(germ_mutated = ifelse(DMP_PATIENT_ID %in% unique(as.character(germ_calls2$DMP_PATIENT_ID)), 1, 0)) %>%
  mutate(germ_cnv = ifelse(DMP_PATIENT_ID %in% unique(as.character(germ_cnv_df_2$DMP_PATIENT_ID)), 1, 0)) %>%
  mutate(germ_event = ifelse(
    germ_mutated == 1 | germ_cnv == 1, 1, 0)) 
dim(data_9)

## Adding number of mutations:
n_germ_muts__df = germ_calls2 %>%
  count(DMP_PATIENT_ID) %>% rename(n_germ_muts = n)
n_germ_cnvs__df = germ_cnv_df_2 %>%
  count(DMP_PATIENT_ID) %>% rename(n_germ_cnvs = n)

data_10 <- data_9 %>%
  left_join(n_germ_muts__df, by='DMP_PATIENT_ID') %>%
  mutate(n_germ_muts = ifelse(is.na(n_germ_muts), 0, n_germ_muts)) %>%
  left_join(n_germ_cnvs__df, by='DMP_PATIENT_ID') %>%
  mutate(n_germ_cnvs = ifelse(is.na(n_germ_cnvs), 0, n_germ_cnvs)) %>%
  mutate(n_germ_events = n_germ_muts+n_germ_cnvs) %>%
  mutate(germ_mutANDcnv=fifelse(
    germ_mutated==1 & germ_cnv==1, 1, 0))
dim(data_10)

# Checking...
#
# data_10 %>% select(germ_event, germ_mutated, germ_cnv, germ_mutANDcnv, n_germ_events, n_germ_muts, n_germ_cnvs) %>% apply(2, table, useNA='always')


## Annotating binary tags for germ-event genes:
germ_events_genes <- sort(unique(c(
  germ_calls2 %>% select(Hugo_Symbol) %>% unique %>% pull %>% as.character,
  germ_cnv_df_2 %>% select(Gene) %>% unique %>% pull %>% as.character
)))

for (i in 1:length(germ_events_genes)) {
  gene <- germ_events_genes[i]
  print(gene)
  
  ## Gene mutation:
  patients_with_mut.gene <- germ_calls2 %>%
    filter(Hugo_Symbol==gene) %>% select(DMP_PATIENT_ID) %>%
    unique() %>% pull() %>% as.character()
  data_10 <- data_10 %>% mutate(gene=fifelse(DMP_PATIENT_ID%in%patients_with_mut.gene, 1, 0))
  colnames(data_10)[ncol(data_10)] <- paste0('g.', gene)
  
  ## Gene CNV:
  patients_with_cnv_gene <- germ_cnv_df_2 %>%
    filter(Gene==gene) %>% select(DMP_PATIENT_ID) %>%
    unique() %>% pull() %>% as.character()
  data_10 <- data_10 %>% mutate(gene=fifelse(DMP_PATIENT_ID%in%patients_with_cnv_gene, 1, 0))
  colnames(data_10)[ncol(data_10)] <- paste0('cnv.', gene)
  
  ## Both:
  data_10$neither = data_10[,(ncol(data_10)-1)]+data_10[,ncol(data_10)]
  colnames(data_10)[ncol(data_10)] <- paste0('gMutORCNV.', gene)
}
data_11 <- data_10
dim(data_11)

## g.CHEK2 VUS mutation 'p.I157T':
## g.CHEK2_p.I157T:
patients_with_mut.gene <- germ_calls2 %>%
  filter(Hugo_Symbol=='CHEK2' & HGVSp_Short=='p.I157T') %>% select(DMP_PATIENT_ID) %>%
  unique() %>% pull() %>% as.character()
data_11 <- data_11 %>%
  mutate(g.CHEK2_p.I157T=fifelse(DMP_PATIENT_ID%in%patients_with_mut.gene, 1, 0))

## OTHER THAN g.CHEK2_p.I157T:
data_11 <- data_11 %>%
  mutate(gMutORCNV.CHEK2_not_p.I157T = ifelse(
    gMutORCNV.CHEK2==1 & g.CHEK2_p.I157T==0, 1, 0))

# Checking...
#
# data_11 %>% select(g.CHEK2, cnv.CHEK2, gMutORCNV.CHEK2, g.CHEK2_p.I157T, gMutORCNV.CHEK2_not_p.I157T) %>% apply(2, table, useNA='always')









### Annotating two.Hit in germ+CH (germ SNV or CNV + CH variant or mCA)
df_ano <- data_11
df_germ.vars <- germ_calls2
df_ch.calls <- calls_1

##
germMut.genes.cols <- colnames(df_ano)[grep('g.', fixed=TRUE, colnames(df_ano))]
germCNV.genes.cols <- colnames(df_ano)[grep('cnv.', fixed=TRUE, colnames(df_ano))]

germ.genes.cols = sort(unique(c(gsub('g.', '', germMut.genes.cols), gsub('cnv.', '', germCNV.genes.cols))))

ch.genes.cols <- colnames(df_ano)[grep('ch.DNMT3A', colnames(df_ano))[1]:grep('ch.RHOA', colnames(df_ano))[1]]

mCA.cols <- paste0('mCA.chr',1:23)
##

potential_twoHits_genes <- sort(germ.genes.cols[germ.genes.cols %in% gsub('ch.', '',ch.genes.cols)])

df_gene_chrom <- rbind( as.data.frame(potential_twoHits_genes) %>%
                          rename(Gene = potential_twoHits_genes) %>% left_join(
                            df_germ.vars %>%
                              mutate(Gene = Hugo_Symbol) %>%
                              mutate(Chrom = paste0('chr', Chromosome)) %>%
                              # mutate(Chrom = Chromosome) %>%
                              select(Gene, Chrom), by='Gene') %>% unique() %>%
                          filter(!is.na(Chrom) & Chrom != 'chr-'),
                        
                        as.data.frame(potential_twoHits_genes) %>%
                          rename(Gene = potential_twoHits_genes) %>% left_join(
                            df_ch.calls %>%
                              mutate(Gene = Hugo_Symbol) %>%
                              mutate(Chrom = paste0('chr', Chromosome)) %>%
                              # mutate(Chrom = Chromosome) %>%
                              select(Gene, Chrom), by='Gene') %>% unique() %>%
                          filter(!is.na(Chrom) & Chrom != 'chr-')
                        
) %>% unique() %>% as.data.frame()
rownames(df_gene_chrom) <- 1:nrow(df_gene_chrom)

df_gene_chrom$Gene <- as.character(df_gene_chrom$Gene)
df_gene_chrom$Chrom <- as.character(df_gene_chrom$Chrom)


for (i in 1:nrow(df_gene_chrom) ){
  GENE <- as.character(df_gene_chrom$Gene[i])
  # GENE <- 'ALK'
  
  germ.gene <- paste0('g.', GENE)
  germ.CNV.gene <- paste0('cnv.', GENE)
  ch.gene <- paste0('ch.', GENE)
  print(paste('Germ gene is', GENE))
  print(paste('CH gene is', ch.gene))
  
  # potential_mCA_cols <- paste0('mCA.', df_gene_chrom[which(df_gene_chrom$Gene == GENE),'Chrom'], c('_del', '_amp', '_loh'))
  potential_mCA_cols <- paste0('mCA.', df_gene_chrom[which(df_gene_chrom$Gene == GENE),'Chrom'])
  somatic_cols_gene <- paste0('ch.', GENE)
  # somatic_cols_gene <- paste0('ch.',c(GENE, 
  # mCA.cols[mCA.cols %in% potential_mCA_cols]))
  
  germline_cols_gene = c(germ.gene, germ.CNV.gene)
  somatic_cols_gene <- c(somatic_cols_gene,
                         mCA.cols[mCA.cols %in% potential_mCA_cols])
  
  cols_interest <- c(germline_cols_gene, somatic_cols_gene)
  
  dd <- df_ano %>% select(all_of(cols_interest))
  dd <- dd %>%
    mutate(all = dd %>% apply(1,sum, na.rm=T)) %>%
    mutate(germline_profile = dd %>% select(all_of(germline_cols_gene)) %>% apply(1,sum, na.rm=T)) %>%
    mutate(ch_profile = dd %>% select(somatic_cols_gene) %>% apply(1,sum, na.rm=T)) %>%
    mutate(g.OR.ch = fifelse(all>0 & !is.na(all), 1, 0)) %>%
    mutate(g.AND.ch = fifelse((get(germ.gene)==1 & get(ch.gene)==1), 1, 0)) %>%
    mutate(two.Hits = fifelse((germline_profile>0 & ch_profile>0), 1, 0)) 
  
  dd[which(is.na(dd$two.Hits)), 'two.Hits'] <- 0
  
  
  df_ano = cbind(df_ano, dd %>% select(g.OR.ch, g.AND.ch, two.Hits))
  colnames(df_ano)[ncol(df_ano)-2] <- paste0('g.OR.ch.', GENE)
  colnames(df_ano)[ncol(df_ano)-1] <- paste0('g.AND.ch.', GENE)
  colnames(df_ano)[ncol(df_ano)] <- paste0('two.Hits.', GENE)
}
two.Hits_genes <- colnames(df_ano)[grep('two.Hits.', colnames(df_ano))]

data_15 = df_ano %>%
  mutate(two.hits = df_ano %>% select(all_of(two.Hits_genes)) %>%
           apply(1,sum, na.rm=T)) %>%
  mutate(two.hits_bin = fifelse(two.hits > 0 , 1, 0)) %>%
  select(-two.hits) %>%
  rename(two.hits=two.hits_bin)
dim(data_15)





## Selecting major mainTypes
data_ano <- data_15

top_types = 20 # Top cancer types

selected_mainTypes_in_data <- data_ano %>% count(mainType) %>% arrange(desc(n)) %>% slice(1:top_types) %>%select(mainType) %>% pull() %>% as.character()

not_selected_mainTypes_in_data <- data_ano %>% count(mainType) %>% filter(! mainType %in% selected_mainTypes_in_data) %>% select(mainType) %>% pull() %>% as.character()



## Annotating 'Major_mainType' of cancer  
data_ano <- data_ano %>% mutate(Major_mainType=fifelse(mainType%in%selected_mainTypes_in_data, mainType, 'Other cancer type'))



## Binary annotation per patient and by 'Major_mainType':
selected_mainTypes_in_data_and_Other = c(selected_mainTypes_in_data, 'Other cancer type') # Major mainTypes vector

for (i in 1:length(selected_mainTypes_in_data_and_Other)){
  x <- selected_mainTypes_in_data_and_Other[i]
  print(x)
  data_ano <- data_ano %>% mutate(term=fifelse(Major_mainType==x, 1, 0))
  colnames(data_ano)[ncol(data_ano)] <- as.character(x)
}


data_16 = data_ano
dim(data_16)







## OTHER CH features:
genes_DDR <- c('PPM1D', 'TP53', 'CHEK2', 'ATM')
genes_splicing <- c('U2AF1', 'SRSF2', 'SF3B1')
genes_epigenetic_modif <- c('DNMT3A', 'ASXL1', 'TET2')
dnmt3a_ch <- c('DNMT3A')


list_of_CH_types = list(genes_DDR, genes_splicing, genes_epigenetic_modif, dnmt3a_ch)
names(list_of_CH_types) = c('DDR_CH', 'Splicing_CH',
                            'EpigModif_CH', 'DNMT3A_CH')


for (i in 1:length(list_of_CH_types)){
  # i = 1
  column_names = paste0('ch.OR.mCA.', list_of_CH_types[[i]])
  ch_types_name = names(list_of_CH_types)[i]
  
  data_16$new = data_16 %>%
    select(all_of(column_names)) %>%
    apply(1, sum)
  data_16 = data_16 %>% mutate(new = ifelse(new>0, 1, 0))
  colnames(data_16)[ncol(data_16)] = ch_types_name
  
}
data_17 = data_16 %>%
  mutate(Other_CH=ifelse(functional_ch==1 & ch.DNMT3A!=1, 1, 0))
data_17 %>% select(all_of(c(names(list_of_CH_types), 'Other_CH', 'functional_ch', 'ch.DNMT3A'))) %>% apply(2, table, useNA='always')


## HighVAF_CH vs LowVAF_CH
high.vaf.ch.pts = impact_ch %>%
  filter(n.VAF >= 0.10) %>%
  select(Tumor_Sample_Barcode) %>%
  unique %>% pull %>% as.character


data_18 = data_17 %>%
  mutate(HighVAF_CH = ifelse(DMP_ASSAY_ID_Normal %in% high.vaf.ch.pts, 1, 0)) %>%
  mutate(functionalCH.OR.mCA = ifelse(
    (functional_ch==1 | ch_cnv==1) & !is.na(ch_cnv), 1, ifelse(
      is.na(functional_ch) | is.na(ch_cnv), NA, 0))) %>%
  mutate(functionalCH.AND.mCA = ifelse(
    (functional_ch==1 & ch_cnv==1) & !is.na(ch_cnv), 1, ifelse(
      is.na(functional_ch) | is.na(ch_cnv), NA, 0)))
#  %>%
# mutate(LowVAF_CH = ifelse(DMP_ASSAY_ID_Normal %in% low.vaf.ch.pts, 1, 0))

colnames(data_18)

## Adding column for numercial value in MajorCancerType:
MajorCancerTypes_cols <- colnames(data_18)[grep('Non Small Cell Lung Cancer', colnames(data_18)):grep('Other cancer type', colnames(data_18))]

MajorCancerTypes_df <- as.data.frame(cbind(MajorCancerTypes_cols, as.numeric((length(MajorCancerTypes_cols)-1):0)))
colnames(MajorCancerTypes_df) <- c('MajorCancerType', 'numerical_value')

data_18$Major.Cancer.Type <- NA
data_18$Major.Cancer.Type_value <- 0
for (i in MajorCancerTypes_cols){
  print(i)
  
  data_18[which(data_18[,i]==1),'Major.Cancer.Type'] = i
  data_18[which(data_18[,i]==1),'Major.Cancer.Type_value'] = as.numeric(as.character(MajorCancerTypes_df[which(as.character(MajorCancerTypes_df$MajorCancerType) == i),'numerical_value']))
}
data_19 = data_18

## Adding functionalCH.OR.mCA for each cancer type:
for (i in 1:length(MajorCancerTypes_cols)){
  CANCER_TYPE = MajorCancerTypes_cols[i]
  
  data_19 = data_19 %>%
    mutate(new = ifelse(get(CANCER_TYPE)==1 & functional_ch==1, 1, 0))
  colnames(data_19)[ncol(data_19)] = paste0('functional_ch__', CANCER_TYPE)
}


## Annotating those patients with full annotation (columns: c('age_bin', 'Sex', 'ancestry_label', 'smoking_bin', 'Major.Cancer.Type', 'binary_therapy'))

data_20 = data_19 %>%
  mutate(fully_annotated = ifelse(
    !is.na(age_bin) &
      !is.na(Sex) &
      !is.na(ancestry_label) &
      !is.na(smoking_bin) &
      !is.na(Major.Cancer.Type) &
      !is.na(binary_therapy),
    1, 0))







####################################################################################
############## OTHER COLUMNS ############## 

df_ano <- data_20

##
patients_in_df_ano = unique(as.character(df_ano$DMP_PATIENT_ID))

germ.genes.cols <- colnames(df_ano)[grep('gMutORCNV.', colnames(df_ano), fixed=T)]
germ.genes.cols = unique(gsub('gMutORCNV.', '', germ.genes.cols))

gMut.cols = paste0('g.', germ.genes.cols)
gCNV.cols = paste0('cnv.', germ.genes.cols)
g_mut.OR.cnv_genes = paste0('gMutORCNV.', germ.genes.cols)
all.germ.cols = c(g_mut.OR.cnv_genes, gMut.cols, gCNV.cols)

ch.genes.cols <- colnames(df_ano)[grep('ch.DNMT3A', colnames(df_ano))[1]:grep('ch.HIST1H3E', colnames(df_ano))[1]]

mca_cols = paste0('mCA.chr', seq(1:23))
mca_genes = colnames(df_ano)[grep('mCA.gene.', colnames(df_ano), fixed=T)]

two.Hits_genes <- colnames(df_ano)[grep('two.Hits.', colnames(df_ano))]

MajorCancerTypes_cols <- colnames(
  df_ano)[grep('Non.Small.Cell.Lung.Cancer', colnames(df_ano))[1]:grep('Other.cancer.type', colnames(df_ano))[1]]
##



## Annotating samples with mutation in genes related to CH-PD:
# genes extracted from --> 'We annotated variants as oncogenic in myeloid disease (CH-myeloid-PD) if they were in a gene hypothesized to drive myeloid/hematologic malignancies (Supplementary Table 5)' [Bolton et al 2020]
ch.pd_genes <- c('NOTCH1', 'ASXL1', 'DNMT3A', 'FLT3', 'NPM1', 'U2AF1', 'NRAS', 'ETV6', 'GNAS', 'BRAF', 'EZH2', 'RUNX1', 'KIT', 'MPL', 'SF3B1', 'CBL', 'PTEN', 'IDH1', 'TP53', 'IDH2', 'PTPN11', 'GATA1', 'GATA2', 'STAG2', 'BCOR', 'WT1', 'JAK2', 'JAK3', 'KRAS', 'TET2', 'SRSF2', 'PPM1D')


v_vector <- c()
for (i in 1:length(ch.pd_genes)){
  gene <- ch.pd_genes[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- g_mut.OR.cnv_genes[grep(gene, g_mut.OR.cnv_genes, fixed=TRUE)]
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}

# all_pts_ch.pd_genes <- as.vector(do.call(cbind, l_list))
all_pts_ch.pd_genes <- unique(v_vector)

df_ano <- df_ano %>% mutate(ch.pd_germ.gene = fifelse(
  DMP_PATIENT_ID %in% all_pts_ch.pd_genes, 1, 0))



## Annotating samples with mutations in genes DDR genes:

ddr_genes <- c('CHEK2', 'ATM', 'TP53')

v_vector <- c()
for (i in 1:length(ddr_genes)){
  gene <- ddr_genes[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- g_mut.OR.cnv_genes[grep(gene, g_mut.OR.cnv_genes, fixed=TRUE)]
    c_patients <- df_ano %>% filter(get(germ.gene[1])==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}

# all_pts_ddr_genes <- as.vector(do.call(cbind, l_list))
all_pts_ddr_genes <- unique(v_vector)

df_ano <- df_ano %>% mutate(DDR_germ.gene = fifelse(
  DMP_PATIENT_ID %in% all_pts_ddr_genes, 1, 0))


## Annotating samples with mutations in genes Splicing genes:

splicing_genes <- c('U2AF1', 'SRSF2', 'SF3B1')

v_vector <- c()
for (i in 1:length(splicing_genes)){
  gene <- splicing_genes[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- g_mut.OR.cnv_genes[grep(gene, g_mut.OR.cnv_genes, fixed=TRUE)]
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}

# all_pts_splicing_genes <- as.vector(do.call(cbind, l_list))
all_pts_splicing_genes <- unique(v_vector)

df_ano <- df_ano %>% mutate(Splicing_germ.gene = fifelse(
  DMP_PATIENT_ID %in% all_pts_splicing_genes, 1, 0))


## Annotating samples with mutations in genes EpigModif genes:

epigmodif_genes <- c('DNMT3A', 'ASXL1', 'TET2')

v_vector <- c()
for (i in 1:length(epigmodif_genes)){
  gene <- epigmodif_genes[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- g_mut.OR.cnv_genes[grep(gene, g_mut.OR.cnv_genes, fixed=TRUE)]
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}

# all_pts_epigmodif_genes <- as.vector(do.call(cbind, l_list))
all_pts_epigmodif_genes <- unique(v_vector)

df_ano <- df_ano %>% mutate(EpigModif_germ.gene = fifelse(
  DMP_PATIENT_ID %in% all_pts_epigmodif_genes, 1, 0))


## Annotating samples with mutatation in genes related to hematological malignancies development:
# genes extracted from --> https://ascopubs.org/doi/10.1200/JCO.2016.70.8644?url_ver=Z39.88-2003.rfr_id=ori:rid:crossref.org.rfr_dat=cr_pub%20%200pubmed
heme_genes_asco <- c('FANCA', 'FANCB', 'FANCC', 'FANCCD2', 'FANCE', 'FANCF', 'FANCG', 'FANCI', 'FANCJ', 'BRIP1', 'BACH1', 'FANCL', 'FANCM', 'FANCN', 'PALB2', 'FANCO', 'RAD51C', 'FANCP', 'SLX4', 'FANCQ', 'ERCC4', 'FANCR', 'RAD51', 'FANCS', 'BRCA1', 'BRCA2', 'FANCT', 'UBE2T', 'FANCU', 'XRCC2', 'FANCV', 'REV7', 'ATM', 'NBN', 'BLM', 'MLH1', 'MSH2', 'PMS2', 'EPCAM', 'RUNX1', 'CEBP', 'GATA2', 'PAX5', 'ETV6', 'GATA1', 'DKC1', 'TERC', 'TERT', 'NOLA3', 'NOP10', 'NOLA2', 'NHP2', 'TINF2', 'WRAP53', 'TCAB1', 'CTC1', 'RTEL1', 'ACD', 'TPP1', 'PARN', 'NAF1', 'STN1', 'TP53', 'PTPN11', 'CBL', 'ELANE', 'ANKRD26', 'WAS', 'HAX1', 'DDX41', 'SMD9L', 'SRP72', 'NF1')


## heme_genes_asco
v_vector <- c()
for (i in 1:length(heme_genes_asco)){
  gene <- heme_genes_asco[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- paste0('gMutORCNV.', gene)
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}
all_pts_heme_genes <- unique(v_vector)
df_ano <- df_ano %>% mutate(heme_genes_asco = fifelse(
  DMP_PATIENT_ID %in% all_pts_heme_genes, 1, 0))





### 
### From Ahmet's mail: Ozge shared germline hem genes with us ('Heme genes.xlsx')
adult_MDS.AML <- c('ANKRD26', 'CEPBA', 'DDX41', 'ETV6', 'GATA2', 'MBD4', 'RTEL1', 'RUNX1', 'SRP72', 'TERC', 'TERT', 'TP53')
bone_marrow_failure <- c('ACD', 'ANKRD26', 'BRCA1', 'BRCA2', 'BRIP1', 'CEPBA', 'CTC1', 'DDX41', 'DKC1', 'ERCC4', 'FANCA', 'FANCB', 'FANCC', 'FANCD2', 'FANCE', 'FANCF', 'FANCG', 'FANCI', 'FANCJ', 'FANCL', 'FANCM', 'FANCN', 'FANCO', 'FANCP', 'FANCQ', 'FANCR', 'FANCT', 'FANCU', 'GATA1', 'GATA2', 'MAD2L2', 'MPL', 'NAF1', 'NHEJ1', 'NHP2', 'NOP10', 'PALB2', 'PARN', 'POT1', 'RAD51', 'RAD51C', 'REV7', 'RTEL1', 'RUNX1', 'SAMD9', 'SAMD9L', 'SLX4', 'SRP72', 'STN1', 'TERC', 'TERT', 'TINF2', 'TP53', 'UBE2T', 'USB1', 'WRAP53', 'XRCC2')
all <- c('ETV6', 'IKZF1', 'PAX5', 'RUNX1', 'TP53')
lymphoma <- c('CASP10', 'FADD', 'FAS', 'FASLG', 'ITK', 'MAGT1 ', 'NPAT', 'PRF1', 'SH2D1A', 'TP53', 'UNC13D')
heme.malignancy <- c('ARID5B', 'ASXL1', 'ATM', 'BLM', 'BRCA1', 'BRCA2', 'BRIP1', 'BTK', 'CASP8', 'CBL', 'CDKN2A', 'CEBPA', 'CHEK2', 'CREBBP', 'CSF3R', 'EP300', 'ERG', 'ETV6', 'FANCA', 'FANCC', 'FANCD2', 'FAS', 'GATA1', 'GATA2', 'GATA3', 'IKZF1', 'IKZF3', 'JAK2', 'JAK3', 'KMT2A', 'MLH1', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'NBN', 'NF1', 'NPM1', 'PALB2', 'PAX5', 'PMS2', 'POT1', 'PTPN11', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAF1', 'RTEL1', 'RUNX1', 'SETD2', 'SETD5', 'SH2B3', 'SP140', 'STAT5B', 'TERT', 'TET2', 'TP53', 'TYK2', 'WHSC1')

## adult_MDS.AML
v_vector <- c()
for (i in 1:length(adult_MDS.AML)){
  gene <- adult_MDS.AML[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- paste0('gMutORCNV.', gene)
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}
all_pts_heme_genes <- unique(v_vector)
df_ano <- df_ano %>% mutate(adult_MDS.AML = fifelse(
  DMP_PATIENT_ID %in% all_pts_heme_genes, 1, 0))

## bone_marrow_failure
v_vector <- c()
for (i in 1:length(bone_marrow_failure)){
  gene <- bone_marrow_failure[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- paste0('gMutORCNV.', gene)
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}
all_pts_heme_genes <- unique(v_vector)
df_ano <- df_ano %>% mutate(bone_marrow_failure = fifelse(
  DMP_PATIENT_ID %in% all_pts_heme_genes, 1, 0))

## all
v_vector <- c()
for (i in 1:length(all)){
  gene <- all[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- paste0('gMutORCNV.', gene)
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}
all_pts_heme_genes <- unique(v_vector)
df_ano <- df_ano %>% mutate(all = fifelse(
  DMP_PATIENT_ID %in% all_pts_heme_genes, 1, 0))

## lymphoma
v_vector <- c()
for (i in 1:length(lymphoma)){
  gene <- lymphoma[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- paste0('gMutORCNV.', gene)
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}
all_pts_heme_genes <- unique(v_vector)
df_ano <- df_ano %>% mutate(lymphoma = fifelse(
  DMP_PATIENT_ID %in% all_pts_heme_genes, 1, 0))

## heme.malignancy
v_vector <- c()
for (i in 1:length(heme.malignancy)){
  gene <- heme.malignancy[i]
  print(gene)
  if (length(grep(gene, g_mut.OR.cnv_genes, fixed=TRUE))>0){
    print('TRUE')
    germ.gene <- paste0('gMutORCNV.', gene)
    c_patients <- df_ano %>% filter(get(germ.gene)==1) %>%
      select(DMP_PATIENT_ID) %>% unique %>% pull %>% as.character
    v_vector<- c(v_vector, c_patients)
  }
}
all_pts_heme_genes <- unique(v_vector)
df_ano <- df_ano %>% mutate(heme.malignancy = fifelse(
  DMP_PATIENT_ID %in% all_pts_heme_genes, 1, 0))



data_21 = df_ano










#############################################################
###################### EXPORTING ############################
d_asso = data_21
colnames(d_asso)

## Exporting 'd_asso':
NAME = paste0(out_dir,
              'CH.GERM_association.matrix___', Sys.Date(), '.tsv')
write.table(d_asso, NAME, col.names=T, row.names=F, quote=F, sep='\t')

NAME




### Export germline columns:
germ.genes.cols <- colnames(d_asso)[grep('gMutORCNV.', colnames(d_asso), fixed=T)]
germ.genes.cols = unique(gsub('gMutORCNV.', '', germ.genes.cols))
germ.genes.cols = germ.genes.cols[-which(germ.genes.cols=='CHEK2_not_p.I157T')]

gMut.cols = paste0('g.', germ.genes.cols)
gCNV.cols = paste0('cnv.', germ.genes.cols)
g_mut.OR.cnv_genes = paste0('gMutORCNV.', germ.genes.cols)
all.germ.cols = c(g_mut.OR.cnv_genes, gMut.cols, gCNV.cols,
                  'g.CHEK2_p.I157T', 'gMutORCNV.CHEK2_not_p.I157T')

write.table(as.data.frame(all.germ.cols), paste0(out_dir, 'all_germ_cols.txt'), col.names=F, row.names=F, quote=T, sep='\t')


