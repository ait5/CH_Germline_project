# CH_Germline_project

## Files tree
    .
    ├── CH_Germline_project_local.Rproj
    ├── Results
    ├── code
    ├── data
    ├── data_mined
    └── manuscript_figures

## Workflow

1. Run `preparing_matrix_for_associations.R` to:  
    - Create main data matrix: `data_mined/CH.GERM_association.matrix___2022-09-08.tsv`  
    - Annotate copy number variants between g.CNVs and ch.mCAs:  
        - `data_mined/germ_cnvs_and_mCA.conflict_2022-09-08.txt`
        - `data_mined/mCA_calls_and_germ.conflict_2022-09-08.txt`  
  
  
2. Run NextFlow pipelines to compute associations (in Juno):  
    2.1. `code/inJuno/NF_main.associations_inJuno/pipeline.nf` will run `script_main.associations.R` to compute main associations (e.g. functional_ch ~ germ_event); results will set in `Results/Main_associations/tables_main/`  
    2.2. `code/inJuno/NF_main.associations.by.gene_inJuno/pipeline.nf` will run `script_main.associations_by.gene.R` to compute main associations by gene (e.g. functional_ch ~ ATM); results will set in `Results/Main_associations_by.gene/tables/`  
    2.3. `code/inJuno/NF_main.associations.by.gene_inJuno/pipeline.nf` will run `script_main.associations_by.gene.R` to compute main associations by gene (e.g. functional_ch ~ ATM); results will set in `Results/Main_associations_by.gene/tables/`  
  
  
3. Run scripts to merge associations:  
    3.1. `merge_main.associations.R` will merge associations from `Results/Main_associations/tables_main/` into `MainAssociations_2022-09-08.tsv`  
    3.2. `merge_main.associations_by.gene.R` will merge associations from `Results/Main_associations_by.gene/tables/` into `MainAssociations_by.gene_2022-09-08.tsv`  
  
  
4. Run `manuscript_figures/all_figures.R` to plot figures for manuscript in `manuscript_figures/plots_paper/`  


## Data in repository

- `Cohort_50k_CNVs_90Genes_deletions_only_pathogenic_only_v1.txt`:
List of CNV deletions in the 90 TSG studied (output from Miika's germline pipeline)

- `DUPLI_filtered_Cohort_50k_Pathogenic_90Genes_v2.txt`:
List of rare and Pathogenic or Likely-Pathogenic SNVs and Indels calls in the 90 TSG studied (output from Miika's germline pipeline) and filtered from misscalled CH variants (essentially, looking for duplicated calls between germline calls and CH calls ['code/running_duplicates_script.R' from Sebas])

- `germ_90_genes.txt`:
List of TSGs included in the germline assessment for rare and P/LP variants

- `FILTERED_50k.CH.calls_2022-02-16_anonymized.tsv_2022-08-23.txt`:
CH calls from 50K study removing duplicates from multiple sample by patients 

- `List_of_anonymized_Samples.txt`:
List of 50k patients included in the germline calling pipeline

- `ch.50k_clinical_annots_by_SAMPLE_2021-12-06_anonymized.tsv`:
Clincal data annotation for patients in the 50K CH characterization (i.e. dmps id, age, inferred ancestry, smoking history, therapy status, Teng's mCA calls, clinical follow-up from Gao et al,  etc)

- `calls_mCA.50k.annotations_2022-03-07_anonymized.tsv`:
List of copy number alterations calls from running FACETS-CH (from Gao et.al 2020)

- `mCA.50k.annotations_2022-03-07_anonymized.tsv`:
Data matrix with annotations for mCA calls

- `oncotree_conversion_codes.txt`:
OncoTree codes and conversion columns to different levels in OncoTree classification

- `refGene_hg19.txt`:
Genomic regions for genes in hg19 genomic assembly (from UCSC Genome Browser database)



## Tree (level 2)
    .
    ├── CH_Germline_project_local.Rproj
    ├── Results
    │   ├── Main_associations
    │   ├── Main_associations_by.gene
    │   └── tables_gene_by_gene
    ├── code
    │   ├── inJuno
    │   ├── merge_main.associations.R
    │   ├── merge_main.associations_by.gene.R
    │   └── preparing_matrix_for_associations.R
    ├── data
    │   ├── Cohort_50k_CNVs_90Genes_deletions_only_pathogenic_only_v1.txt
    │   ├── DUPLI_filtered_Cohort_50k_Pathogenic_90Genes_v2.txt
    │   ├── FILTERED_50k.CH.calls_2022-02-16_anonymized.tsv_2022-08-23.txt
    │   ├── List_of_anonymized_Samples.txt
    │   ├── calls_mCA.50k.annotations_2022-03-07_anonymized.tsv
    │   ├── ch.50k_clinical_annots_by_SAMPLE_2021-12-06_anonymized.tsv
    │   ├── germ_90_genes.txt
    │   ├── mCA.50k.annotations_2022-03-07_anonymized.tsv
    │   ├── oncotree_conversion_codes.txt
    │   └── refGene_hg19.txt
    ├── data_mined
    │   ├── CH.GERM_association.matrix___2022-09-08.tsv
    │   ├── all_germ_cols.txt
    │   ├── germ_cnvs_and_mCA.conflict_2022-09-08.txt
    │   └── mCA_calls_and_germ.conflict_2022-09-08.txt
    └── manuscript_figures
        ├── all_figures.R
        ├── figure_scripts
        ├── funs
        ├── plots_paper
        └── script_all_figures.R


## Full tree:
    ├── CH_Germline_project_local.Rproj
    ├── Results
    │   ├── Main_associations
    │   │   ├── MainAssociations_2022-09-08.tsv
    │   │   └── tables_main
    │   │       ├── MainAssociations_DDR_germ.gene_2022-09-08.tsv
    │   │       ├── MainAssociations_adult_MDS.AML_2022-09-08.tsv
    │   │       ├── MainAssociations_all_2022-09-08.tsv
    │   │       ├── MainAssociations_bone_marrow_failure_2022-09-08.tsv
    │   │       ├── MainAssociations_ch.pd_germ.gene_2022-09-08.tsv
    │   │       ├── MainAssociations_germ_cnv_2022-09-08.tsv
    │   │       ├── MainAssociations_germ_event_2022-09-08.tsv
    │   │       ├── MainAssociations_germ_mutated_2022-09-08.tsv
    │   │       ├── MainAssociations_heme.malignancy_2022-09-08.tsv
    │   │       ├── MainAssociations_heme_genes_asco_2022-09-08.tsv
    │   │       └── MainAssociations_lymphoma_2022-09-08.tsv
    │   ├── Main_associations_by.gene
    │   │   ├── MainAssociations_by.gene_2022-09-08.tsv
    │   │   └── tables
    │   │       ├── MainAssociations_by._cnv.ATM_.tsv
    │   │       ├── MainAssociations_by._cnv.BAP1_.tsv
    │   │       ├── MainAssociations_by._cnv.BARD1_.tsv
    │   │       ├── MainAssociations_by._cnv.BMPR1A_.tsv
    │   │       ├── MainAssociations_by._cnv.BRCA1_.tsv
    │   │       ├── MainAssociations_by._cnv.BRCA2_.tsv
    │   │       ├── MainAssociations_by._cnv.BRIP1_.tsv
    │   │       ├── MainAssociations_by._cnv.CDC73_.tsv
    │   │       ├── MainAssociations_by._cnv.CDH1_.tsv
    │   │       ├── MainAssociations_by._cnv.CHEK2_.tsv
    │   │       ├── MainAssociations_by._cnv.EPCAM_.tsv
    │   │       ├── MainAssociations_by._cnv.FANCA_.tsv
    │   │       ├── MainAssociations_by._cnv.FANCC_.tsv
    │   │       ├── MainAssociations_by._cnv.FH_.tsv
    │   │       ├── MainAssociations_by._cnv.FLCN_.tsv
    │   │       ├── MainAssociations_by._cnv.HOXB13_.tsv
    │   │       ├── MainAssociations_by._cnv.MEN1_.tsv
    │   │       ├── MainAssociations_by._cnv.MLH1_.tsv
    │   │       ├── MainAssociations_by._cnv.MSH2_.tsv
    │   │       ├── MainAssociations_by._cnv.MSH3_.tsv
    │   │       ├── MainAssociations_by._cnv.NBN_.tsv
    │   │       ├── MainAssociations_by._cnv.NF1_.tsv
    │   │       ├── MainAssociations_by._cnv.PALB2_.tsv
    │   │       ├── MainAssociations_by._cnv.PHOX2B_.tsv
    │   │       ├── MainAssociations_by._cnv.PMS2_.tsv
    │   │       ├── MainAssociations_by._cnv.POLD1_.tsv
    │   │       ├── MainAssociations_by._cnv.POLE_.tsv
    │   │       ├── MainAssociations_by._cnv.RAD51B_.tsv
    │   │       ├── MainAssociations_by._cnv.RAD51C_.tsv
    │   │       ├── MainAssociations_by._cnv.RAD51D_.tsv
    │   │       ├── MainAssociations_by._cnv.RAD51_.tsv
    │   │       ├── MainAssociations_by._cnv.RB1_.tsv
    │   │       ├── MainAssociations_by._cnv.SDHA_.tsv
    │   │       ├── MainAssociations_by._cnv.SDHC_.tsv
    │   │       ├── MainAssociations_by._cnv.SDHD_.tsv
    │   │       ├── MainAssociations_by._cnv.SMARCA4_.tsv
    │   │       ├── MainAssociations_by._cnv.SMARCB1_.tsv
    │   │       ├── MainAssociations_by._cnv.TGFBR1_.tsv
    │   │       ├── MainAssociations_by._cnv.TMEM127_.tsv
    │   │       ├── MainAssociations_by._g.APC_.tsv
    │   │       ├── MainAssociations_by._g.ATM_.tsv
    │   │       ├── MainAssociations_by._g.BAP1_.tsv
    │   │       ├── MainAssociations_by._g.BARD1_.tsv
    │   │       ├── MainAssociations_by._g.BLM_.tsv
    │   │       ├── MainAssociations_by._g.BMPR1A_.tsv
    │   │       ├── MainAssociations_by._g.BRCA1_.tsv
    │   │       ├── MainAssociations_by._g.BRCA2_.tsv
    │   │       ├── MainAssociations_by._g.BRIP1_.tsv
    │   │       ├── MainAssociations_by._g.CDH1_.tsv
    │   │       ├── MainAssociations_by._g.CDKN2A_.tsv
    │   │       ├── MainAssociations_by._g.CHEK2_.tsv
    │   │       ├── MainAssociations_by._g.CHEK2_p.I157T_.tsv
    │   │       ├── MainAssociations_by._g.DICER1_.tsv
    │   │       ├── MainAssociations_by._g.EGFR_.tsv
    │   │       ├── MainAssociations_by._g.ERBB2_.tsv
    │   │       ├── MainAssociations_by._g.ERCC3_.tsv
    │   │       ├── MainAssociations_by._g.FANCA_.tsv
    │   │       ├── MainAssociations_by._g.FANCC_.tsv
    │   │       ├── MainAssociations_by._g.FH_.tsv
    │   │       ├── MainAssociations_by._g.FLCN_.tsv
    │   │       ├── MainAssociations_by._g.HOXB13_.tsv
    │   │       ├── MainAssociations_by._g.LZTR1_.tsv
    │   │       ├── MainAssociations_by._g.MEN1_.tsv
    │   │       ├── MainAssociations_by._g.MITF_.tsv
    │   │       ├── MainAssociations_by._g.MLH1_.tsv
    │   │       ├── MainAssociations_by._g.MSH2_.tsv
    │   │       ├── MainAssociations_by._g.MSH3_.tsv
    │   │       ├── MainAssociations_by._g.MSH6_.tsv
    │   │       ├── MainAssociations_by._g.MUTYH_.tsv
    │   │       ├── MainAssociations_by._g.NBN_.tsv
    │   │       ├── MainAssociations_by._g.NF1_.tsv
    │   │       ├── MainAssociations_by._g.NTHL1_.tsv
    │   │       ├── MainAssociations_by._g.PALB2_.tsv
    │   │       ├── MainAssociations_by._g.PHOX2B_.tsv
    │   │       ├── MainAssociations_by._g.PMS2_.tsv
    │   │       ├── MainAssociations_by._g.PTCH1_.tsv
    │   │       ├── MainAssociations_by._g.PTEN_.tsv
    │   │       ├── MainAssociations_by._g.RAD51B_.tsv
    │   │       ├── MainAssociations_by._g.RAD51C_.tsv
    │   │       ├── MainAssociations_by._g.RAD51D_.tsv
    │   │       ├── MainAssociations_by._g.RAD51_.tsv
    │   │       ├── MainAssociations_by._g.RB1_.tsv
    │   │       ├── MainAssociations_by._g.RECQL_.tsv
    │   │       ├── MainAssociations_by._g.RET_.tsv
    │   │       ├── MainAssociations_by._g.RTEL1_.tsv
    │   │       ├── MainAssociations_by._g.SDHAF2_.tsv
    │   │       ├── MainAssociations_by._g.SDHA_.tsv
    │   │       ├── MainAssociations_by._g.SDHB_.tsv
    │   │       ├── MainAssociations_by._g.SDHC_.tsv
    │   │       ├── MainAssociations_by._g.SDHD_.tsv
    │   │       ├── MainAssociations_by._g.SMAD3_.tsv
    │   │       ├── MainAssociations_by._g.SMAD4_.tsv
    │   │       ├── MainAssociations_by._g.SMARCA4_.tsv
    │   │       ├── MainAssociations_by._g.SMARCB1_.tsv
    │   │       ├── MainAssociations_by._g.STK11_.tsv
    │   │       ├── MainAssociations_by._g.SUFU_.tsv
    │   │       ├── MainAssociations_by._g.TERT_.tsv
    │   │       ├── MainAssociations_by._g.TGFBR1_.tsv
    │   │       ├── MainAssociations_by._g.TMEM127_.tsv
    │   │       ├── MainAssociations_by._g.TP53_.tsv
    │   │       ├── MainAssociations_by._g.TRIP13_.tsv
    │   │       ├── MainAssociations_by._g.TSC2_.tsv
    │   │       ├── MainAssociations_by._g.VHL_.tsv
    │   │       ├── MainAssociations_by._g.YAP1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.APC_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.ATM_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BAP1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BARD1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BLM_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BMPR1A_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BRCA1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BRCA2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.BRIP1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.CDC73_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.CDH1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.CDKN2A_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.CHEK2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.CHEK2_not_p.I157T_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.DICER1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.EGFR_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.EPCAM_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.ERBB2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.ERCC3_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.FANCA_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.FANCC_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.FH_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.FLCN_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.HOXB13_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.LZTR1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MEN1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MITF_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MLH1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MSH2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MSH3_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MSH6_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.MUTYH_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.NBN_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.NF1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.NTHL1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.PALB2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.PHOX2B_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.PMS2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.POLD1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.POLE_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.PTCH1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.PTEN_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RAD51B_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RAD51C_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RAD51D_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RAD51_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RB1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RECQL_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RET_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.RTEL1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SDHAF2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SDHA_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SDHB_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SDHC_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SDHD_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SMAD3_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SMAD4_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SMARCA4_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SMARCB1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.STK11_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.SUFU_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TERT_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TGFBR1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TMEM127_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TP53_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TRIP13_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TSC1_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.TSC2_.tsv
    │   │       ├── MainAssociations_by._gMutORCNV.VHL_.tsv
    │   │       └── MainAssociations_by._gMutORCNV.YAP1_.tsv
    │   └── tables_gene_by_gene
    │       ├── cnv.ALK_ass.with_CH.genes.tsv
    │       ├── cnv.APC_ass.with_CH.genes.tsv
    │       ├── cnv.ATM_ass.with_CH.genes.tsv
    │       ├── cnv.BAP1_ass.with_CH.genes.tsv
    │       ├── cnv.BARD1_ass.with_CH.genes.tsv
    │       ├── cnv.BLM_ass.with_CH.genes.tsv
    │       ├── cnv.BMPR1A_ass.with_CH.genes.tsv
    │       ├── cnv.BRCA1_ass.with_CH.genes.tsv
    │       ├── cnv.BRCA2_ass.with_CH.genes.tsv
    │       ├── cnv.BRIP1_ass.with_CH.genes.tsv
    │       ├── cnv.CDC73_ass.with_CH.genes.tsv
    │       ├── cnv.CDH1_ass.with_CH.genes.tsv
    │       ├── cnv.CDK4_ass.with_CH.genes.tsv
    │       ├── cnv.CDKN2A_ass.with_CH.genes.tsv
    │       ├── cnv.CHEK2_ass.with_CH.genes.tsv
    │       ├── cnv.CTR9_ass.with_CH.genes.tsv
    │       ├── cnv.DICER1_ass.with_CH.genes.tsv
    │       ├── cnv.EGFR_ass.with_CH.genes.tsv
    │       ├── cnv.EPCAM_ass.with_CH.genes.tsv
    │       ├── cnv.ERBB2_ass.with_CH.genes.tsv
    │       ├── cnv.ERCC3_ass.with_CH.genes.tsv
    │       ├── cnv.ETV6_ass.with_CH.genes.tsv
    │       ├── cnv.FANCA_ass.with_CH.genes.tsv
    │       ├── cnv.FANCC_ass.with_CH.genes.tsv
    │       ├── cnv.FH_ass.with_CH.genes.tsv
    │       ├── cnv.FLCN_ass.with_CH.genes.tsv
    │       ├── cnv.HOXB13_ass.with_CH.genes.tsv
    │       ├── cnv.KIT_ass.with_CH.genes.tsv
    │       ├── cnv.LZTR1_ass.with_CH.genes.tsv
    │       ├── cnv.MAX_ass.with_CH.genes.tsv
    │       ├── cnv.MEN1_ass.with_CH.genes.tsv
    │       ├── cnv.MET_ass.with_CH.genes.tsv
    │       ├── cnv.MITF_ass.with_CH.genes.tsv
    │       ├── cnv.MLH1_ass.with_CH.genes.tsv
    │       ├── cnv.MSH2_ass.with_CH.genes.tsv
    │       ├── cnv.MSH3_ass.with_CH.genes.tsv
    │       ├── cnv.MSH6_ass.with_CH.genes.tsv
    │       ├── cnv.MUTYH_ass.with_CH.genes.tsv
    │       ├── cnv.NBN_ass.with_CH.genes.tsv
    │       ├── cnv.NF1_ass.with_CH.genes.tsv
    │       ├── cnv.NF2_ass.with_CH.genes.tsv
    │       ├── cnv.NTHL1_ass.with_CH.genes.tsv
    │       ├── cnv.PALB2_ass.with_CH.genes.tsv
    │       ├── cnv.PHOX2B_ass.with_CH.genes.tsv
    │       ├── cnv.PMS2_ass.with_CH.genes.tsv
    │       ├── cnv.POLD1_ass.with_CH.genes.tsv
    │       ├── cnv.POLE_ass.with_CH.genes.tsv
    │       ├── cnv.PTCH1_ass.with_CH.genes.tsv
    │       ├── cnv.PTEN_ass.with_CH.genes.tsv
    │       ├── cnv.RAD51B_ass.with_CH.genes.tsv
    │       ├── cnv.RAD51C_ass.with_CH.genes.tsv
    │       ├── cnv.RAD51D_ass.with_CH.genes.tsv
    │       ├── cnv.RAD51_ass.with_CH.genes.tsv
    │       ├── cnv.RB1_ass.with_CH.genes.tsv
    │       ├── cnv.RECQL_ass.with_CH.genes.tsv
    │       ├── cnv.RET_ass.with_CH.genes.tsv
    │       ├── cnv.RTEL1_ass.with_CH.genes.tsv
    │       ├── cnv.RUNX1_ass.with_CH.genes.tsv
    │       ├── cnv.SDHAF2_ass.with_CH.genes.tsv
    │       ├── cnv.SDHA_ass.with_CH.genes.tsv
    │       ├── cnv.SDHB_ass.with_CH.genes.tsv
    │       ├── cnv.SDHC_ass.with_CH.genes.tsv
    │       ├── cnv.SDHD_ass.with_CH.genes.tsv
    │       ├── cnv.SMAD3_ass.with_CH.genes.tsv
    │       ├── cnv.SMAD4_ass.with_CH.genes.tsv
    │       ├── cnv.SMARCA4_ass.with_CH.genes.tsv
    │       ├── cnv.SMARCB1_ass.with_CH.genes.tsv
    │       ├── cnv.SMARCE1_ass.with_CH.genes.tsv
    │       ├── cnv.STK11_ass.with_CH.genes.tsv
    │       ├── cnv.SUFU_ass.with_CH.genes.tsv
    │       ├── cnv.TERT_ass.with_CH.genes.tsv
    │       ├── cnv.TGFBR1_ass.with_CH.genes.tsv
    │       ├── cnv.TMEM127_ass.with_CH.genes.tsv
    │       ├── cnv.TP53_ass.with_CH.genes.tsv
    │       ├── cnv.TRIP13_ass.with_CH.genes.tsv
    │       ├── cnv.TSC1_ass.with_CH.genes.tsv
    │       ├── cnv.TSC2_ass.with_CH.genes.tsv
    │       ├── cnv.VHL_ass.with_CH.genes.tsv
    │       ├── cnv.WT1_ass.with_CH.genes.tsv
    │       ├── cnv.YAP1_ass.with_CH.genes.tsv
    │       ├── g.ALK_ass.with_CH.genes.tsv
    │       ├── g.APC_ass.with_CH.genes.tsv
    │       ├── g.ATM_ass.with_CH.genes.tsv
    │       ├── g.BAP1_ass.with_CH.genes.tsv
    │       ├── g.BARD1_ass.with_CH.genes.tsv
    │       ├── g.BLM_ass.with_CH.genes.tsv
    │       ├── g.BMPR1A_ass.with_CH.genes.tsv
    │       ├── g.BRCA1_ass.with_CH.genes.tsv
    │       ├── g.BRCA2_ass.with_CH.genes.tsv
    │       ├── g.BRIP1_ass.with_CH.genes.tsv
    │       ├── g.CDC73_ass.with_CH.genes.tsv
    │       ├── g.CDH1_ass.with_CH.genes.tsv
    │       ├── g.CDK4_ass.with_CH.genes.tsv
    │       ├── g.CDKN2A_ass.with_CH.genes.tsv
    │       ├── g.CHEK2_ass.with_CH.genes.tsv
    │       ├── g.CHEK2_p.I157T_ass.with_CH.genes.tsv
    │       ├── g.CTR9_ass.with_CH.genes.tsv
    │       ├── g.DICER1_ass.with_CH.genes.tsv
    │       ├── g.EGFR_ass.with_CH.genes.tsv
    │       ├── g.EPCAM_ass.with_CH.genes.tsv
    │       ├── g.ERBB2_ass.with_CH.genes.tsv
    │       ├── g.ERCC3_ass.with_CH.genes.tsv
    │       ├── g.ETV6_ass.with_CH.genes.tsv
    │       ├── g.FANCA_ass.with_CH.genes.tsv
    │       ├── g.FANCC_ass.with_CH.genes.tsv
    │       ├── g.FH_ass.with_CH.genes.tsv
    │       ├── g.FLCN_ass.with_CH.genes.tsv
    │       ├── g.HOXB13_ass.with_CH.genes.tsv
    │       ├── g.KIT_ass.with_CH.genes.tsv
    │       ├── g.LZTR1_ass.with_CH.genes.tsv
    │       ├── g.MAX_ass.with_CH.genes.tsv
    │       ├── g.MEN1_ass.with_CH.genes.tsv
    │       ├── g.MET_ass.with_CH.genes.tsv
    │       ├── g.MITF_ass.with_CH.genes.tsv
    │       ├── g.MLH1_ass.with_CH.genes.tsv
    │       ├── g.MSH2_ass.with_CH.genes.tsv
    │       ├── g.MSH3_ass.with_CH.genes.tsv
    │       ├── g.MSH6_ass.with_CH.genes.tsv
    │       ├── g.MUTYH_ass.with_CH.genes.tsv
    │       ├── g.NBN_ass.with_CH.genes.tsv
    │       ├── g.NF1_ass.with_CH.genes.tsv
    │       ├── g.NF2_ass.with_CH.genes.tsv
    │       ├── g.NTHL1_ass.with_CH.genes.tsv
    │       ├── g.PALB2_ass.with_CH.genes.tsv
    │       ├── g.PHOX2B_ass.with_CH.genes.tsv
    │       ├── g.PMS2_ass.with_CH.genes.tsv
    │       ├── g.POLD1_ass.with_CH.genes.tsv
    │       ├── g.POLE_ass.with_CH.genes.tsv
    │       ├── g.PTCH1_ass.with_CH.genes.tsv
    │       ├── g.PTEN_ass.with_CH.genes.tsv
    │       ├── g.RAD51B_ass.with_CH.genes.tsv
    │       ├── g.RAD51C_ass.with_CH.genes.tsv
    │       ├── g.RAD51D_ass.with_CH.genes.tsv
    │       ├── g.RAD51_ass.with_CH.genes.tsv
    │       ├── g.RB1_ass.with_CH.genes.tsv
    │       ├── g.RECQL_ass.with_CH.genes.tsv
    │       ├── g.RET_ass.with_CH.genes.tsv
    │       ├── g.RTEL1_ass.with_CH.genes.tsv
    │       ├── g.RUNX1_ass.with_CH.genes.tsv
    │       ├── g.SDHAF2_ass.with_CH.genes.tsv
    │       ├── g.SDHA_ass.with_CH.genes.tsv
    │       ├── g.SDHB_ass.with_CH.genes.tsv
    │       ├── g.SDHC_ass.with_CH.genes.tsv
    │       ├── g.SDHD_ass.with_CH.genes.tsv
    │       ├── g.SMAD3_ass.with_CH.genes.tsv
    │       ├── g.SMAD4_ass.with_CH.genes.tsv
    │       ├── g.SMARCA4_ass.with_CH.genes.tsv
    │       ├── g.SMARCB1_ass.with_CH.genes.tsv
    │       ├── g.SMARCE1_ass.with_CH.genes.tsv
    │       ├── g.STK11_ass.with_CH.genes.tsv
    │       ├── g.SUFU_ass.with_CH.genes.tsv
    │       ├── g.TERT_ass.with_CH.genes.tsv
    │       ├── g.TGFBR1_ass.with_CH.genes.tsv
    │       ├── g.TMEM127_ass.with_CH.genes.tsv
    │       ├── g.TP53_ass.with_CH.genes.tsv
    │       ├── g.TRIP13_ass.with_CH.genes.tsv
    │       ├── g.TSC1_ass.with_CH.genes.tsv
    │       ├── g.TSC2_ass.with_CH.genes.tsv
    │       ├── g.VHL_ass.with_CH.genes.tsv
    │       ├── g.WT1_ass.with_CH.genes.tsv
    │       ├── g.YAP1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.ALK_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.APC_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.ATM_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BAP1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BARD1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BLM_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BMPR1A_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BRCA1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BRCA2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.BRIP1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CDC73_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CDH1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CDK4_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CDKN2A_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CHEK2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CHEK2_not_p.I157T_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.CTR9_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.DICER1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.EGFR_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.EPCAM_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.ERBB2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.ERCC3_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.ETV6_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.FANCA_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.FANCC_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.FH_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.FLCN_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.HOXB13_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.KIT_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.LZTR1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MAX_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MEN1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MET_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MITF_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MLH1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MSH2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MSH3_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MSH6_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.MUTYH_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.NBN_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.NF1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.NF2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.NTHL1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.PALB2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.PHOX2B_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.PMS2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.POLD1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.POLE_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.PTCH1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.PTEN_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RAD51B_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RAD51C_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RAD51D_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RAD51_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RB1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RECQL_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RET_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RTEL1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.RUNX1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SDHAF2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SDHA_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SDHB_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SDHC_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SDHD_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SMAD3_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SMAD4_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SMARCA4_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SMARCB1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SMARCE1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.STK11_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.SUFU_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TERT_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TGFBR1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TMEM127_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TP53_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TRIP13_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TSC1_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.TSC2_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.VHL_ass.with_CH.genes.tsv
    │       ├── gMutORCNV.WT1_ass.with_CH.genes.tsv
    │       └── gMutORCNV.YAP1_ass.with_CH.genes.tsv
    ├── code
    │   ├── inJuno
    │   │   ├── NF_main.associations.by.gene_inJuno
    │   │   │   ├── all_germ_cols.txt
    │   │   │   ├── nextflow.config
    │   │   │   ├── pipeline.nf
    │   │   │   ├── script_main.associations_by.gene.R
    │   │   │   ├── v2_pipeline.nf
    │   │   │   └── vNF_script_BY.germ.GENE_by_ch.gene.R
    │   │   └── NF_main.associations_inJuno
    │   │       ├── all_germ_conditions.txt
    │   │       ├── nextflow.config
    │   │       ├── pipeline.nf
    │   │       └── script_main.associations.R
    │   ├── merge_main.associations.R
    │   ├── merge_main.associations_by.gene.R
    │   └── preparing_matrix_for_associations.R
    ├── data
    │   ├── Cohort_50k_CNVs_90Genes_deletions_only_pathogenic_only_v1.txt
    │   ├── DUPLI_filtered_Cohort_50k_Pathogenic_90Genes_v2.txt
    │   ├── FILTERED_50k.CH.calls_2022-02-16_anonymized.tsv_2022-08-23.txt
    │   ├── List_of_anonymized_Samples.txt
    │   ├── calls_mCA.50k.annotations_2022-03-07_anonymized.tsv
    │   ├── ch.50k_clinical_annots_by_SAMPLE_2021-12-06_anonymized.tsv
    │   ├── germ_90_genes.txt
    │   ├── mCA.50k.annotations_2022-03-07_anonymized.tsv
    │   ├── oncotree_conversion_codes.txt
    │   └── refGene_hg19.txt
    ├── data_mined
    │   ├── CH.GERM_association.matrix___2022-09-08.tsv
    │   ├── all_germ_cols.txt
    │   ├── germ_cnvs_and_mCA.conflict_2022-09-08.txt
    │   └── mCA_calls_and_germ.conflict_2022-09-08.txt
    └── manuscript_figures
        ├── all_figures.R
        ├── figure_scripts
        │   ├── figure_1.R
        │   ├── figure_2.R
        │   ├── figure_3.R
        │   ├── supp_figure_1.R
        │   ├── supp_figure_2.R
        │   ├── supp_figure_3.R
        │   ├── supp_figure_4.R
        │   ├── supp_figure_5.R
        │   ├── supp_figure_6.R
        │   ├── supp_figure_7.R
        │   ├── supp_figure_8.R
        │   └── supp_figure_9.R
        ├── funs
        │   ├── fun_CH.PD_stacked.R
        │   ├── fun_CH.rates.by.AgeBin.in.gGene.R
        │   ├── fun_CH.rates.by.CONFOUNDER.R
        │   ├── fun_CH.rates.by.Cancer.Types.R
        │   ├── fun_CH.rates.by.OncotreeCode.R
        │   ├── fun_DDR_mutual_exclusivity.R
        │   ├── fun_barplot.AGE.ch.R
        │   ├── fun_barplot.AGE.ch_BY.GENE.R
        │   ├── fun_barplot.N.muts.R
        │   ├── fun_barplot_mcas.R
        │   ├── fun_boxplot.R
        │   ├── fun_clonal_fraction_by_vaf.R
        │   ├── fun_general.stacked.bp.R
        │   ├── fun_glm_odds_cancers.R
        │   ├── fun_ht.gene.v.gene.R
        │   ├── fun_ht.gene.v.gene_DDR.R
        │   ├── fun_ht.main.associations.by.gene.R
        │   ├── fun_ht.main.associations.by.gene_mCAs.R
        │   ├── fun_ht_mutual.exclusivity.R
        │   ├── fun_mCA.types.barplot.R
        │   ├── fun_main.associations_germ.event.R
        │   ├── fun_main.associations_germ.event_mCAs.R
        │   ├── fun_main.plot.by.gene.R
        │   ├── fun_plot_confounders_by_data.R
        │   ├── fun_plot_odds_confounders.R
        │   ├── fun_summary_table_confounders.R
        │   ├── fun_surv.plot.R
        │   ├── fun_table_to_export.R
        │   ├── fun_top.mutated.genes.R
        │   ├── fun_vaf.ch.plot.R
        │   ├── fun_venn.diagram.R
        │   ├── v2_fun_ht.gene.v.gene_DDR.R
        │   └── v2fun_bp.R
        ├── plots_paper
        │   ├── Fig.1
        │   │   ├── Figure_1_2022-09-09.png
        │   │   ├── General_CH.fraction_2022-09-09.png
        │   │   ├── General_Germ.fraction_2022-09-09.png
        │   │   ├── barplot_by.age.bin_CH.v.Germ_2022-09-09.png
        │   │   ├── p_GERMLINE_Main_associations2022-09-09.png
        │   │   ├── top.genes_germ~CH.assoc_2022-09-09.png
        │   │   ├── top.genes_germ~CH.assoc_2022-09-09.tsv
        │   │   └── venn_CH.v.Germline_2022-09-09.png
        │   ├── Fig.2
        │   │   ├── Figure_2_2022-09-09.png
        │   │   ├── TABLE_top19_ht_Gene_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.tsv
        │   │   ├── ht_mutExcl_DDR_Germ.v.CH__2022-09-09.png
        │   │   └── top19_ht_Gene_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.png
        │   ├── Fig.3
        │   │   ├── Figure_3_2022-09-09.png
        │   │   ├── HR_forest_2022-09-09.png
        │   │   └── sur.plot_GroupsSelected_2022-09-09.png
        │   ├── Supp_Fig.1
        │   │   ├── CH-PDvsCH.fraction_2022-09-09.png
        │   │   ├── CH.muts_per.patient_2022-09-09.png
        │   │   ├── Germ.events_per.patient_2022-09-09.png
        │   │   ├── Supp_Figure_1_2022-09-09.png
        │   │   ├── top10_CH-genes_2022-09-09.png
        │   │   ├── top10_Germ.genes_EVENTS_2022-09-09.png
        │   │   ├── vaf_CH.PD_mutations_2022-09-09.png
        │   │   ├── vaf_CH.mutations_2022-09-09.png
        │   │   └── vaf_Germ.variants_2022-09-09.png
        │   ├── Supp_Fig.2
        │   │   ├── Supp_Figure_2_2022-09-09.png
        │   │   ├── Table_Freq_UniOdds_MultiOdds_2022-09-08.png
        │   │   ├── Table_Freq_UniOdds_MultiOdds_2022-09-09.html
        │   │   ├── Table_Freq_UniOdds_MultiOdds_2022-09-09.png
        │   │   ├── tbl_confounders_Multi.ODDs__2022-09-09.tsv
        │   │   ├── tbl_confounders_Uni.ODDs__2022-09-09.tsv
        │   │   └── tbl_confounders_notCH.VS.CH_freqs__2022-09-09.tsv
        │   ├── Supp_Fig.3
        │   │   ├── Supp_Figure_3_2022-09-09.png
        │   │   ├── TABLE_ht_Gene_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.tsv
        │   │   ├── ht_Gene_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.png
        │   │   └── ht_byGene_Main_associations__(Germ\ events)__glm__2022-09-09.png
        │   ├── Supp_Fig.4
        │   │   ├── All_Supp_Figure_4_2022-09-09.png
        │   │   ├── Supp_Figure_4_2022-09-09.png
        │   │   ├── gMutORCNV.APC_Supp_Figure_4_2022-09-09.png
        │   │   ├── gMutORCNV.APC_table_2022-09-09.tsv
        │   │   ├── gMutORCNV.ATM_Supp_Figure_4_2022-09-09.png
        │   │   ├── gMutORCNV.ATM_table_2022-09-09.tsv
        │   │   ├── gMutORCNV.CHEK2_Supp_Figure_4_2022-09-09.png
        │   │   └── gMutORCNV.CHEK2_table_2022-09-09.tsv
        │   ├── Supp_Fig.5
        │   │   ├── Supp_Figure_5_2022-09-09.png
        │   │   ├── TABLE_ht_Gene.v.Cancer__(SNVs\ +\ CNVs)__glm__2022-09-09.tsv
        │   │   ├── ht_Gene.v.Cancer__(SNVs\ +\ CNVs)__glm__2022-09-09.png
        │   │   └── ht_byGene_CancerTypes__(Germ\ events)__glm__2022-09-09.png
        │   ├── Supp_Fig.6
        │   │   ├── ATM.Germ.muts__ALL.png
        │   │   ├── ATM.Germ.muts__CH.png
        │   │   ├── CH_muts__Functional_or_not_2022-09-09.png
        │   │   ├── DATA_ATM.maf_CH__Germ.vars.txt
        │   │   ├── DATA_ATM.maf__Germ.vars.txt
        │   │   ├── Supp_Figure_13_2022-09-09.png
        │   │   └── wt.VS.g.ATM_age.bins2022-09-09.png
        │   ├── Supp_Fig.7
        │   │   ├── Supp_Figure_4_2022-09-09.png
        │   │   ├── TABLE_v2_ht_DDR_Gene_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.tsv
        │   │   ├── ht_mutExcl_DDR_CH.v.CH__2022-09-09.png
        │   │   └── v2_ht_DDR_Gene_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.png
        │   ├── Supp_Fig.8
        │   │   ├── Supp_Figure_5_2022-09-09.png
        │   │   ├── TABLE_ht_Gene.v.mCA.genes__(SNVs\ +\ CNVs)__glm__2022-09-09.tsv
        │   │   ├── TABLE_ht_by.mCAs_Main_associations__(Germ\ events)__glm__2022-09-09.tsv
        │   │   ├── TABLE_ht_mCA_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.tsv
        │   │   ├── ht_Gene.v.mCA.genes__(SNVs\ +\ CNVs)__glm__2022-09-09.png
        │   │   ├── ht_by.mCAs_Main_associations__(Germ\ events)__glm__2022-09-09.png
        │   │   ├── ht_mCA_by_Gene__(SNVs\ +\ CNVs)__glm__2022-09-09.png
        │   │   ├── mCA_summary_2022-09-09.png
        │   │   ├── p_GERMLINE_mCA_Main_associations2022-09-09.png
        │   │   ├── venn_CH.v.Germline.v.mCA_2022-09-09.png
        │   │   └── venn_CH.v.mCA_2022-09-09.png
        │   ├── Supp_Fig.9
        │   │   ├── Supp_Figure_11_2022-09-09.png
        │   │   ├── letters_Supp_Figure_11_2022-09-09.png
        │   │   ├── leuk_positive
        │   │   │   ├── All_groups
        │   │   │   │   ├── ATM_sur.plot_2022-09-09.png
        │   │   │   │   ├── BRCA1_sur.plot_2022-09-09.png
        │   │   │   │   ├── BRCA2_sur.plot_2022-09-09.png
        │   │   │   │   ├── CHEK2_sur.plot_2022-09-09.png
        │   │   │   │   └── TP53_sur.plot_2022-09-09.png
        │   │   │   └── Selected_groups
        │   │   │       ├── ATM_sur.plot_GroupsSelected_2022-09-09.png
        │   │   │       ├── BRCA1_sur.plot_GroupsSelected_2022-09-09.png
        │   │   │       ├── BRCA2_sur.plot_GroupsSelected_2022-09-09.png
        │   │   │       ├── CHEK2_sur.plot_GroupsSelected_2022-09-09.png
        │   │   │       └── TP53_sur.plot_GroupsSelected_2022-09-09.png
        │   │   └── sur.plot_All.Variables_2022-09-09.png
        │   └── Tables
        │       ├── Table_90genes_germline_2022-09-08.html
        │       ├── Table_90genes_germline_2022-09-08.png
        │       └── Table_90genes_germline_2022-09-08.tsv
        └── script_all_figures.R
