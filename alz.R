# [1] Systematic review and meta-analysis of human transcriptomics reveals neuroinflammation, deficient energy metabolism, and proteostasis failure across neurodegeneration

# [1a] https://data.mendeley.com/datasets/752nd4w7pd/2

# rm(list = ls(pattern = 'temp.*|test.*'))
source("functions_for_data_prep.R")
source("functions_for_annotation.R")
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/bioinfo_little_helpers.R')
library(magrittr)
library(Hmisc)

save(alz, file = "als")


#########################
### ALZHEIMER DATASET ###

alz <- list()

###################
### GET INPUT 1 ###

alz$inputAlzData$preFull <- purrr::map_dfr(
  .x = list.dirs(path = "Alzheimer's Disease/", full.names = F, recursive = F), 
  .f = function(dataset){
    
    files <- purrr::map_dfr(
      .x = list.files(path = paste0("Alzheimer's Disease/", dataset), pattern =  ".*csv"), 
      .f = function(file_){
        
        file_loaded <- read.csv(file = paste0("Alzheimer's Disease/", dataset, "/", file_))
        
        file_loaded$Probe <- as.character(file_loaded$Probe)
        
        file_loaded$exp_comp <- paste0(dataset, "__",stringr::str_remove(string = file_, pattern = "\\.csv") )
        
        return(file_loaded)
      })
    
    return(files)
  })

alz$inputAlzData$full <- alz$inputAlzData$preFull

alz$inputAlzData$full <- subset(x = alz$inputAlzData$full, subset = alz$inputAlzData$full$pVal < 0.05)

### GET INPUT 1 ###
###################



########################
### GET DESCRIPTIONS ###

alz$inputAlzData$preFullDesc <- data.frame("exp_comp" = unique(alz$inputAlzData$preFull$exp_comp))

alz$inputAlzData$preFullDesc$platform <- stringr::str_extract(string = alz$inputAlzData$preFullDesc$exp_comp, pattern = "GPL.*")

alz$inputAlzData$preFullDesc$species <- "homo"

alz$inputAlzData$preFullDesc$tissue <- stringr::str_extract(string = alz$inputAlzData$preFullDesc$exp_comp, pattern = " - .*,")
alz$inputAlzData$preFullDesc$tissue <- stringr::str_remove(string = alz$inputAlzData$preFullDesc$tissue, pattern = ",")
alz$inputAlzData$preFullDesc$tissue <- stringr::str_remove(string = alz$inputAlzData$preFullDesc$tissue, pattern = " - ")

alz$inputAlzData$preFullDesc$exp_comp <- stringr::str_remove(string = alz$inputAlzData$preFullDesc$exp_comp, pattern = ".*__")
alz$inputAlzData$preFullDesc$exp_comp <- stringr::str_remove(string = alz$inputAlzData$preFullDesc$exp_comp, pattern = ",.*")
alz$inputAlzData$preFullDesc$exp_comp <- stringr::str_replace(string = alz$inputAlzData$preFullDesc$exp_comp, pattern = " DEGs - ", replacement = "__")

alz$inputAlzData$preFullDesc$exp <- as.numeric(c(1:length(alz$inputAlzData$preFullDesc[[1]]))) ### !!! This is actually comparison. Does it matter?

temp <- read.delim("data/tissueDesc.tsv")

alz$inputAlzData$preFullDesc$tissueGroup <- recode_values_based_on_key(
  to_recode_chrvec = alz$inputAlzData$preFullDesc$tissue, 
  replace_this_chrvec = temp$tissue, 
  with_this_chrvec = temp$tissueGroup)

### GET DESCRIPTIONS ###
########################



###################
### GET INPUT 2 ###

alz$inputAlzData$full <- dplyr::select(.data = alz$inputAlzData$full, Probe, Entrez, Symbol, logFC, exp_comp)

alz$inputAlzData$full$exp_comp <- stringr::str_remove(string = alz$inputAlzData$full$exp_comp, pattern = ".*__")
alz$inputAlzData$full$exp_comp <- stringr::str_remove(string = alz$inputAlzData$full$exp_comp, pattern = ",.*")
alz$inputAlzData$full$exp_comp <- stringr::str_replace(string = alz$inputAlzData$full$exp_comp, pattern = " DEGs - ", replacement = "__")

alz$inputAlzData$full$exp <- as.numeric( recode_values_based_on_key(
  to_recode_chrvec = alz$inputAlzData$full$exp_comp, 
  replace_this_chrvec = alz$inputAlzData$preFullDesc$exp_comp, 
  with_this_chrvec = alz$inputAlzData$preFullDesc$exp)
) ### !!! This is actually comparison. Does it matter?

### GET INPUT 2 ###
###################







##############################################
### IDENTIFY PLATFORM FOR PROBE ANNOTATION ###

alz$inputAlzData[['highest_hit_analysis']] <- master_annotator(
  descriptions_df = alz$inputAlzData$preFullDesc, 
  exp_id_col = "exp", 
  des_species_col = "species", 
  input_df = alz$inputAlzData$full, 
  input_id_col = 'Probe',
  PERFORM_B_get_the_highest_hit_returning_id_type = T,
  B_highest_hit_int_Probe_IDs_to_test = 200) ### !!! Tutaj wykorzystuję exp, które jest comp

# highest_hit_analysis <- alz$inputAlzData[['highest_hit_analysis']]
# alz$inputAlzData[['highest_hit_analysis']] <- highest_hit_analysis

### REMOVE EXPERIMENTS WITH IMPROPERLY DETECTED PLATFORMS ###
# exp 3 = gene symbols
# exp 11-13 = transcript id
# exp 14-15 = no microarray in ensembl
# exp 19 = no microarray in ensembl (PLATFORM DATA: GPL2700, Sentrix HumanRef-8 Expression BeadChip, GI_10047089-S; DATA: GI_4585642-S;  FILTER: illumina_humanref_8_v3, ILMN_2416019)
# exp 21 - different probe numbers (PLATFORM DATA: GPL6255, Illumina humanRef-8 v2.0, ILMN_10000; DATA: ILMN_9214, FILTER: illumina_humanref_8_v3)

alz$inputAlzData$improperly_detected_exp <- c(2, 3, 11:15, 19, 21, 25, 26, 30, 31, 40)

alz$inputAlzData$exp_for_probe_analysis$platform_ids_for_probe_annotation <- subset(
  x = alz$inputAlzData[['highest_hit_analysis']][['best_ID']],
  subset = alz$inputAlzData[['highest_hit_analysis']][['best_ID']][['Exp_ID']] %nin%
    alz$inputAlzData[['improperly_detected_exp']])

alz$inputAlzData$exp_for_probe_analysis$platform_ids_for_probe_annotation$platform_to_use <- as.character(alz$inputAlzData$exp_for_probe_analysis$platform_ids_for_probe_annotation[['platform_to_use']])

alz$inputAlzData$exp_for_probe_analysis$input_proper <- subset(
  x = alz$inputAlzData$full, 
  subset = alz$inputAlzData$full[["exp"]] %in% 
    alz$inputAlzData$exp_for_probe_analysis[['platform_ids_for_probe_annotation']][['Exp_ID']])

alz$inputAlzData$exp_NOT_for_probe_analysis$input_proper <- rbind(
  alz$inputAlzData$exp_NOT_for_probe_analysis$input,
  subset(
    x = alz$inputAlzData$full, 
    subset = alz$inputAlzData$full[["exp"]] %nin% 
      alz$inputAlzData$exp_for_probe_analysis[['platform_ids_for_probe_annotation']][['Exp_ID']])
)

length(alz$inputAlzData$exp_NOT_for_probe_analysis$input[[1]]) + length(alz$inputAlzData$exp_for_probe_analysis$input[[1]]) ==
  length(alz$inputAlzData$exp_NOT_for_probe_analysis$input_proper[[1]]) + length(alz$inputAlzData$exp_for_probe_analysis$input_proper[[1]])



alz$inputAlzData$input <- list(
  'data_for_probe_annotation' = alz$inputAlzData$exp_for_probe_analysis$input_proper,
  'platform_ids_for_probe_annotation' = alz$inputAlzData$exp_for_probe_analysis$platform_ids_for_probe_annotation$platform_to_use,
  'data_NOT_for_probe_annotation' = alz$inputAlzData$exp_NOT_for_probe_analysis$input_proper
)

alz$inputAlzData$input$qa$data_for_probe_annotation <- verify_df(df_ = alz$inputAlzData$input[['data_for_probe_annotation']], only_qa = T)
alz$inputAlzData$input$qa$data_NOT_for_probe_annotation <- verify_df(df_ = alz$inputAlzData$input[['data_NOT_for_probe_annotation']], only_qa = T)

### IDENTIFY PLATFORM FOR PROBE ANNOTATION ###
##############################################



################
### ANNOTATE ###

alz$inputAlzData$probes_direct$annotations <- master_annotator(
  descriptions_df = alz[["inputAlzData"]][["preFullDesc"]],
  exp_id_col = "exp",
  des_species_col = "species",
  input_df = alz$inputAlzData$input[['data_for_probe_annotation']],
  input_id_col = 'Probe',
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = alz$inputAlzData$input[['platform_ids_for_probe_annotation']],
  A_D_C_string_separator__ = ', ')



alz$inputAlzData$probes_direct$finalized <- subset(x = alz$inputAlzData$probes_direct$annotations, subset = !is.na(alz$inputAlzData$probes_direct[['annotations']][['external_gene_name']]))

alz$inputAlzData$probes_direct$leftover <- subset(x = alz$inputAlzData$probes_direct$annotations, subset = is.na(alz$inputAlzData$probes_direct[['annotations']][['external_gene_name']]))



length(alz$inputAlzData$input[['data_for_probe_annotation']][[1]]) == length(alz$inputAlzData$probes_direct[['annotations']][[1]]) 

check_was_the_spliting_of_df_by_filtering_ok(
  df_original = alz$inputAlzData$probes_direct[['annotations']], 
  list_df_splited = list(
    alz$inputAlzData$probes_direct[['finalized']], 
    alz$inputAlzData$probes_direct[['leftover']]))



alz$inputAlzData$identifers$input <- rbind(
  dplyr::select(alz$inputAlzData$probes_direct[['leftover']], -'external_gene_name'),
  alz$inputAlzData$input[['data_NOT_for_probe_annotation']])

alz$inputAlzData$identifers$input$Entrez <- as.character(alz$inputAlzData$identifers$input$Entrez)



alz$inputAlzData$identifers$annotations <- master_annotator(
  descriptions_df = alz[["inputAlzData"]][["preFullDesc"]],
  exp_id_col = "exp",
  des_species_col = "species",
  input_df = alz$inputAlzData$identifers$input,
  input_id_col = "Entrez",
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = "entrezgene_id",
  A_D_C_string_separator__ = ', ')



alz$inputAlzData$identifers$finalized <- subset(x = alz$inputAlzData$identifers$annotations, subset = !is.na(alz$inputAlzData$identifers[['annotations']][['external_gene_name']]))

alz$inputAlzData$identifers$leftover <- subset(x = alz$inputAlzData$identifers$annotations, subset = is.na(alz$inputAlzData$identifers[['annotations']][['external_gene_name']]))



length(alz$inputAlzData$identifers$input[[1]]) == length(alz$inputAlzData$identifers[['annotations']][[1]]) 

check_was_the_spliting_of_df_by_filtering_ok(
  df_original = alz$inputAlzData$identifers[['annotations']], 
  list_df_splited = list(
    alz$inputAlzData$identifers[['finalized']], 
    alz$inputAlzData$identifers[['leftover']]))



alz$inputAlzData$identifers[['leftover']]$external_gene_name <- tolower(alz$inputAlzData$identifers[['leftover']]$Symbol)



alz$inputAlzData$annotated <- rbind(
  alz$inputAlzData$probes_direct$finalized,
  alz$inputAlzData$identifers$finalized,
  alz$inputAlzData$identifers$leftover
)

length(alz$inputAlzData$full[[1]]) == length(alz$inputAlzData$annotated[[1]]) 

alz$inputAlzData$annotated$bestName <- purrr::map_chr(.x = alz$inputAlzData$annotated$external_gene_name, .f = select_best_geneName_wrapper_for_single_string)

### ANNOTATE ###
################


alz$inputAlzData$toUse <- dplyr::select(.data = alz$inputAlzData$annotated, Symbol = bestName, logFC, exp_comp)

alz$inputAlzData$toUse <- alz$inputAlzData$toUse %>%
  dplyr::group_by(exp_comp, Symbol) %>%
  dplyr::summarize( logFC_median = median(logFC) )

alz$inputAlzData$toUseWhole <- get_spread_data_and_entries_with_n_or_more_values_per_exp(
  dataset_ = alz$inputAlzData$toUse, 
  spread_key = "exp_comp", 
  spread_value = "logFC_median",
  n_ = 3,
  dataset_name = 'inputAlzData_full', 
  extract_exp_from_exp_and_comp_ = '__',
  gene_name_col_ = "Symbol",
  dir_ = "no",
  inside_dir_name = 'datasets',
  save = F) ### !!! ok, tu jest problem, bo dostajemy dane z 3 albo więcej porównań, nie papierów. Ale czy to ma znaczenie, jeśli klastrowanie nie jest tutaj naszym celem? Also - ta funkcja źle się nazywa. Nie ma znaczenia, bo potem i tak korzystamy tylko z tabelki whole

alz$alzData$toUseWholeNb <- get_number_and_percentage_of_directionality_of_exp_and_comp_first_column_names(
  spread_df_ = alz$inputAlzData$toUseWhole$spread_work_dataset,
  name_of_df_ = "nah",
  gene_col = "Symbol",
  dir_ = "xxx",
  inside_dir_name = 'xxx',
  exp_col = "exp",
  logfc_col = 'logFC_median',
  pattern_to_extract_exp_from_exp_and_comp_ = '__.*')

alz$alzData$toUseWholeNb$pTotal <- NA
alz$alzData$toUseWholeNb$pUp <- NA
alz$alzData$toUseWholeNb$ranking <- NA

for (rowNb in seq_along(alz$alzData$toUseWholeNb[[1]])) {
  alz$alzData$toUseWholeNb$pTotal[[rowNb]] <- binom.test(
    x = alz$alzData$toUseWholeNb$no_of_comps[[rowNb]], 
    n = length( unique(alz$inputAlzData$toUseWhole$work_dataset$exp_comp) ),
    p = 0.05, 
    alternative = "greater")$p.value
  
  alz$alzData$toUseWholeNb$pUp[[rowNb]] <- binom.test(
    x = as.integer(alz$alzData$toUseWholeNb$no_of_comps[[rowNb]] * (alz$alzData$toUseWholeNb$perc_of_upregulated[[rowNb]] / 100)), 
    n = alz$alzData$toUseWholeNb$no_of_comps[[rowNb]],
    p = 0.5, 
    alternative = "two.sided")$p.value
  
  alz$alzData$toUseWholeNb$ranking[[rowNb]] <- ifelse(
    test = alz$alzData$toUseWholeNb$perc_of_upregulated[[rowNb]] >= 50, 
    yes = (alz$alzData$toUseWholeNb$perc_of_upregulated[[rowNb]] * alz$alzData$toUseWholeNb$number_of_exp[[rowNb]]) / 100, 
    no = ( (100 - alz$alzData$toUseWholeNb$perc_of_upregulated[[rowNb]]) * alz$alzData$toUseWholeNb$number_of_exp[[rowNb]] ) / 100
           )
}

alz$alzData$toUseWholeNb$padjTotal <- p.adjust(p = alz$alzData$toUseWholeNb$pTotal, method = "bonferroni")
alz$alzData$toUseWholeNb$padjUp <- p.adjust(p = alz$alzData$toUseWholeNb$pUp, method = "bonferroni")


alz$alzData$toUseWholeNb <- subset(alz$alzData$toUseWholeNb, alz$alzData$toUseWholeNb$pTotal < 0.05)



### ALZHEIMER DATASETS ###
##########################



































###############
### GC DATA ###

alz$inputCompData$gc <- openxlsx::read.xlsx("data/FULL LIST Supplementary data 1 updated list_AMS_v2.xlsx", sheet = "data")

######################
### UPDATE GC DATA ###


###########
### DVV ###
alz$updates$gc$dvv$input <- openxlsx::read.xlsx("data/Glucocorticoids prime the inflammatory response of human hippocampal cells through up-regulation of inflammatory pathways/DVV vs VVV FC 1.2/DVV vs VVV FC 1.2_lista per IPA.xlsx", colNames = F)

colnames(alz$updates$gc$dvv$input) <- c("Symbol", "p", "fc")

alz$updates$gc$dvv$input$Symbol <- tolower(alz$updates$gc$dvv$input$Symbol)
alz$updates$gc$dvv$input$dummy_exp <- 1

alz$updates$gc$dvv$dvv_desc_dummy <- data.frame("Species" = "homo", "dummy_exp" = 1)


alz$updates$gc$dvv$ncbi_annotation_of_symbols_to_gene_ids <- master_annotator(
  descriptions_df = alz$updates$gc$dvv$dvv_desc_dummy,
  exp_id_col = "dummy_exp",
  des_species_col = "Species",
  input_df = alz$updates$gc$dvv$input,
  input_id_col = 'Symbol',
  C_legacy_D_str_identifier_type__ = 'Gene name',
  PERFORM_D_ncbi_annotation = T,
  A_D_C_string_separator__ = ", ")

alz$updates$gc$dvv$input_plus_new_gene_id_from_ncbi_column <- master_annotator(
  descriptions_df = alz$updates$gc$dvv$dvv_desc_dummy,
  exp_id_col ="dummy_exp",
  des_species_col = "Species",
  input_df = alz$updates$gc$dvv$input,
  input_id_col = 'Symbol',
  PERFORM_E_add_new_gene_id_col_originating_from_ncbi = T,
  E_PERFORM_D_output = alz$updates$gc$dvv[['ncbi_annotation_of_symbols_to_gene_ids']],
  A_D_C_string_separator__ = ", ")

alz$updates$gc$dvv$pseudomemoize_external_gene_names_to_gene_id <- master_annotator(
  descriptions_df = alz$updates$gc$dvv$dvv_desc_dummy,
  exp_id_col = 'dummy_exp',
  des_species_col = "Species",
  input_df = alz$updates$gc$dvv[['ncbi_annotation_of_symbols_to_gene_ids']],
  input_id_col = 'Gene_ID',
  A_C_all_Gene_ID_col = 'Gene_ID',
  PERFORM_A_should_i_prepare_dbs_for_pseudomemoization = T,
  A_return_qa_of_pseudomemoization = F,
  A_D_C_string_separator__ = ", ")

alz$updates$gc$dvv$annotationsList <- master_annotator(
  descriptions_df = alz$updates$gc$dvv$dvv_desc_dummy,
  exp_id_col = 'dummy_exp',
  des_species_col = "Species",
  input_df = alz$updates$gc$dvv[['input_plus_new_gene_id_from_ncbi_column']],
  input_id_col = 'Gene_ID_from_ncbi', ### from probe
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = 'Gene_ID',
  A_C_all_Gene_ID_col = 'Gene_ID_from_ncbi',
  C_pseudo_memoized_db = alz$updates$gc$dvv[['pseudomemoize_external_gene_names_to_gene_id']],
  A_D_C_string_separator__ = ", ")

alz$updates$gc$dvv$annotationsDF <- rlist::list.rbind(alz$updates$gc$dvv$annotationsList)

alz$updates$gc$dvv$annotationsDF$direction <- ifelse(
  test = alz$updates$gc$dvv$annotationsDF$fc > 0, 
  yes = "up", 
  no = "down")

alz$updates$gc$dvv$output <- data.frame("Symbol" = NA, "dvvDir" = alz$updates$gc$dvv$annotationsDF$direction)

alz$updates$gc$dvv$output$Symbol <- ifelse(
  test = !is.na(alz$updates$gc$dvv$annotationsDF$external_gene_name), 
  yes = tolower(alz$updates$gc$dvv$annotationsDF$external_gene_name), 
  no = alz$updates$gc$dvv$annotationsDF$Symbol)

alz$updates$gc$dvv$output$Symbol <- purrr::map_chr(.x = alz$updates$gc$dvv$output$Symbol, .f = select_best_geneName_wrapper_for_single_string)

### DVV ###
###########






#################
### GSE132423 ###

alz$updates$gc$GSE132423$input <- readr::read_tsv("data/GSE132423.tsv")

alz$updates$gc$GSE132423$input <- subset(
  alz$updates$gc$GSE132423$input, 
  subset = alz$updates$gc$GSE132423$input$P.Value < 0.05, 
  select = c("ID", "logFC", "Gene.symbol", "Gene.ID"))

alz$updates$gc$GSE132423$input$dummy_exp <- 1



alz$updates$gc$GSE132423$GSE13_desc_dummy <- data.frame("Species" = "mus", "dummy_exp" = 1)



alz$updates$gc$GSE132423$probes_direct$annotations <- master_annotator(
  descriptions_df = alz$updates$gc$GSE132423$GSE13_desc_dummy,
  exp_id_col = "dummy_exp",
  des_species_col = "Species",
  input_df = alz$updates$gc$GSE132423$input,
  input_id_col = 'ID',
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = "illumina_mousewg_6_v2",
  A_D_C_string_separator__ = ', ')


alz$updates$gc$GSE132423$probes_direct$finalized <- subset(x = alz$updates$gc$GSE132423$probes_direct$annotations, subset = !is.na(alz$updates$gc$GSE132423$probes_direct[['annotations']][['external_gene_name']]))

alz$updates$gc$GSE132423$probes_direct$leftover <- subset(x = alz$updates$gc$GSE132423$probes_direct$annotations, subset = is.na(alz$updates$gc$GSE132423$probes_direct[['annotations']][['external_gene_name']]))



length(alz$updates$gc$GSE132423$input[[1]]) == length(alz$updates$gc$GSE132423$probes_direct[['annotations']][[1]]) 

check_was_the_spliting_of_df_by_filtering_ok(
  df_original = alz$updates$gc$GSE132423$probes_direct[['annotations']], 
  list_df_splited = list(
    alz$updates$gc$GSE132423$probes_direct[['finalized']], 
    alz$updates$gc$GSE132423$probes_direct[['leftover']]))



alz$updates$gc$GSE132423$identifers$input <- dplyr::select(alz$updates$gc$GSE132423$probes_direct[['leftover']], -'external_gene_name')

alz$updates$gc$GSE132423$identifers$input$Entrez <- as.character(alz$updates$gc$GSE132423$identifers$input$Gene.ID)




alz$updates$gc$GSE132423$identifers$annotations <- master_annotator(
  descriptions_df = alz$updates$gc$GSE132423$GSE13_desc_dummy,
  exp_id_col = "dummy_exp",
  des_species_col = "Species",
  input_df = alz$updates$gc$GSE132423$identifers$input,
  input_id_col = "Entrez",
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = "entrezgene_id",
  A_D_C_string_separator__ = ', ')



alz$updates$gc$GSE132423$identifers$finalized <- subset(x = alz$updates$gc$GSE132423$identifers$annotations, subset = !is.na(alz$updates$gc$GSE132423$identifers[['annotations']][['external_gene_name']]))

alz$updates$gc$GSE132423$identifers$leftover <- subset(x = alz$updates$gc$GSE132423$identifers$annotations, subset = is.na(alz$updates$gc$GSE132423$identifers[['annotations']][['external_gene_name']]))



length(alz$updates$gc$GSE132423$identifers$input[[1]]) == length(alz$updates$gc$GSE132423$identifers[['annotations']][[1]]) 

check_was_the_spliting_of_df_by_filtering_ok(
  df_original = alz$updates$gc$GSE132423$identifers[['annotations']], 
  list_df_splited = list(
    alz$updates$gc$GSE132423$identifers[['finalized']], 
    alz$updates$gc$GSE132423$identifers[['leftover']]))



alz$updates$gc$GSE132423$identifers[['leftover']]$external_gene_name <- tolower(alz$updates$gc$GSE132423$identifers[['leftover']]$Gene.symbol)



alz$updates$gc$GSE132423$annotated <- rbind(
  alz$updates$gc$GSE132423$probes_direct$finalized,
  dplyr::select(alz$updates$gc$GSE132423$identifers$finalized, -Entrez),
  dplyr::select(alz$updates$gc$GSE132423$identifers$leftover, -Entrez)
)

length(alz$updates$gc$GSE132423$input[[1]]) == length(alz$updates$gc$GSE132423$annotated[[1]]) 

alz$updates$gc$GSE132423$annotated$bestName <- purrr::map_chr(.x = alz$updates$gc$GSE132423$annotated$external_gene_name, .f = select_best_geneName_wrapper_for_single_string)

alz$updates$gc$GSE132423$output <- alz$updates$gc$GSE132423$annotated

alz$updates$gc$GSE132423$output <- subset(x = alz$updates$gc$GSE132423$output, subset = !is.na(alz$updates$gc$GSE132423$output$bestName))

alz$updates$gc$GSE132423$output$GSE13Dir <- ifelse(
  test = alz$updates$gc$GSE132423$output$logFC > 0, 
  yes = "up", 
  no = "down")

alz$updates$gc$GSE132423$output <- dplyr::select(alz$updates$gc$GSE132423$output, Symbol = bestName, GSE13Dir)


alz$updates$gc$GSE132423$output <- tidyr::pivot_wider(data = alz$updates$gc$GSE132423$output, names_from = GSE13Dir, values_from = GSE13Dir)

alz$updates$gc$GSE132423$output$GSE13Dir <- ifelse(
  test = alz$updates$gc$GSE132423$output$up != "NULL" & alz$updates$gc$GSE132423$output$down != "NULL",
  yes = "ambi",
  no = ifelse(
    test = alz$updates$gc$GSE132423$output$up != "NULL" & alz$updates$gc$GSE132423$output$down == "NULL", 
    yes = "up", 
    no = "down"))

alz$updates$gc$GSE132423$output <- subset(
  x = alz$updates$gc$GSE132423$output, 
  subset = alz$updates$gc$GSE132423$output$GSE13Dir != "ambi",
  select = c("Symbol", "GSE13Dir"))

### GSE132423 ###
#################






#################
### GSE128148 ###

alz$updates$gc$GSE128148$input <- readr::read_tsv("data/GSE128148_dea_veh_cort.tsv")

alz$updates$gc$GSE128148$input$dummy_exp <- as.numeric(1)

alz$updates$gc$GSE128148$identifers$annotations <- master_annotator(
  descriptions_df = data.frame("Species" = "mus", "dummy_exp" = 1),
  exp_id_col = "dummy_exp",
  des_species_col = "Species",
  input_df = alz$updates$gc$GSE128148$input,
  input_id_col = "genes",
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = "ensembl_gene_id",
  A_D_C_string_separator__ = ', ')

alz$updates$gc$GSE128148$identifers$annotations$Symbol <- ifelse(
  test = is.na(alz$updates$gc$GSE128148$identifers$annotations$external_gene_name), 
  yes = tolower(alz$updates$gc$GSE128148$identifers$annotations$genes), 
  no = alz$updates$gc$GSE128148$identifers$annotations$external_gene_name)

alz$updates$gc$GSE128148$output <- subset(
  x = alz$updates$gc$GSE128148$identifers$annotations, 
  subset = alz$updates$gc$GSE128148$identifers$annotations$PValue < 0.05)

alz$updates$gc$GSE128148$output$GSE12Dir <- ifelse(
  test = alz$updates$gc$GSE128148$output$logFC > 0, 
  yes = "up", 
  no = "down")

alz$updates$gc$GSE128148$output <- dplyr::select(.data = alz$updates$gc$GSE128148$output, Symbol, GSE12Dir)

### GSE128148 ###
#################


### We are using all.x = T here, because if given gene was not found in original analysis, than it will have at most 3 hits, and that is not enough to be of interest
alz$updates$gc$merged <- purrr::reduce(
  .x = list(
    alz[["inputCompData"]][["gc"]], 
    alz$updates$gc$dvv$output, 
    alz$updates$gc$GSE132423$output, 
    alz$updates$gc$GSE128148$output), 
  .f = merge, 
  by = "Symbol",
  all.x = T)

for (rowNb in seq_along(alz$updates$gc$merged[[1]])) {
  
  alz$updates$gc$merged$updatedTotal[[rowNb]] <- sum(
    !is.na( c(alz$updates$gc$merged$dvvDir[[rowNb]], alz$updates$gc$merged$GSE13Dir[[rowNb]], alz$updates$gc$merged$GSE12Dir[[rowNb]]) )
  )
  
  alz$updates$gc$merged$updatedUp[[rowNb]] <- sum(
    c(alz$updates$gc$merged$dvvDir[[rowNb]] == "up", alz$updates$gc$merged$GSE13Dir[[rowNb]] == "up", alz$updates$gc$merged$GSE12Dir[[rowNb]] == "up"),
    na.rm = T
  )
  
  alz$updates$gc$merged$updatedDown[[rowNb]] <- sum(
    c(alz$updates$gc$merged$dvvDir[[rowNb]] == "down", alz$updates$gc$merged$GSE13Dir[[rowNb]] == "down", alz$updates$gc$merged$GSE12Dir[[rowNb]] == "down"),
    na.rm = T
  )
}

alz$updates$gc$merged$updatedTotal <- as.numeric(alz$updates$gc$merged$updatedTotal)
alz$updates$gc$merged$updatedUp <- as.numeric(alz$updates$gc$merged$updatedUp)
alz$updates$gc$merged$updatedDown <- as.numeric(alz$updates$gc$merged$updatedDown)






alz$updates$gc$merged$gc_down[is.na(alz$updates$gc$merged$gc_down)] <- 0

alz$updates$gc$merged$gc_down <- ifelse(
  test = alz$updates$gc$merged$gc_up == "down", 
  yes = alz$updates$gc$merged$gc_down + 1, 
  no = alz$updates$gc$merged$gc_down)


alz$updates$gc$merged$gc_up[alz$updates$gc$merged$gc_up == "up"] <- 1

alz$updates$gc$merged$gc_up[alz$updates$gc$merged$gc_up %in% c("down", "ambig")] <- 0

alz$updates$gc$merged$gc_up <- as.numeric(alz$updates$gc$merged$gc_up)



alz$updates$gc$merged$gc_totalPubs <- alz$updates$gc$merged$gc_totalPubs + alz$updates$gc$merged$updatedUp + alz$updates$gc$merged$updatedDown

alz$updates$gc$merged$gc_up <- alz$updates$gc$merged$gc_up + alz$updates$gc$merged$updatedUp

alz$updates$gc$merged$gc_down <- alz$updates$gc$merged$gc_down + alz$updates$gc$merged$updatedDown






alz$compData$gc <- alz$updates$gc$merged

alz$compData$gc <- dplyr::filter(alz$compData$gc, !(alz$compData$gc$Symbol == "cd47" & alz$compData$gc$gc_up == 0)) ### Był 1 błąd w porpzednich danych

alz$compData$gc <- subset(
  x = alz$compData$gc, 
  subset = alz$compData$gc$gc_totalPubs > 2, 
  select = c("Symbol", "gc_totalPubs", "gc_up"))

alz$compData$gc$gc_up <- (alz$compData$gc$gc_up / alz$compData$gc$gc_totalPubs) * 100

alz$compData$gc$gc_pTotal <- as.numeric(NA)
alz$compData$gc$gc_pUp <- as.numeric(NA)
alz$compData$gc$gc_ranking <- as.numeric(NA)

for (rowNb in seq_along(alz$compData$gc[[1]])) {
  alz$compData$gc$gc_pTotal[[rowNb]] <- binom.test(
    x = alz$compData$gc$gc_totalPubs[[rowNb]], 
    n = 20,
    p = 0.05, 
    alternative = "greater")$p.value
  
  alz$compData$gc$gc_pUp[[rowNb]] <- binom.test(
    x = as.integer(alz$compData$gc$gc_totalPubs[[rowNb]] * (alz$compData$gc$gc_up[[rowNb]] / 100)), 
    n = alz$compData$gc$gc_totalPubs[[rowNb]],
    p = 0.5, 
    alternative = "two.sided")$p.value
  
  alz$compData$gc$gc_ranking[[rowNb]] <- ifelse(
    test = alz$compData$gc$gc_up[[rowNb]] >= 50, 
    yes = (alz$compData$gc$gc_up[[rowNb]] * alz$compData$gc$gc_totalPubs[[rowNb]]) / 100, 
    no = ( (100 - alz$compData$gc$gc_up[[rowNb]]) * alz$compData$gc$gc_totalPubs[[rowNb]] ) / 100 )
}

### UPDATE GC DATA ###
######################

### GC DATA ###
###############









###############
### MERGING ###

alz$all <- purrr::map2(
  .x = alz$alzData, 
  .y = names(alz$alzData),
  .f = function(dataset, name){
    all <- purrr::reduce(
    .x = rlist::list.insert(alz$compData, 1, name = dataset), 
    .f = merge, 
    by = "Symbol", 
    all.x = T)
    
    all <- dplyr::select(.data = all, dplyr::everything(), -isch_p, -age1_mCorP, -age1_fCorP, -age1_mCorPadj, -age1_fCorPadj, -padjUp)
  
    all[all == "NA"] <- NA
  
    return(all)
})


alz$allFiltered <- purrr::map(
  .x = alz$all,
  .f = function(dataset){
    
    dplyr::filter(.data = dataset, dataset$gc_pTotal < 0.05 & dataset$padjTotal < 0.05)
  })

alz$allFiltered <- rlist::list.clean(alz$allFiltered, fun = function(x){ length(x[[1]])  == 0 })

  
openxlsx::write.xlsx(x = alz$allFiltered, file = "results.xlsx")


alz$other$significantAd <- purrr::map(
  .x = alz$all,
  .f = function(dataset){
    
    dplyr::filter(.data = dataset, dataset$padjTotal < 0.05)
  })

openxlsx::write.xlsx(x = alz$other$significantAd$toUseWholeNb, file = "significantAd_toUseWholeNb.xlsx")


### MERGING ###
###############





######################
### FOR ENRICHMENT ###
# 0.125 oznacza, że 4/4 są w jedną stronę

alz$other$enrich <- purrr::map2(
  .x = alz$allFiltered, 
  .y = names(alz$allFiltered),
  .f = function(dataset, name){
    
    dir = "enrichr/inputs/"
    
    all <- dplyr::select(.data = dataset, Symbol)
    readr::write_tsv(x = all, file = paste0(dir, name, ".txt"))
    
    return_ <- list("all" = all)
    
    return_ <- rlist::list.clean(.data = return_, fun = function(x){length(x[[1]]) == 0})

    return(return_)
  })

### FOR ENRICHMENT ###
######################











#############
### OTHER ###

alz$other$genesOfInterest$input <- subset(x = alz[["allFiltered"]][["toUseWholeNb"]], subset = alz[["allFiltered"]][["toUseWholeNb"]]$Symbol %in% c("aldoa", "arl4d", "btg1", "dab2", "dio2", "egr1", "ezr", "fam107a", "g3bp2", "gap43", "grina", "homer1", "hopx", "map1b", "mical2", "msn", "ndufb8", "nedd4l", "nrxn3", "nuak1", "prkx", "rcan2", "rgs4", "scamp2", "sesn1", "sparc", "timm17a", "tjp2", "trim36", "zfp36l1"), select = c("Symbol", "no_of_comps", "ranking", "gc_totalPubs", "gc_ranking", "perc_of_upregulated"))

openxlsx::write.xlsx(alz$other$genesOfInterest$input, "genesOfInterest.xlsx")



alz$other$nbOfAdRelevantGenes <- subset(x = alz$alzData$toUseWholeNb, subset = alz$alzData$toUseWholeNb$padjTotal < 0.05)
  


### ENRICHR PLAY ###

alz$other$compareTFs$toUseWholeNb <- read.delim(file = "enrichr/outputs_whole/toUseWholeNb_output.tsv", dec = ",")

alz$other$compareTFs$toUseWholeNb <- subset(alz$other$compareTFs$toUseWholeNb, subset = alz$other$compareTFs$toUseWholeNb$Database %in% c("ARCHS4_IDG_Coexp", "ARCHS4_Kinases_Coexp", "ARCHS4_TFs_Coexp", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "Enrichr_Submissions_TF-Gene_Coocurrence", "ESCAPE", "Genome_Browser_PWMs", "KEA_2015", "miRTarBase_2017", "NURSA_Human_Endogenous_Complexome", "PPI_Hub_Proteins", "TargetScan_microRNA", "TargetScan_microRNA_2017", "TF_Perturbations_Followed_by_Expression", "TF-LOF_Expression_from_GEO", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019"))

alz$other$compareTFs$toUseWholeNb$Term2 <- tolower( ifelse(
  test = alz$other$compareTFs$toUseWholeNb$Database %in% c("ARCHS4_IDG_Coexp", "ARCHS4_Kinases_Coexp", "ARCHS4_TFs_Coexp", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TF_Perturbations_Followed_by_Expression", "TF-LOF_Expression_from_GEO", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019"), 
  yes = stringr::str_remove(string = alz$other$compareTFs$toUseWholeNb$Term, pattern = " .*"),
  no =  ifelse(
    
    test = alz$other$compareTFs$toUseWholeNb$Database %in% c("ESCAPE"), 
    yes = stringr::str_remove(
      string = stringr::str_remove(
        string = alz$other$compareTFs$toUseWholeNb$Term, 
        pattern = "-.*"),
      pattern = "^CHiP |^Protein "), 
    no = ifelse(
      
      test = alz$other$compareTFs$toUseWholeNb$Database %in% c("Genome_Browser_PWMs"), 
      yes = stringr::str_remove(
        string = stringr::str_remove(
          string = alz$other$compareTFs$toUseWholeNb$Term, 
          pattern = ".*\\$"), 
        pattern = " .*"), 
      no = ifelse(
        
        test = alz$other$compareTFs$toUseWholeNb$Database %in% c("miRTarBase_2017", "TargetScan_microRNA_2017"), 
        yes = stringr::str_remove(
          string = stringr::str_remove(
            string = alz$other$compareTFs$toUseWholeNb$Term, 
            pattern = "-3p$|-5p$"), 
          pattern = "^[a-z]{3}-"),
        no = ifelse(
          
          test = alz$other$compareTFs$toUseWholeNb$Database %in% c("NURSA_Human_Endogenous_Complexome"), 
          yes = stringr::str_remove(
            string = stringr::str_remove(
              string = alz$other$compareTFs$toUseWholeNb$Term, 
              pattern = "\\)$"),
            pattern = ".*\\("), 
          no = ifelse(
            
            test = alz$other$compareTFs$toUseWholeNb$Database %in% c("TargetScan_microRNA"), 
            yes = stringr::str_remove(
                string = stringr::str_remove(
                  string = stringr::str_remove(
                    string = alz$other$compareTFs$toUseWholeNb$Term, 
                    pattern = "^[A-Z]{7},"),
                  pattern = ",.*"),
                pattern = "-3[p|P]$|-5[p|P]$"), 
            no = alz$other$compareTFs$toUseWholeNb$Term))))))
)



alz$other$compareTFs$toUseWholeNb_gc <- merge(x = alz$other$compareTFs$toUseWholeNb, y = alz$allFiltered$toUseWholeNb, by.x = "Term2", by.y = "Symbol")[,c(1:5, 8:17)]


### ENRICHR PLAY ###








### ALL GCS VS AD-RELATED GCS ENRICHMENT ###

alz$other$allGcVsAdRelatedGc$input_all <- dplyr::filter(.data = alz$compData$gc, alz$compData$gc$gc_pTotal < 0.05)

readr::write_tsv(
  x = as.data.frame( alz$other$allGcVsAdRelatedGc$all$Symbol ), 
  file = "enrichr/gc_nonAdSelected.txt",
  col_names = F)



alz$other$allGcVsAdRelatedGc$enrichr_all <- read.delim(file = "enrichr/gc_nonAdSelected/annotationenrichr.tsv")

alz$other$allGcVsAdRelatedGc$enrichr_all$Overlap_n <- as.numeric(stringr::str_remove(string = alz$other$allGcVsAdRelatedGc$enrichr_all$Overlap, pattern = " .*"))

alz$other$allGcVsAdRelatedGc$enrichr_all$n_of_genes <- length(alz$other$allGcVsAdRelatedGc$input_all[[1]])

alz$other$allGcVsAdRelatedGc$enrichr_all$bern_prob <- alz$other$allGcVsAdRelatedGc$enrichr_all$Overlap_n / alz$other$allGcVsAdRelatedGc$enrichr_all$n_of_genes

alz$other$allGcVsAdRelatedGc$enrichr_all$Term <- paste(alz$other$allGcVsAdRelatedGc$enrichr_all$Term, alz$other$allGcVsAdRelatedGc$enrichr_all$Database, sep = "__")

colnames(alz$other$allGcVsAdRelatedGc$enrichr_all)[2:13] <- stringr::str_c(colnames(alz$other$allGcVsAdRelatedGc$enrichr_all)[2:13], "_all")




alz$other$allGcVsAdRelatedGc$input_AdRelatedGc <- alz[["allFiltered"]][["toUseWholeNb"]]

alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc <- read.delim(file = "enrichr/outputs_whole/toUseWholeNb_output.tsv")

alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$Overlap_n <- as.numeric(stringr::str_remove(string = alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$Overlap, pattern = " .*"))

alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$n_of_genes <- length(alz$other$allGcVsAdRelatedGc$input_AdRelatedGc[[1]])

alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$bern_prob <- alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$Overlap_n / alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$n_of_genes

alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$Term <- paste(alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$Term, alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc$Database, sep = "__")

colnames(alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc)[2:13] <- stringr::str_c(colnames(alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc)[2:13], "_adgc")




alz$other$allGcVsAdRelatedGc$enrichr_merge <- merge(x = alz$other$allGcVsAdRelatedGc$enrichr_all, y = alz$other$allGcVsAdRelatedGc$enrichr_AdRelatedGc, by = "Term", all = T)

alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets <- subset(x = alz$other$allGcVsAdRelatedGc$enrichr_merge, subset = !is.na(alz$other$allGcVsAdRelatedGc$enrichr_merge$bern_prob_all) & !is.na(alz$other$allGcVsAdRelatedGc$enrichr_merge$bern_prob_adgc))
  
alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInoneDataset <- subset(x = alz$other$allGcVsAdRelatedGc$enrichr_merge, subset = is.na(alz$other$allGcVsAdRelatedGc$enrichr_merge$bern_prob_all) | is.na(alz$other$allGcVsAdRelatedGc$enrichr_merge$bern_prob_adgc))


alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets$p_difference_of_binom <- NA

for (rowNb in seq_along(alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets[[1]])) {
  alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets$p_difference_of_binom[[rowNb]] <- prop.test(x = c(alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets$Overlap_n_all[[rowNb]], alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets$Overlap_n_adgc[[rowNb]]),
            n = c(alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets$n_of_genes_all[[rowNb]], alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets$n_of_genes_adgc[[rowNb]]), 
            alternative = "two.sided")$p.value
}

alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasetsClean <- dplyr::select(alz$other$allGcVsAdRelatedGc$enrichr_mergeTermInBothDatasets, Term, Adjusted.P.value_all, bern_prob_all, Adjusted.P.value_adgc, bern_prob_adgc, p_difference_of_binom)






prop.test(x, n, p = NULL,
          alternative = c("two.sided", "less", "greater"),
          conf.level = 0.95, correct = TRUE)



fisher.test(rbind(c(1556,2455-1556), c(1671,2730-1671)), alternative="less")

### ALL GCS VS AD-RELATED GCS ENRICHMENT ###



### COMPARE CLUSTERS TO NEW WHOLE GENES ###

alz$other$clust_vs_wholeGenes <- purrr::map(
  .x = alz$other$genesOfInterest$clust_enrich,
  .f = function(cluster, name){
    
    merge(x = data.frame("Symbol" = na.exclude(cluster)), y = alz$allFiltered$toUseWholeNb, by = "Symbol")
    
  })

### COMPARE CLUSTERS TO NEW WHOLE GENES ###





### OTHER ###
#############




