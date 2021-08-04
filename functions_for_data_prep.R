# INPUT: single string containing multiple gene names separeted by separator
select_best_geneName_wrapper_for_single_string <- function(string_to_be_vectorised, separator = opts$str_sep)
{
  assertthat::assert_that(!is.na(string_to_be_vectorised))
  
  vectorised_string <- as.character(stringr::str_split(
    string = string_to_be_vectorised,
    pattern = separator,
    simplify = T
  ))
  
  the_best_name_ <- select_best_geneName(char_vec = vectorised_string)
  
  return(the_best_name_)
}


















# INPUT: char_vec - includes all gene names returned by ensembl for given gene. Hence, the function should be iterated over list of vectors, each vector/list element for single gene. regex_to_detect_bad_names_with - genes can have 4 numbers: Olr1237 OUTPUT: the_best_name - char_vec of length 1
select_best_geneName <- function(char_vec, regex_to_detect_bad_names_with = '(\\d{5})|(^gm\\d)', regex_to_detect_less_bad_names_with = '(rik)|(loc)|(gm)')
{
  char_vec <- tolower(char_vec)
  
  regex_to_detect_bad_names_with <- tolower(regex_to_detect_bad_names_with)
  
  regex_to_detect_less_bad_names_with <- tolower(regex_to_detect_less_bad_names_with)
  
  log_is_this_name_bad <- stringr::str_detect(
    string = char_vec, 
    pattern = regex_to_detect_bad_names_with)
  
  good_names <- subset(x = char_vec, subset = !log_is_this_name_bad)
  
  bad_names <- subset(x = char_vec, subset = log_is_this_name_bad)
  
  log_is_this_name_less_bad <- stringr::str_detect(
    string = bad_names,
    pattern = regex_to_detect_less_bad_names_with)
  
  less_bad_names <- subset(x = bad_names, subset = log_is_this_name_less_bad)
  
  if(length(good_names) != 0)
  {
    the_best_name <- good_names[[1]]
  }
  else if(length(less_bad_names) != 0)
  {
    the_best_name <- less_bad_names[[1]]
  }
  else if(length(bad_names) != 0)
  {
    the_best_name <- bad_names[[1]]
  }
  else
  {
    the_best_name <- NA
  }
  
  return(the_best_name)
}























get_subset_vector_for_entries_with_n_or_more_values_per_exp <- function(work_dataset_, exp_col, gene_name_col, n, extract_exp_from_exp_and_comp)
{
  work_dataset_[[exp_col]] <- stringr::str_remove(string = work_dataset_[[exp_col]], pattern = extract_exp_from_exp_and_comp)
  
  temp <- unique(subset(
    x = work_dataset_, 
    select = c(exp_col, gene_name_col)))

  temp$present <- T
  
  spread_temp <- tidyr::spread(
    data = temp, 
    key = eval(parse(text = exp_col)), 
    value = present)
  
  entries_with_n_or_more_values <- purrrlyr::by_row(
    .d = spread_temp[,-1], 
    .collate = "rows", 
    ..f = function(x) {
      sum(!is.na(x)) >= n}
  )$.out
  
  return(entries_with_n_or_more_values)
}























are_any_cols_empty_in_3_exps_up_data_wrapper <- function(spread_work_dataset_3_exps_up_)
{
  spread_work_dataset_3_exps_up_pseudomatrix <- spread_work_dataset_3_exps_up_[,2:ncol(spread_work_dataset_3_exps_up_)]
  
  spread_work_dataset_3_exps_up_pseudomatrix[!is.na(spread_work_dataset_3_exps_up_pseudomatrix)] <- 1
  
  spread_work_dataset_3_exps_up_pseudomatrix[is.na(spread_work_dataset_3_exps_up_pseudomatrix)] <- 0

  return_ <- t(
    purrr::map_df(
      .x = spread_work_dataset_3_exps_up_pseudomatrix, 
      .f = sum))
  
  logic_return_ <- return_ == 0
  
  return_ <- subset(return_, subset = logic_return_)
  
  return(return_)
}




















get_spread_data_and_entries_with_n_or_more_values_per_exp <- function(dataset_, spread_key, spread_value, extract_exp_from_exp_and_comp_ = '_.*', gene_name_col_ = opts[['ext_gene_name_col']], n_, dataset_name, dir_ = opts$data[['folder']], inside_dir_name = 'datasets', save = T)
{
  library(tidyr)

  work_dataset_ <- dataset_
  
  spread_work_dataset_ <- tidyr::spread(
    data = dataset_, 
    key = {{spread_key}}, 
    value = {{spread_value}})
  

  
  bool_entries_with_3_or_more_values_ <- get_subset_vector_for_entries_with_n_or_more_values_per_exp(
    work_dataset_ = work_dataset_,
    exp_col = spread_key,
    gene_name_col = gene_name_col_,
    n = n_,
    extract_exp_from_exp_and_comp = extract_exp_from_exp_and_comp_)
  
  spread_work_dataset_3_exps_up_ <- subset(
    x = spread_work_dataset_,
    subset = bool_entries_with_3_or_more_values_)

  work_dataset_3_exps_up_ <- subset(
    x = work_dataset_,
    subset = work_dataset_[[gene_name_col_]] %in% spread_work_dataset_3_exps_up_[[gene_name_col_]] )

  are_any_cols_empty_in_3_exps_up_data_ <- are_any_cols_empty_in_3_exps_up_data_wrapper(
    spread_work_dataset_3_exps_up_ = spread_work_dataset_3_exps_up_)
  
  spread_work_dataset_3_exps_up_ <- subset(
    x = spread_work_dataset_3_exps_up_, 
    select = colnames(spread_work_dataset_3_exps_up_) %nin%
      rownames(are_any_cols_empty_in_3_exps_up_data_))
  
  if (bazar::is.empty(work_dataset_3_exps_up_)) {
    work_dataset_3_exps_up_ <- NA
  }
  if (bazar::is.empty(spread_work_dataset_3_exps_up_)) {
    spread_work_dataset_3_exps_up_ <- NA
  }
  if (bazar::is.empty(are_any_cols_empty_in_3_exps_up_data_)) {
    are_any_cols_empty_in_3_exps_up_data_ <- NA
  }
  

  if (save == T) 
    {
    dir.create(paste0(dir_, '/', inside_dir_name))
    
    purrr::walk2(
      .x = list(work_dataset_, spread_work_dataset_, work_dataset_3_exps_up_, spread_work_dataset_3_exps_up_), 
      .y = list('work_dataset_', 'spread_work_dataset_', 'work_dataset_3_exps_up_', 'spread_work_dataset_3_exps_up_'),
      .f = function(dataset, name){
        write.table(
          x = dataset, 
          file = paste0(dir_, '/', inside_dir_name, '/', dataset_name, '_', name, '.tsv'), 
          sep = '\t', 
          dec = ',', 
          row.names = F)
      })
  }

  

  return(list('work_dataset' = dataset_, 'spread_work_dataset' = spread_work_dataset_, 'work_dataset_3_exps_up' = work_dataset_3_exps_up_, 'spread_work_dataset_3_exps_up' = spread_work_dataset_3_exps_up_, 'are_any_cols_empty_in_3_exps_up_data' = are_any_cols_empty_in_3_exps_up_data_))
}




















write_table_for_clustering <- function(dataset_, dataset_name, gene_name_col = opts[['ext_gene_name_col']], dir_ = opts$data[['folder']], inside_dir_name = 'for_clustering')
{
  dir.create(paste0(dir_, '/', inside_dir_name))
  
  for_clustering <- as.data.frame(dataset_)
  rownames(for_clustering) <- for_clustering[[gene_name_col]]
  for_clustering[[gene_name_col]] <- NULL
  for_clustering <- as.matrix(for_clustering)
  for_clustering[is.na(for_clustering)] <- 0
  
  
  
  write.table(
    x = as.data.frame(for_clustering), 
    file = paste0(dir_, '/', inside_dir_name, '/for_clustering_', dataset_name,'.tsv'),
    sep = '\t',
    row.names = T,
    col.names = NA,
    dec = '.'
  ) # For cluster-based hierachical clustering
  
  return(T)
}



















get_dataset_sets_wrapper <- function(subgroups_, dataset_set_name, whole_dataset)
{
  set_of_datasets <- list()
  
  for (subgroup_name in names(subgroups_)) {
    
    for (variable_name in names(subgroups_[[subgroup_name]])) 
      {
      dataset_name <- paste0(
        subgroup_name, 
        '_', 
        variable_name)
      
      if (!is.na(subgroups_[[subgroup_name]][[variable_name]])) 
        {
        if (!is.na(subgroups_[[subgroup_name]][[variable_name]][[dataset_set_name]][[1]])) 
          {
          set_of_datasets[[dataset_name]] <- subgroups_[[subgroup_name]][[variable_name]][[dataset_set_name]]
        }
        
      }
    }
  }
  
  set_of_datasets <- rlist::list.append(
    set_of_datasets, 
    'whole_dataset' = whole_dataset[[dataset_set_name]])
  
  return(set_of_datasets)
}







# get_dataset_sets_wrapper <- function(data_, dataset_set_name)
# {
#   set_of_datasets <- list()
#   
#   for (subgroup_name in names(data_[['subgroups']][['data']])) {
#     
#     for (variable_name in names(data_[['subgroups']][['data']][[subgroup_name]])) 
#     {
#       dataset_name <- paste0(
#         subgroup_name, 
#         '_', 
#         variable_name)
#       
#       if (!is.na(data_[['subgroups']][['data']][[subgroup_name]][[variable_name]])) 
#       {
#         if (!is.na(data_[['subgroups']][['data']][[subgroup_name]][[variable_name]][[dataset_set_name]][[1]])) 
#         {
#           set_of_datasets[[dataset_name]] <- data_[['subgroups']][['data']][[subgroup_name]][[variable_name]][[dataset_set_name]]
#         }
#         
#       }
#     }
#   }
#   
#   set_of_datasets <- rlist::list.append(
#     set_of_datasets, 
#     'whole_dataset' = data_$whole_dataset[[dataset_set_name]])
#   
#   return(set_of_datasets)
# }





























get_number_and_percentage_of_directionality_of_exp_and_comp_first_column_names <- function(spread_df_, name_of_df_, gene_col = opts[['ext_gene_name_col']], dir_ = opts$data[['folder']], inside_dir_name = 'exp_and_comp_numbers', exp_col = opts$check_cols_for_input[['Exp.']], logfc_col = opts$check_cols_for_input[['logFC']], pattern_to_extract_exp_from_exp_and_comp_ = '_.*')
{
  dir.create(paste0(dir_, '/', inside_dir_name))
  
  validate_col_types(
    df_ = spread_df_, 
    col_names_list = list(gene_col), 
    col_types_list = list('character'))
  
  test_array <- as.matrix(spread_df_[, colnames(spread_df_) %nin% gene_col]) # Check if all other columns beside first are numeric
  
  spread_df_[is.na(spread_df_)] <- 0
  
  
  
  comparison <- get_number_and_percentage_of_directionality_of_comp_first_column_names(
    spread_dataset_ = spread_df_, 
    gene_col_ = gene_col)
  
  logic_comparison <- comparison$no_of_comps != 0 ### !!! how is this even possible? debug
  
  comparison <- subset(x = comparison, subset = logic_comparison)

  
  
  experiment <- get_number_and_percentage_of_directionality_of_exp_first_column_names(spread_data_ = spread_df_, gene_col_ = gene_col, exp_col_ = exp_col, logfc_col_ = logfc_col, pattern_to_extract_exp_from_exp_and_comp = pattern_to_extract_exp_from_exp_and_comp_)
  
  exp_and_comp <- merge(comparison, experiment, by = gene_col, all = T)
  
  # write.table(
  #   x = exp_and_comp, 
  #   file = paste0(dir_, '/', inside_dir_name, '/exp_and_comp_nb_and_perc_', name_of_df_, '.tsv'), 
  #   sep = '\t', 
  #   dec = ',', 
  #   row.names = F)
  
  
  
  return(exp_and_comp)
}


































get_number_and_percentage_of_directionality_of_comp_first_column_names <- function(spread_dataset_, gene_col_)
{
  pseudomatrix <- as.data.frame(spread_dataset_[, colnames(spread_dataset_) %nin% gene_col_])

  no_of_comps <- purrrlyr::by_row(
    .d = pseudomatrix, 
    .collate = "rows", 
    ..f = function(x) {
      sum(x != 0)}
  )[['.out']]
  
  temp_2 <- purrrlyr::by_row(
    .d = pseudomatrix, 
    .collate = "rows", 
    ..f = function(x) {
      sum(x > 0)}
  )[['.out']]
  
  perc_of_upregulated <- (temp_2/no_of_comps) * 100
  
  spread_dataset_ <- cbind(
    spread_dataset_[gene_col_], 
    no_of_comps, 
    perc_of_upregulated)
  
  return(spread_dataset_)
}



























get_number_and_percentage_of_directionality_of_exp_first_column_names <- function(spread_data_, gene_col_, exp_col_, logfc_col_, pattern_to_extract_exp_from_exp_and_comp)
{
  spread_data_[spread_data_ == 0] <- NA

  temp_ <- tidyr::gather(
    data = spread_data_, 
    key = !!exp_col_, 
    value = !!logfc_col_, 
    na.rm = T, 
    -gene_col_)

  temp_ <- temp_[,-which(colnames(temp_) == logfc_col_)]
  
  temp_[[exp_col_]] <- stringr::str_remove(
    string = temp_[[exp_col_]], 
    pattern = pattern_to_extract_exp_from_exp_and_comp)

  temp_ <- temp_ %>%
    unique() %>%
    dplyr::group_by(eval(parse(text = gene_col_))) %>%
    dplyr::summarise(number_of_exp = dplyr::n())
  
  colnames(temp_)[[1]] <- gene_col_

  return(temp_)
}






























get_number_of_genes_detected_in_given_number_of_experiments <- function(dataset_, number_of_times_detected_col, cols_to_remove)
{
  dataset_ <- dataset_ %>%
    dplyr::group_by(eval(parse(text = number_of_times_detected_col))) %>%
    dplyr::mutate(nb_of_genes_detected_in_this_nb_of_papers = dplyr::n()) %>%
    dplyr::mutate(percent_of_genes_detected_in_this_nb_of_papers = (dplyr::n()/length(dataset_[[1]])) )
  
  dataset_[['eval(parse(text = number_of_times_detected_col))']] <- NULL
  
  dataset_ <- subset(
    x = dataset_, 
    select = colnames(dataset_) %nin% cols_to_remove)
  
  dataset_ <- unique(dataset_)
  
  dataset_ <- dataset_[order(dataset_[[number_of_times_detected_col]]),]
  
  return(dataset_)
}


































collapse_highest_numbers_of_experiments_for_plotting_nb_of_genes_detected_in_given_nb_of_exps <- function(data_, nb_of_exp_to_collapse_higher_data_into, number_of_times_detected_col_, nb_of_genes_detected_in_this_nb_of_papers_col = 'nb_of_genes_detected_in_this_nb_of_papers', percent_of_genes_detected_in_this_nb_of_papers_col = 'percent_of_genes_detected_in_this_nb_of_papers')
{
  assertthat::assert_that(any(data_[[number_of_times_detected_col_]] == nb_of_exp_to_collapse_higher_data_into), msg = 'nb_of_exp_to_collapse_higher_data_into value is not in data')
  
  data_ <- as.data.frame(data_)
  
  data_ <- data_[order(data_[[number_of_times_detected_col_]]),]
  
  data_$rownb <- seq_along(data_[[1]])
  
  rownames(data_) <- data_$rownb
  
  collapse_entry <- data_$rownb[data_[[number_of_times_detected_col_]] == nb_of_exp_to_collapse_higher_data_into]
  
  final_entry <- max(data_$rownb)
  
  data_[[nb_of_genes_detected_in_this_nb_of_papers_col]][collapse_entry] <- sum(
    data_[[nb_of_genes_detected_in_this_nb_of_papers_col]][collapse_entry:final_entry])
  
  data_[[percent_of_genes_detected_in_this_nb_of_papers_col]][collapse_entry] <- sum(data_[[percent_of_genes_detected_in_this_nb_of_papers_col]][collapse_entry:final_entry]) # probably dont need this at all, could order it through ggplot
  
  data_ <- data_[-c((collapse_entry + 1):final_entry),]
  
  data_$rownb <- NULL

  return(data_)
}































get_percentage_of_our_data_are_protein_coding_genes <- function(species, spread_gene_names_col, biomart_protein_database = 'ensembl_peptide_id')
{
  normalized_species <- normalize_species_names(species__vec = species)
  
  biomaRt_ <- set_mart_to_be_used(species_ = normalized_species)
  
  bm_results_all <- biomaRt::getBM(
    attributes = c(biomart_protein_database, "external_gene_name"),
    uniqueRows = T,
    mart = biomaRt_)
  
  logic_bm_results_subset <- bm_results_all[[biomart_protein_database]] != ''
  
  bm_results_subset <- subset(x = bm_results_all, subset = logic_bm_results_subset)
  
  bm_results_subset_all <- data.frame(
    'external_gene_name' = tolower(unique(bm_results_subset[['external_gene_name']])), 
    'protein' = T)
  
  our_genes <- data.frame(
    'external_gene_name' = spread_gene_names_col, 
    'in_our_data' = T)
  
  mergement <- merge(x = our_genes, y = bm_results_subset_all, by = 'external_gene_name')
  
  return(length(mergement[[1]])/length(spread_gene_names_col))
}
























map_data_from_list_of_gene_groups <- function(gene_groups, gather_dataset, gene_name_col = opts[['ext_gene_name_col']], columns_to_analyze_charvec, element_name_with_gene_names_within_gene_group_element, descriptions, exp_and_comp_id_for_data_and_descriptions = opts$data[['exp_and_comp_col']])
{
  mapped_data <- purrr::map2(
    .x = gene_groups, 
    .y = names(gene_groups), 
    .f = function(group, group_name){

      logic_group_expression <- gather_dataset[[gene_name_col]] %in% group[[element_name_with_gene_names_within_gene_group_element]]
      
      group_expression <- subset(x = gather_dataset, subset = logic_group_expression)
      
      categories <- purrr::map(
        .x = columns_to_analyze_charvec,
        .f = function(column)
          {
          temp <- group_expression

          temp[[exp_and_comp_id_for_data_and_descriptions]] <- recode_values_based_on_key(
            to_recode_chrvec = temp[[exp_and_comp_id_for_data_and_descriptions]], 
            replace_this_chrvec = descriptions[[exp_and_comp_id_for_data_and_descriptions]], 
            with_this_chrvec = descriptions[[column]])
          
          logic_temp <- !is.na(temp[[exp_and_comp_id_for_data_and_descriptions]])
          
          temp <- subset(x = temp, subset = logic_temp)
          
          return(temp)
        })
      
      names(categories) <- columns_to_analyze_charvec

      return(categories)
    })
  
  names(mapped_data) <- names(gene_groups)  
  
  return(mapped_data)
}



















print_figures_for_mapped_groups_wrapper <- function(map_data_from_list_of_gene_groups_output, logFC_col = opts$data[['logfc_median_col']], data_folder = opts$data[['folder']], gene_group_folder = opts$data[['folder_gene_groups']], exp_and_comp_col_ = opts$data[['exp_and_comp_col']], gene_name_col = opts[['ext_gene_name_col']])
{
  library(ggplot2)
  
  dir.create(paste0(
    data_folder, 
    '/', 
    gene_group_folder,
    '/figures'))
  
  purrr::walk2(
    .x = map_data_from_list_of_gene_groups_output,
    .y = names(map_data_from_list_of_gene_groups_output),
    .f = function(gene_group, gene_group_name){
      
      dir.create(paste0(
        data_folder, 
        '/', 
        gene_group_folder,
        '/figures/',
        gene_group_name))
      
      purrr::walk2(
        .x = gene_group,
        .y = names(gene_group),
        .f = function(column, column_name){
          
          nb_of_genes <- length(unique(column[[gene_name_col]]))
          
          if (nb_of_genes < 7){
            column %>%
              ggplot(aes(
                x = eval(parse(text = exp_and_comp_col_)),
                y = eval(parse(text = logFC_col)),
                colour = eval(parse(text = gene_name_col)),
                shape = eval(parse(text = gene_name_col)) )) +
              geom_jitter(
                width = 0.2,
                height = 0) +
              labs(title = tolower(column_name)) +
              xlab(NULL) +
              ylab(logFC_col) +
              scale_shape_discrete(name = 'gene name') +
              scale_colour_discrete(name = 'gene name') +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            
          } else {
            column %>%
              ggplot(aes(
                x = eval(parse(text = exp_and_comp_col_)),
                y = eval(parse(text = logFC_col)) )) +
              geom_jitter(
                width = 0.2,
                height = 0,
                size = 2/log2(nb_of_genes + 1)) +
              labs(title = tolower(column_name)) +
              xlab(NULL) +
              ylab(logFC_col) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
          }

          ggsave(paste0(
            data_folder,
            '/',
            gene_group_folder,
            '/figures/',
            gene_group_name,
            '/',
            column_name,
            '.jpg'))
        })
    })
}




















print_tables_and_figures_for_preety_clusters_wrapper <- function(spread_matrixi_list, cluster_name, folder = opts$data[['folder']], subfolder = opts$data[['folder_clusters']], additional_info_for_file, sep_ = '\t', dec_ = '.', comparisons_to_print)
{
  subdir <- paste0(folder, '/', subfolder)
  dir.create(folder)
  dir.create(subdir)
  
  to_print <- spread_matrixi_list[names(spread_matrixi_list) %in% comparisons_to_print]
  
  purrr::walk2(
    .x = to_print, 
    .y = names(to_print), 
    .f = function(spread_matrix, spread_matrix_name)
      {
      assertthat::assert_that(is.matrix(spread_matrix), msg = paste0(
        'cluster ', 
        cluster_name, 
        ' needs to be in matrix format'))
      
      filename <- paste0(
        subdir,
        '/',
        cluster_name, 
        '_',
        spread_matrix_name,
        '_preety_cluster.')
      
      if (
        length(as.data.frame(spread_matrix)[[1]]) > 1 & 
        length(colnames(as.data.frame(spread_matrix))) > 1) { ### !!! ncol?
        
        tiff(
          paste0(filename, 'tiff'), 
          width = 1920,  
          height = 1080)
        
        gplots::heatmap.2(
          x = spread_matrix, 
          trace = "none", 
          dendrogram = 'column', 
          col = colorRamps::matlab.like,
          lwid = c(0.2,4),
          margins = c(6,8)
        )
        
        dev.off()
      }
      
      write.table(
        x = as.data.frame(spread_matrix), 
        file = paste0(filename, 'tsv'),
        sep = sep_,
        row.names = T,
        col.names = NA,
        dec = dec_) # Table for clustering
    })
}

















prepare_spread_data_for_getting_preety_clusters <- function(genesets_to_analyze_group, gene_name_col_in_geneset_to_analyze_group, logfc_cutoff, ratio_of_0_values_cutoff, spread_datasets, ext_gene_name_col_ = opts[['ext_gene_name_col']])
{
  
  spread_groups <- purrr::map2(
    .x = spread_datasets, 
    .y = names(spread_datasets), 
    .f = function(spread_group, spread_group_name){
      
      # if (spread_group_name == "Drug_ketamine") {
      #   browser()
      # }
      
      logic_spread_group <- spread_group[[ ext_gene_name_col_ ]] %in% genesets_to_analyze_group[[gene_name_col_in_geneset_to_analyze_group]] ### !!! does 'genes' column surely include the updated genes, just as they are in dataset_annotated element of group?
      
      spread_group <- as.data.frame(subset(x = spread_group, subset = logic_spread_group))
      
      if (bazar::is.empty(spread_group)) {
        return(NA)
      }
  
      rownames(spread_group) <- spread_group[[ ext_gene_name_col_ ]]
      
      spread_group[[ ext_gene_name_col_ ]] <- NULL
      
      spread_group[is.na(spread_group)] <- 0
      
      cols_to_kill <- as.numeric()
  
      for (col_nb in seq_along(colnames(spread_group))) {
  
        spread_group[[col_nb]][ abs(spread_group[col_nb]) < logfc_cutoff ] <- 0
        
        number_of_0 <- length(spread_group[[col_nb]][ spread_group[col_nb] == 0 ])
        
        number_of_all <- length(spread_group[[col_nb]])
        
        assertthat::assert_that(ratio_of_0_values_cutoff >= 0 & ratio_of_0_values_cutoff <= 1, msg = 'ratio_of_0_values_cutoff needs to be between 0 and 1')
        

        
        if ( (number_of_0/number_of_all) > ratio_of_0_values_cutoff ) {
          cols_to_kill <- c(cols_to_kill, col_nb)
        }
      }
      
      if (!bazar::is.empty(cols_to_kill)) {
        spread_group <- subset(spread_group, select = -cols_to_kill)
      }

      if (bazar::is.empty(spread_group) | all(colnames(spread_group) == '0')) {
        return(NA)
      }
      
      spread_group <- as.matrix(spread_group)
      
      return(spread_group)
  })
  
  return(spread_groups)

}





















extract_genesets_from_datasets <- function(groups_to_analyze)
{
  return_ <- purrr::map2(
    .x = groups_to_analyze,
    .y = names(groups_to_analyze),
    .f = function(group, group_name){
      
      if (!is.null(group[['dataset']])) 
        {
        assertthat::assert_that(
          length(group[['dataset']][[ group[['gene_name_col']] ]]) == length(unique(group[['dataset']][[ group[['gene_name_col']] ]])), 
          msg = paste0('gene names are not unique in dataset ', group_name))
        
        group[['genes']] <- group[['dataset']][[ group[['gene_name_col']] ]]
        
        return(group)
        
      } else return(group)
    })
  
  return(return_)
}




















ncbi_annotation_for_gene_vectors <- function(gene_name_vec, species, str_sep_ = opts$str_sep, Probe_col = opts$check_cols_for_input[['Probe']], external_gene_name_ = opts[['ext_gene_name_col']], columns = list(
  'ENSG_col' = opts$check_cols_for_input$ENSG_, 
  'ENST_col' = opts$check_cols_for_input$ENST_, 
  'Gene_ID_col' = opts$check_cols_for_input$Gene_ID, 
  'NM_col' = opts$check_cols_for_input$NM_, 
  'Accession_col' = opts$check_cols_for_input$Accession, 
  'Unigene_col' = opts$check_cols_for_input$Unigene, 
  'NR_col' = opts$check_cols_for_input$NR_, 
  'XM_col' = opts$check_cols_for_input$XM_, 
  'XR_col' = opts$check_cols_for_input$XR_))
{
  dummy_description <- data.frame(
    'des_species_col' = species,
    'exp_id_col' = 1,
    stringsAsFactors = F)

  dummy_input <- data.frame(
    'input_id_col' = gene_name_vec, 
    'exp_id_col' = 1,
    stringsAsFactors = F)
  
  dummy_input[[Probe_col]] <- as.character(NA)

  

  message('Translating all gene names in the data to gene ids using ncbi... ')
  ncbi_annotation_of_symbols_to_gene_ids <- master_annotator(
    descriptions_df = dummy_description,
    exp_id_col = 'exp_id_col',
    des_species_col = 'des_species_col',
    input_df = dummy_input,
    input_id_col = 'input_id_col',
    C_legacy_D_str_identifier_type__ = 'Gene name',
    PERFORM_D_ncbi_annotation = T, 
    A_D_C_string_separator__ = str_sep_)
  message('... translated.')
  
  
  
  message('Adding ncbi-derived gene ids to original data... ')
  input_plus_new_gene_id_from_ncbi_column <- master_annotator(
    descriptions_df = dummy_description,
    exp_id_col = 'exp_id_col',
    des_species_col = 'des_species_col',
    input_df = dummy_input,
    input_id_col = 'input_id_col',
    PERFORM_E_add_new_gene_id_col_originating_from_ncbi = T, 
    E_PERFORM_D_output = ncbi_annotation_of_symbols_to_gene_ids,
    A_D_C_string_separator__ = str_sep_)
  message('.. added.')
  
  
  
  message('Creating database for annotation...')
  pseudomemoize_external_gene_names_to_gene_id <- master_annotator(
    descriptions_df = ncbi_annotation_of_symbols_to_gene_ids,
    exp_id_col = 'dummy_exp',
    des_species_col = 'organism',
    input_df = ncbi_annotation_of_symbols_to_gene_ids,
    input_id_col = columns[['Gene_ID_col']],
    A_C_all_Gene_ID_col = columns[['Gene_ID_col']],
    PERFORM_A_should_i_prepare_dbs_for_pseudomemoization = T,
    A_return_qa_of_pseudomemoization = F,
    A_D_C_string_separator__ = str_sep_)
  message('... created.')
  
  
  
  message('Starting annotation...')
  result_list <- master_annotator(
    descriptions_df = dummy_description,
    exp_id_col = 'exp_id_col',
    des_species_col = 'des_species_col',
    input_df = input_plus_new_gene_id_from_ncbi_column,
    input_id_col = Probe_col,
    PERFORM_C_annotation = T,
    C_legacy_D_str_identifier_type__ = columns[['Gene_ID_col']],
    A_C_all_Gene_ID_col = 'Gene_ID_from_ncbi',
    C_pseudo_memoized_db = pseudomemoize_external_gene_names_to_gene_id,
    A_D_C_string_separator__ = str_sep_)
  message('... annotated.')
  
  result_df <- result_list[[1]]

  result_df <- subset(
  x = result_df, 
  select = c('input_id_col', external_gene_name_, 'identifed_identifer'))
  
  for (row_nb in seq_along(result_df[[1]])) {
    
    if (is.na(result_df[[external_gene_name_]][[row_nb]])) 
      {
      result_df[[external_gene_name_]][[row_nb]] <- result_df[['input_id_col']][[row_nb]]
    }
  }
  
  result_df[[external_gene_name_]] <- tolower(result_df[[external_gene_name_]])

  return(result_df)
}














merge_group_and_numbers_wrapper <- function(group_, analysis_name, ext_gene_name_col = opts[['ext_gene_name_col']], dataset_to_merge_group_with, cols_to_keep_c = c('no_of_comps', 'perc_of_upregulated', 'number_of_exp'), dir_ =  opts$data[['folder']], inside_dir_name = 'exp_and_comp_numbers_groups')
{
  group_and_all_perc_data <- merge(
    x = data.frame('gene_name' = group_[['genes']]), 
    y = dataset_to_merge_group_with, 
    by.x = 'gene_name',
    by.y = ext_gene_name_col,
    all.x = T)
  
  colnames(group_and_all_perc_data)[colnames(group_and_all_perc_data) == 'gene_name'] <- ext_gene_name_col


  group_and_all_perc_data <- subset(
    group_and_all_perc_data, 
    select = c(ext_gene_name_col, cols_to_keep_c))
  
  dir.create(paste0(dir_, '/', inside_dir_name))
  
  write.table(
    x = group_and_all_perc_data, 
    file = paste0(dir_, '/', inside_dir_name, '/exp_and_comp_nb_and_perc_', analysis_name, '.tsv'), 
    sep = '\t', 
    dec = ',', 
    row.names = F)
  

return(group_and_all_perc_data)
}























get_specific_subset_of_genes_from_spread_data_using_descriptions <- function(spread_data, descriptions_, col_with_covariate_in_desc, covariates_name, gene_name_col, exp_name_col)
{
  exp_with_covariate <- subset(descriptions_, subset = normalize_species_names(descriptions_[[col_with_covariate_in_desc]]) %in% covariates_name)[[exp_name_col]]
  
  spread_data_cov <- subset(spread_data, select = colnames(spread_data) %in% c(exp_with_covariate, gene_name_col))
  
  return(spread_data_cov[[gene_name_col]])
}





















are_genes_in_group_expressed_in_the_same_direction <- function(spread_data_list_from_genegroup, name_, ext_gene_name_col = opts[['ext_gene_name_col']], logFC_col_ = opts$check_cols_for_input[['logFC']], exp_and_comp_col_= opts$data[['exp_and_comp_col']], dir = opts$data[['folder']], inner_dir = 'stability')
{
  purrr::map2(.x = spread_data_list_from_genegroup, .y = names(spread_data_list_from_genegroup), .f = function(spread_data_from_genegroup, spread_data_from_genegroup_name)
  {
    temp <- as.data.frame(t(spread_data_from_genegroup))
    
    temp_names <- as.data.frame(rownames(temp), stringsAsFactors = F)
    
    temp_2 <- cbind(temp_names, temp)
    
    colnames(temp_2) <- stringr::str_replace(
      string = colnames(temp_2), 
      pattern = 'rownames\\(temp\\)', 
      replacement = ext_gene_name_col)
    
    temp_3 <- get_number_and_percentage_of_directionality_of_exp_and_comp_first_column_names(
      spread_df_ = temp_2, 
      gene_col = ext_gene_name_col, 
      name_of_df_ = name_, 
      dir_ = dir, 
      inside_dir_name = inner_dir, 
      exp_col = exp_and_comp_col_, 
      logfc_col = logFC_col_)
    
    mean_ <- mean(temp_3[['perc_of_upregulated']])
    
    sd_ <- sd(temp_3[['perc_of_upregulated']])
    
    return_ <- list('mean' = mean_, 'sd' = sd_)
    
    temp_3[['number_of_exp']] <- NULL
    
    colnames(temp_3)[which(colnames(temp_3) == 'no_of_comps')] <- 'no_of_genes_from_the_group'
    
    return_[[spread_data_from_genegroup_name]] <- temp_3
    
    return(return_)
  })
  

}






















get_mean_and_sd_of_perc_of_upregulation_for_all_comparisons <- function(group_with_numbers_, percentage_analyses, ext_gene_name_col_ = opts[['ext_gene_name_col']], perc_of_upregulated_col = 'perc_of_upregulated')
{
  all_comparisons_for_this_gene_group <- purrr::map2(
    .x = percentage_analyses, 
    .y = names(percentage_analyses), 
    .f = function(comparison, comparison_name){
      
      colnames(comparison)[colnames(comparison) == perc_of_upregulated_col] <- paste0(perc_of_upregulated_col, '_comp')
      
      group_plus_comp <- merge(
        x = group_with_numbers_, 
        y = comparison, 
        by = ext_gene_name_col_, 
        all.x = T) 
      
      mean_perc_of_upregulated <- mean(group_plus_comp[[ paste0(perc_of_upregulated_col, '_comp') ]])
      sd_perc_of_upregulated <- sd(group_plus_comp[[ paste0(perc_of_upregulated_col, '_comp') ]])
      
      return(paste0('Mean: ', mean_perc_of_upregulated, ', sd: ', sd_perc_of_upregulated))
    })
}






















get_individual_genes_in_group_from_all_comps_gathered <- function(groups_, full_work_dataset, descriptions_ = data$descriptions, ext_gene_name_col = opts[['ext_gene_name_col']], genes_list_name_in_group_dataset = 'genes', exp_and_comp_col_ = opts$data[['exp_and_comp_col']])
{
  full_work_dataset <- merge(x = full_work_dataset, y = descriptions_, by = exp_and_comp_col_, all.x = T)
  
  purrr::map(
    .x = groups_,
    .f = function(group){

      individual_groups <- purrr::map(
        .x = group[[genes_list_name_in_group_dataset]], 
        .f = function(gene){
          logic_ <- full_work_dataset[[ext_gene_name_col]] == gene
          
          gene_subset <- subset(x = full_work_dataset, subset = logic_)
          
          return(gene_subset)
        })
      
      names(individual_groups) <- group[[genes_list_name_in_group_dataset]]
      
      return(individual_groups)
    })
}























is_characteristic_enriched_wrapper <- function(gene_group, gene_group_names, full_dataset, column_of_interest, is_this_value_in_column_of_interest_enriched, descriptions_, exp_and_comp_col_ = opts$data[['exp_and_comp_col']])
{
  logic_only_comps_present_in_actual_dataset <- descriptions_[[exp_and_comp_col_]] %in% unique(full_dataset[[exp_and_comp_col_]])
  
  only_comps_present_in_actual_dataset <- subset(descriptions_, subset = logic_only_comps_present_in_actual_dataset)
  
  number_of_all_comparisons <- length(unique(only_comps_present_in_actual_dataset[[exp_and_comp_col_]]))
  
  
  
  logic_only_comps_with_values_of_interest <- only_comps_present_in_actual_dataset[[column_of_interest]] == is_this_value_in_column_of_interest_enriched
  
  only_comps_with_values_of_interest <- subset(only_comps_present_in_actual_dataset, subset = logic_only_comps_with_values_of_interest)
  
  number_of_comparisons_with_values_of_interest <- length(unique(only_comps_with_values_of_interest[[exp_and_comp_col_]]))
  
  
  
  number_of_comparisons_with_values_of_interest_to_of_all_comparisons <- number_of_comparisons_with_values_of_interest/number_of_all_comparisons
  
  
  
  number_of_values_of_interest_for_genes_in_group <- purrr::map(
    .x = gene_group,
    .f = function(gene_table){
      
      nb_of_comps_in_which_the_gene_was_detected <- length(unique(gene_table[[exp_and_comp_col_]]))
      
      
      
      logic_gene_group_with_values_of_interest <- gene_table[[column_of_interest]] == is_this_value_in_column_of_interest_enriched
      
      gene_group_with_values_of_interest <- subset(gene_table, subset = logic_gene_group_with_values_of_interest)
      
      nb_of_comps_with_value_of_interest_in_which_the_gene_was_detected <- length(unique(gene_group_with_values_of_interest[[exp_and_comp_col_]]))

      nb_of_comps_with_value_of_interest_in_which_the_gene_was_detected_to_nb_of_comps_in_which_the_gene_was_detected <- nb_of_comps_with_value_of_interest_in_which_the_gene_was_detected/nb_of_comps_in_which_the_gene_was_detected 

      return(nb_of_comps_with_value_of_interest_in_which_the_gene_was_detected_to_nb_of_comps_in_which_the_gene_was_detected)
    })
  
  names(number_of_values_of_interest_for_genes_in_group) <- gene_group_names

  
  
  mean_nb_of_comps_with_value_of_interest_in_which_genes_were_detected_to_nb_of_comps_in_which_genes_were_detected <- mean(na.omit(as.numeric(number_of_values_of_interest_for_genes_in_group)))
  
  sd_nb_of_comps_with_value_of_interest_in_which_genes_were_detected_to_nb_of_comps_in_which_genes_were_detected <- sd(na.omit(as.numeric(number_of_values_of_interest_for_genes_in_group)))

  
  
  return(list(
    'number_of_comparisons_with_values_of_interest_to_of_all_comparisons' = number_of_comparisons_with_values_of_interest_to_of_all_comparisons,
    'nb_of_comps_with_value_of_interest_in_which_the_gene_was_detected_to_nb_of_comps_in_which_the_gene_was_detected' = number_of_values_of_interest_for_genes_in_group, 
    'mean_nb_of_comps_with_value_of_interest_in_which_genes_were_detected_to_nb_of_comps_in_which_genes_were_detected' = mean_nb_of_comps_with_value_of_interest_in_which_genes_were_detected_to_nb_of_comps_in_which_genes_were_detected, 
    'sd_nb_of_comps_with_value_of_interest_in_which_genes_were_detected_to_nb_of_comps_in_which_genes_were_detected' = sd_nb_of_comps_with_value_of_interest_in_which_genes_were_detected_to_nb_of_comps_in_which_genes_were_detected)) 
}




















get_enrichment_for_values_in_columns_for_groups_wrapper <- function(named_list_of_cols_and_values_to_analyze_, full_dataset_, individual_genes_in_group_from_all_comps_gathered_, descriptions__ = data$descriptions)
{
  purrr::map(
    .x = named_list_of_cols_and_values_to_analyze_,
    .f = function(cols_of_interest_){
      
      data_on_column <- purrr::map(
        .x = cols_of_interest_[['values_of_interest']], 
        .f = function(value) {
          
          data_on_value_within_column <- purrr::map(
            .x = individual_genes_in_group_from_all_comps_gathered_,
            .f = function(genes_in_group){
              
              genes_ <- is_characteristic_enriched_wrapper(
                gene_group = genes_in_group, ### !!! here goes entire gene group, and not individual genes
                gene_group_names = names(genes_in_group),
                full_dataset = full_dataset_,
                is_this_value_in_column_of_interest_enriched = value, 
                column_of_interest = cols_of_interest_[['col_name']],
                descriptions_ = descriptions__) 

              return(genes_)
            })

          return(data_on_value_within_column)
        })
      
      names(data_on_column) <- cols_of_interest_[['values_of_interest']]
      return(data_on_column)
    })
}

















get_number_of_times_the_gene_detected_in_paper <- function(individual_genes_in_group_from_all_comps_gathered_, exp_and_comp_col_ = opts$data[['exp_and_comp_col']], logfc_col = opts$data[['logfc_median_col']], pattern_to_extract_exp_from_exp_and_comp)
{
  purrr::map(
    .x = individual_genes_in_group_from_all_comps_gathered_,
    .f = function(gene_group){
      
      genes_ <- purrr::map(
        .x = gene_group,
        .f = function(gene_) {
          
          logic_temp <- gene_[[logfc_col]] != 0
          
          temp_ <- subset(x = gene_, subset = logic_temp)
          
          experiment <- stringr::str_remove(temp_[[exp_and_comp_col_]], pattern = pattern_to_extract_exp_from_exp_and_comp) 
          
          result_df <- as.data.frame(table(experiment))
          
          result_df[['experiment']] <- as.character(result_df[['experiment']])
          
          return(result_df)
        })
      names(genes_) <- names(gene_group)
      
      return(genes_)
    })
}




















get_number_of_exps_for_paper_per_gene_wrapper <- function(how_many_times_were_genes_detected_in_paper, experiment_to_check_number_of_comparisons_in, experiment_name_col_name = 'experiment', frequency_col_name = 'Freq')
{
  numbers_how_many_times_gene_detected_in_exp_per_paper <- list()

  for (experiment in experiment_to_check_number_of_comparisons_in) {
    
    numbers_how_many_times_gene_detected_in_exp_per_paper[[experiment]]$numbers <-
      purrr::map(
        .x = how_many_times_were_genes_detected_in_paper,
        .f = function(gene) {

          logic_how_many_times_was_the_gene_detected_in_paper <- as.character(gene[[experiment_name_col_name]]) == experiment

          temp <- subset(gene, subset = logic_how_many_times_was_the_gene_detected_in_paper)

          if (length(temp[[1]]) == 1) {
            
            return( temp[[frequency_col_name]] )
            
          } else if (length(temp[[1]]) == 0) {
            
            return(0)
            
          } else stop('ERROR: duplicated rows in numbers_how_many_times_gene_detected_in_exp_per_paper')
        })

    
    numbers_how_many_times_gene_detected_in_exp_per_paper[[experiment]]$mean <- mean(as.numeric(numbers_how_many_times_gene_detected_in_exp_per_paper[[experiment]][['numbers']]))
    
    numbers_how_many_times_gene_detected_in_exp_per_paper[[experiment]]$sd <- sd(as.numeric(numbers_how_many_times_gene_detected_in_exp_per_paper[[experiment]][['numbers']]))
  }

  return(numbers_how_many_times_gene_detected_in_exp_per_paper)
}


















get_exps_in_which_gene_was_present_for_given_value_of_interest_wrapper <- function(gene_group_, logFC_col = opts$data[['logfc_median_col']], exp_col = opts$data[['exp_and_comp_col']], gene_name_col = opts[['ext_gene_name_col']], descriptions_ = data$descriptions, list_for_subsetting_name_is_colname_value_is_value_of_interest)
{
  gene_non_zero_df <- purrr::map2(
    .x = gene_group_,
    .y = names(gene_group_),
    .f = function(gene, gene_name)
      {
      logic <- gene[[logFC_col]] != 0

      gene_subset <- subset(gene, subset = logic)
      
      exp_col_ <- gene_subset[[exp_col]]
      
      return_ <- data.frame(
        'experiment' = as.character(exp_col_), 
        'gene_name' = as.character(rep(gene_name, length(exp_col_))) )

      return(return_)
    })

  
  gene_non_zero_df_bind <- rlist::list.rbind(gene_non_zero_df)
  
  gene_non_zero_df_bind$present <- T
  
  gene_non_zero_df_bind_spread <- tidyr::spread(
    gene_non_zero_df_bind, 
    key = gene_name,
    value = present)
  
  gene_non_zero_df_bind_spread$end_of_genes <- NA
  
  gene_non_zero_df_bind_spread_merge <- merge(
    x = gene_non_zero_df_bind_spread, 
    y = descriptions_, 
    by.x = 'experiment', 
    by.y = exp_col, 
    all.x = T)

  
  
  gene_non_zero_df_bind_spread_merge_for_values_of_interest <- purrr::map(
    .x = list_for_subsetting_name_is_colname_value_is_value_of_interest, 
    .f = function(values_of_interest)
      {
      logic_dfs_ <- gene_non_zero_df_bind_spread_merge[[ values_of_interest[['col_name']] ]] %in%
        values_of_interest[['values_of_interest']]
      
      dfs_ <- subset(gene_non_zero_df_bind_spread_merge, subset = logic_dfs_)

      cols_with_genes <- dfs_[,(which(colnames(dfs_) == 'experiment') + 1):(which(colnames(dfs_) == 'end_of_genes') - 1)]
      
      nb_of_values_of_intr <- purrr::map(
        .x = cols_with_genes, 
        .f = function(gene) {
          length(subset(gene, gene == T))})
      
      
      
      return(list(
        'dfs_' = dfs_, 
        'col_with_genes' = cols_with_genes, 
        'nb_of_values_of_intr' = nb_of_values_of_intr, 
        'mean_nb_of_values_of_intr' = mean(as.integer(nb_of_values_of_intr)),
        'sd_nb_of_values_of_intr' = sd(as.integer(nb_of_values_of_intr))))
    })

  return(gene_non_zero_df_bind_spread_merge_for_values_of_interest)
}
