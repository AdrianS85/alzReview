master_annotator <- function(
  descriptions_df, 
  des_species_col, 
  input_df,
  input_id_col, # used in: get_the_highest_hit_returning_id_type, annotate_identifiers_to_geneID, annotate_no,
  exp_id_col, 
  PERFORM_A_should_i_prepare_dbs_for_pseudomemoization = F,
  A_return_qa_of_pseudomemoization = F,
  PERFORM_B_get_the_highest_hit_returning_id_type = F,
  B_highest_hit_int_Probe_IDs_to_test = 200,
  PERFORM_C_annotation = F,
  C_pseudo_memoized_db = NULL,
  A_C_all_ENSG_col = NA, A_C_all_ENST_col = NA, A_C_all_Gene_ID_col = NA, A_C_all_NM_col = NA, A_C_all_Accession_col = NA, A_C_all_Unigene_col = NA, A_C_all_NR_col = NA, A_C_all_XM_col = NA, A_C_all_XR_col = NA, 
  PERFORM_D_ncbi_annotation = F,
  A_D_C_string_separator__,
  C_legacy_D_str_identifier_type__, #Generally: 'Gene name' for annotate_identifiers_to_geneID. We can pass two types of data into str_identifier_name: 1) one-element string vector, which is the same name as filter name in biomartr-ensembl database. 2) vector of strings with filter name for each experiment dataset in the list of experiments )
  PERFORM_E_add_new_gene_id_col_originating_from_ncbi = F,
  E_PERFORM_D_output
  )
{
  check_nb_of_options_selected(
    A = PERFORM_A_should_i_prepare_dbs_for_pseudomemoization,
    B = PERFORM_B_get_the_highest_hit_returning_id_type,
    C = PERFORM_C_annotation,
    D = PERFORM_D_ncbi_annotation,
    E = PERFORM_E_add_new_gene_id_col_originating_from_ncbi)
  
  library(magrittr)
  library(Hmisc)
  
  message('master_annotator: validating description and input column types')
  validate_col_types(df_ = descriptions_df, col_names_list = list(exp_id_col, des_species_col), col_types_list = list('numeric', 'character'))
  
  validate_col_types(df_ = input_df, col_names_list = list(exp_id_col, input_id_col), col_types_list = list('numeric', 'character'))

  
  message('master_annotator: ordering descriptions and input by exp_id_col')
  descriptions_df <- descriptions_df[order(descriptions_df[[exp_id_col]]),]
  input_df <- input_df[order(input_df[[exp_id_col]]),]

  
  
  message('master_annotator: subsetting exp_id_col, des_species_col from descriptions')
  exp_species__ <- descriptions_df %>%
    dplyr::select(exp_id_col, des_species_col) %>%
    unique()
  

  ### !!! this seems to be bugged. If I dont first set subset_ value, than for 1 row exp_species subset returns 0 rows. I cannot replicate this behavior outside the function
  subset_ <- exp_species__[[exp_id_col]] %in% unique(input_df[[exp_id_col]])
  exp_species__ <- subset(x = exp_species__, subset = subset_)

  
  message('master_annotator: normalizing species names')
  exp_species__[[des_species_col]] <- normalize_species_names(exp_species__[[des_species_col]])
  
  
  message('master_annotator: splitting input by experiment column')
  LIST_DATA__ <- split(input_df, f = input_df[[exp_id_col]])

  if (PERFORM_A_should_i_prepare_dbs_for_pseudomemoization == T) {
    message('master_annotator: constructing gene name database')
    pseudo_memoization_db <- prepare_db_for_pseudo_memoization(
      LIST_DATA___ = LIST_DATA__, 
      des_species_vec = exp_species__[[des_species_col]],
      ENSG_ = A_C_all_ENSG_col, ENST_ = A_C_all_ENST_col, Gene_ID = A_C_all_Gene_ID_col, NM_ = A_C_all_NM_col, Accession = A_C_all_Accession_col, Unigene = A_C_all_Unigene_col, NR_ = A_C_all_NR_col, XM_ = A_C_all_XM_col, XR_ = A_C_all_XR_col, 
      return_qa_of_pseudomemoization_ = A_return_qa_of_pseudomemoization,
      string_separator__ = A_D_C_string_separator__)
    
    memoized_db <- pseudo_memoize(pseudo_memoization_db_ = pseudo_memoization_db)
    
    return(memoized_db)
    
  } else if (PERFORM_B_get_the_highest_hit_returning_id_type == T) {
    message('master_annotator: identifying best ensembl-based platforms to be used for analysis')
    highest_hits <- get_the_highest_hit_returning_id_type(
      exp_species_ = exp_species__, 
      des_exp_id_col_ = exp_id_col,
      des_species_col_ = des_species_col,
      list_LIST_DATA_ = LIST_DATA__,
      input_probe_col_ = input_id_col,
      int_Probe_IDs_to_test = B_highest_hit_int_Probe_IDs_to_test)
    
    return(highest_hits)
    
  } else if (PERFORM_C_annotation == T) {
    message('master_annotator: annotating')

    annotation <- annotate_now(
        list_LIST_DATA_ = LIST_DATA__,
        str_identifier_type_ = C_legacy_D_str_identifier_type__,
        str_vector_of_species_names__ = exp_species__[[des_species_col]], 
        input_id_col_ = input_id_col,
        ENSG_col = A_C_all_ENSG_col, ENST_col = A_C_all_ENST_col, Gene_ID_col = A_C_all_Gene_ID_col, NM_col = A_C_all_NM_col, Accession_col = A_C_all_Accession_col, Unigene_col = A_C_all_Unigene_col, NR_col = A_C_all_NR_col, XM_col = A_C_all_XM_col, XR_col = A_C_all_XR_col,
        pseudo_memoized_db_ = C_pseudo_memoized_db,
        str_sep_ = A_D_C_string_separator__)
    
    return(annotation)
    
  } else if (PERFORM_D_ncbi_annotation == T) {
    message('master_annotator: translating all gene names in the data to gene ids using ncbi')
    annotated_with_ncbi <- annotate_identifiers_to_geneID(
      descriptions_df = exp_species__, 
      desc_exp_id_col_ = exp_id_col, 
      desc_species_col_ = des_species_col, 
      input_df_ = input_df, 
      input_exp_id_col_ = exp_id_col, 
      input_id_col_ = input_id_col,
      string_separator_ = A_D_C_string_separator__,
      chr_gene_identifier_type = C_legacy_D_str_identifier_type__# should be 'Gene name' for gene name
    )
    
    return(annotated_with_ncbi)
    
  } else if (PERFORM_E_add_new_gene_id_col_originating_from_ncbi == T) {
    message('master_annotator: adding ncbi-derived gene ids to original data')
    df_with_new_col <- add_new_gene_id_col_originating_from_ncbi_annotation(
      descriptions_df = exp_species__, 
      desc_exp_id_col = exp_id_col, 
      organism_col = des_species_col,
      input_df = input_df, 
      input_exp_id_col = exp_id_col,
      input_input_id_col = input_id_col,
      perform_ncbi_annotation_output = E_PERFORM_D_output,
      str_sep = A_D_C_string_separator__)
    
    return(df_with_new_col)
  } else stop('No analysis type was selected. See PERFORM_ arguments for available analyses.')
}




















get_the_highest_hit_returning_id_type <- function(
  exp_species_, 
  des_exp_id_col_,
  des_species_col_,
  list_LIST_DATA_,
  input_probe_col_,
  int_Probe_IDs_to_test = 200)
{
  SHORT_LIST_DATA <- lapply(X = list_LIST_DATA_, FUN = function(x){ x[1:int_Probe_IDs_to_test,] })
  ANNOT_SHORT_LIST_DATA <- rep(list(list()), times = length(SHORT_LIST_DATA))
  all_ID_annotations <- rep(list(list()), times = length(SHORT_LIST_DATA))
  ### Prepare secondary input ###
  
  
  for(n in seq_along(SHORT_LIST_DATA))
  {
    # Here we establish which mart(species) we are using in this given dataset based on "exp_species" vector
    usedMart_ <- set_mart_to_be_used(species_ = exp_species_[[des_species_col_]][n], int_loop = n)

    potental_identifiers <- get_the_potental_identifiers(usedMart___ = usedMart_, int_loop_nb = n)
    
    message('Starting to annotate the data')
    
    for (m in seq_along(potental_identifiers))
    {
      message( paste0('Annotating data for step ', n, ' and ', m, '...'))
      
      ANNOT_SHORT_LIST_DATA[[n]][[m]] <- biomaRt::getBM(
        attributes = c(potental_identifiers[[m]], "external_gene_name"),
        filters = potental_identifiers[[m]],
        values = SHORT_LIST_DATA[[n]][[input_probe_col_]],
        uniqueRows = T,
        mart = usedMart_)
      
      message( paste0('Data for step for ', n, ' and ', m, ' annotated') )
    }
    message( 'Data annotated' )
    
    usedMart_ <- NULL
    
    
    
    message( 'Starting all_ID_annotations step' )
    
    # Here we save number of annotations from each ID
    for(k in seq_along(ANNOT_SHORT_LIST_DATA[[n]]))
    {
      message( paste0('Starting all_ID_annotations step for ', n, ' and ', k) )
      all_ID_annotations[[n]][[k]] <- length(ANNOT_SHORT_LIST_DATA[[n]][[k]][[1]])
      message( paste0('Competed all_ID_annotations step for ', n, ' and ', k) )
    }
    
    message( 'Competed whole all_ID_annotations step' )
  }
  
  names(all_ID_annotations) <- exp_species_[[des_exp_id_col_]]
  
  # Lists inside main list are changed into dfs (vectors) as 'which' function demands it
  df_all_ID_annotations <- lapply(all_ID_annotations, FUN = unlist)
  
  # Here we will be returing results of appropriate microarray search
  HIGHEST_HIT_LIST <- rep(list(list()), times = length(SHORT_LIST_DATA))
  
  # Here we are getting all of the highest yielding IDs
  for(n in seq_along(list_LIST_DATA_)){
    for(m in seq_along(which(df_all_ID_annotations[[n]] == max(df_all_ID_annotations[[n]])))){
      HIGHEST_HIT_LIST[[n]][[m]] <- ANNOT_SHORT_LIST_DATA[[n]][[(which(df_all_ID_annotations[[n]] == max(df_all_ID_annotations[[n]]))[m])]]
    } }
  
  names(HIGHEST_HIT_LIST) <- exp_species_[[des_exp_id_col_]]
  
  rm(df_all_ID_annotations)
  
  
  proper_length_vector_for_checking_annotation_percentage <- get_proper_length_vector_for_checking_annotation_percentage(list_LIST_DATA__ = list_LIST_DATA_, int_Probe_IDs_to_test_ = int_Probe_IDs_to_test) 
  
  check_annotation_percentage <- purrr::map2(
    .x = HIGHEST_HIT_LIST, 
    .y = proper_length_vector_for_checking_annotation_percentage, 
    .f = function(.x, .y) {  
      perc <- (length(unique(.x[[1]][[1]])) / .y) * 100
      
      if (perc < 50) {
        warning(sprintf('\nWARNING!\nWARNING!\nWARNING!\nplatform %s has annotation lower than 50 percent, it may be an artifact', colnames(.x[[1]])[[1]]))
      }
      
      return(perc)
      }
  )
  
  check_annotation_percentage <- data.frame(exp_species_[[des_exp_id_col_]], rlist::list.rbind(check_annotation_percentage))
  colnames(check_annotation_percentage) <- c('Exp_ID', 'highest_annotated_identifier_percentages')
  
  # Here we simply copy/establish names of features, that were highest by themselves. The conditions ask: 1) is there at least a single hit with highest number (I dont know if there can be 0 though...) 2) Is the first (and though each) highest hit list has at least single hit?
  NAMES_HIGHEST_HIT_LIST <- lapply(X = HIGHEST_HIT_LIST, FUN = set_0_hit_annotations_to_na)
  NAMES_HIGHEST_HIT_LIST <- data.frame(exp_species_[[des_exp_id_col_]], rlist::list.rbind(NAMES_HIGHEST_HIT_LIST))
  colnames(NAMES_HIGHEST_HIT_LIST) <- c('Exp_ID', 'platform_to_use')
  
  NAMES_HIGHEST_HIT_LIST <- merge(x = NAMES_HIGHEST_HIT_LIST, y = check_annotation_percentage, by = 'Exp_ID')
  
  NAMES_HIGHEST_HIT_LIST <- merge(x = NAMES_HIGHEST_HIT_LIST, y = exp_species_, by.x = 'Exp_ID', by.y = des_exp_id_col_)
  
  names(ANNOT_SHORT_LIST_DATA) <- exp_species_[[des_exp_id_col_]]

  
  result_list <- list('best_ID' = NAMES_HIGHEST_HIT_LIST, 'best_ID_annotation' = HIGHEST_HIT_LIST, 'annotations_from_every_ID' = ANNOT_SHORT_LIST_DATA, 'perc_annotations_from_every_ID' = all_ID_annotations)

  return(result_list)
}

















set_mart_to_be_used <- function(species_, int_loop = 1)
{
  normalized_species_names <- normalize_species_names(return_normalized_names = T)

  if (species_ %nin% c(normalized_species_names[['mouse']], normalized_species_names[['rat']], normalized_species_names[['human']], normalized_species_names[['sheep']], normalized_species_names[['saimiri']])) {
    stop(paste0('set_mart_to_be_used does not recognize species id, which should be one of those: ', normalized_species_names[['mouse']], normalized_species_names[['rat']], normalized_species_names[['human']], normalized_species_names[['sheep']], normalized_species_names[['saimiri']]))
  }
  
  message( paste0('Setting mart for step ', int_loop, ' with species ', species_, "...") )
  
  done <- F
  counter <- 0
  
  while (done == F) {
    
    tryCatch( {
      if (species_ == normalized_species_names[['mouse']]) {
        usedMart__ <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
      } else if (species_ == normalized_species_names[['rat']]) {
        usedMart__ <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl")
      } else if (species_ == normalized_species_names[['human']]) {
        usedMart__ <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
      } else if (species_ == normalized_species_names[['sheep']]) {
        usedMart__ <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "oaries_gene_ensembl")
      } else if (species_ == normalized_species_names[['saimiri']]) {
        usedMart__ <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "sbboliviensis_gene_ensembl")
      }
      
      done <- T
    },
    error=function(e) {
      Sys.sleep(5)
      warning('SL certificate problem: unable to get local issuer certificate, trying again in 5 s')
    })
  }
  


  message( paste0('Mart set as ', usedMart__@dataset, ' for step ', int_loop) )
  
  return(usedMart__)
}






















get_the_potental_identifiers <- function(usedMart___, int_loop_nb)
{
  message( paste0('Extracting potental_identifiers for step ', int_loop_nb) )
  
  potental_identifiers <- c(usedMart___@filters[grep(pattern = "^ensembl(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^refseq(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^affy(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^agilent(.*)", usedMart___@filters[[1]]) , 1], 
                            usedMart___@filters[grep(pattern = "^illumina(.*)", usedMart___@filters[[1]]) , 1])
  
  message( paste0('Potental identifiers for step ', int_loop_nb, 'extracted') )
  
  return(potental_identifiers)
}
























get_proper_length_vector_for_checking_annotation_percentage <- function(list_LIST_DATA__, int_Probe_IDs_to_test_)
{
  temp <- list()
  for (n in seq_along(list_LIST_DATA__))
  {
    if(length(list_LIST_DATA__[[n]][[1]]) <= int_Probe_IDs_to_test_)
    {
      temp[n] <- length(list_LIST_DATA__[[n]][[1]])
    }
    else
    {
      temp[n] <- int_Probe_IDs_to_test_
    }
  }
  return(rlist::list.rbind(temp))
}





















set_0_hit_annotations_to_na <- function(list_of_dataframes) 
{ 
  if(length(list_of_dataframes) >= 1 && length(list_of_dataframes[[1]][[1]]) > 0 )
  {
    list_of_dataframes <- names(list_of_dataframes[[1]][1])
  }
  else
  {
    list_of_dataframes <- NA
  }
}



















































# INPUT: list. Defaults to writing lenghts of object on first level of depth of the list
write_lenghts_of_list_objects <- function(list_, string_name_of_the_file, int_length_at_this_depth = 1)
{
  temp <- lapply(
    X = list_, 
    FUN = function(x){ 
      length(x[[int_length_at_this_depth]]) })
  temp2 <- as.data.frame( rlist::list.rbind(temp) )
  write.table(temp2, string_name_of_the_file, sep = '\t')
  rm(temp, temp2)
} 





















annotate_now <- function(
  list_LIST_DATA_,
  input_id_col_,
  str_identifier_type_,
  str_vector_of_species_names__,
  ENSG_col, ENST_col, Gene_ID_col, NM_col, Accession_col, Unigene_col, NR_col, XM_col, XR_col,
  pseudo_memoized_db_,
  str_sep_)
{
  if (!is.null(pseudo_memoized_db_)) {
    for (n in seq_along(list_LIST_DATA_)) # loop through experiments
    {
      identifiers_used_for_annotation <- return_all_usable_id_types(
          ENSG_ = ENSG_col,
          ENST_ = ENST_col,
          Gene_ID = Gene_ID_col,
          NM_ = NM_col,
          Accession = Accession_col,
          Unigene = Unigene_col,
          XM_ = XM_col,
          XR_ = XR_col,
          NR_ = NR_col,
          list_LIST_DATA_ = list_LIST_DATA_[[n]])
      
      message('annotate_now: working on step ', n, '; experiment ', names(list_LIST_DATA_)[[n]])
      message(paste0('annotate_now: identifiers used for annotation in this step are: ', identifiers_used_for_annotation))
      message('annotate_now: performing pseudo-memoized database workflow')
      
      if (is.null(list_LIST_DATA_[[n]]$external_gene_name)) {
        list_LIST_DATA_[[n]]$external_gene_name <- NA
      }
      if (is.null(list_LIST_DATA_[[n]]$identifed_identifer)) {
        list_LIST_DATA_[[n]]$identifed_identifer <- NA
      }
      if (is.null(list_LIST_DATA_[[n]]$unidentifed_identifers)) {
        list_LIST_DATA_[[n]]$unidentifed_identifers <- NA
      }
      
      for (gene_row_nb in seq_along(list_LIST_DATA_[[n]][[1]]))
      {
        if (is.na(list_LIST_DATA_[[n]]$external_gene_name[[gene_row_nb]]))
        {
          for (identifiers_of_given_type in identifiers_used_for_annotation)
          {
            
            ### !!! 280521
            pseudo_memoized_db_[[str_vector_of_species_names__[[n]]]][[identifiers_of_given_type$filter_name]][[1]] <- tolower(
              pseudo_memoized_db_[[str_vector_of_species_names__[[n]]]][[identifiers_of_given_type$filter_name]][[1]]) ### !!! 280521
            
            if (is.na(list_LIST_DATA_[[n]]$external_gene_name[[gene_row_nb]]))
            {
              ids_to_check <- stringr::str_split(
                string = list_LIST_DATA_[[n]][[gene_row_nb, identifiers_of_given_type$col_name]],
                pattern = str_sep_,
                simplify = T
              )
              
              for (id in ids_to_check)
              {
                if (is.na(list_LIST_DATA_[[n]]$external_gene_name[[gene_row_nb]])  & !is.na(id))
                {
                  db_to_search <-
                    pseudo_memoized_db_[[str_vector_of_species_names__[[n]]]][[identifiers_of_given_type$filter_name]]
                  
                  getBM_result_logic <- stringr::str_detect(
                    string = db_to_search[[1]],
                    pattern = stringr::fixed(tolower(id))) ### added to lower !!! 280521
                  
                  getBM_result <- subset(
                    x = db_to_search, 
                    subset = getBM_result_logic)
                  
                  if (length(getBM_result[[1]]) == 0) {
                    if (is.na(list_LIST_DATA_[[n]]$unidentifed_identifers[[gene_row_nb]])) {
                      list_LIST_DATA_[[n]]$unidentifed_identifers[[gene_row_nb]] <- id
                      
                    } else {
                      list_LIST_DATA_[[n]]$unidentifed_identifers[[gene_row_nb]] <-
                        paste0(list_LIST_DATA_[[n]]$unidentifed_identifers[[gene_row_nb]],
                               str_sep_,
                               id)
                    }
                    
                  } else {
                    getBM_result_ext <-
                      paste(getBM_result$external_gene_name,
                            collapse = str_sep_)
                    
                    list_LIST_DATA_[[n]]$external_gene_name[[gene_row_nb]] <-
                      getBM_result_ext
                    list_LIST_DATA_[[n]]$identifed_identifer[[gene_row_nb]] <-
                      id
                  }
                } else
                  break
              }
            } else
              break
          }
        } else
          break
      }
    }
    
    return(list_LIST_DATA_)
    
  } else {
    ###################
    ### LEGACY PART ###
    ###################
    ## !! This part doesnt return undientfied_indetifiers (column reporting what was checked, but not found) - for now its ok, but probably to look at in the future
    message('annotate_now: performing legacy workflow')
    
    ANNOT_LIST_DATA <- list()
    
    if (length(str_identifier_type_) != 1)
    {
      identifiers_used_for_annotation <- str_identifier_type_
    } else
    {
      identifiers_used_for_annotation <-
        set_identifiers_used_for_annotation_if_not_probeID(str_identifier_type = str_identifier_type_,
                                                           list_LIST_DATA = list_LIST_DATA_)
    }
    
    for (n in seq_along(list_LIST_DATA_))
    {
      if (!is.na(identifiers_used_for_annotation[[n]])) {
        usedMart_ <- set_mart_to_be_used(species_ = str_vector_of_species_names__[n], int_loop = n)
        
        message( paste0(
            'Annotating experiment ',
            names(list_LIST_DATA_[n]),
            ' in step ',
            n,
            ' with identifier ',
            identifiers_used_for_annotation[[n]]
          ) )

        
        ANNOT_LIST_DATA[[n]] <- biomaRt::getBM(
          attributes = c(
            identifiers_used_for_annotation[[n]],
            "external_gene_name"),
          filters = identifiers_used_for_annotation[[n]],
          values = list_LIST_DATA_[[n]][[input_id_col_]],
          uniqueRows = F,
          mart = usedMart_)

      } else {
        ANNOT_LIST_DATA[[n]] <- tibble::tibble('id_not_found' = list_LIST_DATA_[[n]][[input_id_col_]])
        ANNOT_LIST_DATA[[n]]$external_gene_name <- NA
      }
    }
    
    FINAL_ANNOT_LIST_DATA <- purrr::pmap(
      .l = list(
        list_LIST_DATA_,
        ANNOT_LIST_DATA,
        identifiers_used_for_annotation
      ),
      .f = function(.x, .y, .z)
      {
        if (!is.na(.z)) {
          merge(
            x = .x,
            y = .y,
            by = input_id_col_,
            by.y = .z,
            all.x = T
          )
        } else {
          merge(
            x = .x,
            y = .y,
            by = input_id_col_,
            by.y = 'id_not_found',
            all.x = T
          )
        }
      }
    )

    DF_FINAL_ANNOT_LIST_DATA <-
      rlist::list.rbind(FINAL_ANNOT_LIST_DATA)
    
    DF_FINAL_ANNOT_LIST_DATA <- unique(DF_FINAL_ANNOT_LIST_DATA)
    
    DF_FINAL_ANNOT_LIST_DATA <- subset(DF_FINAL_ANNOT_LIST_DATA, DF_FINAL_ANNOT_LIST_DATA[["external_gene_name"]] != "")
    
    DF_FINAL_ANNOT_LIST_DATA <-
      collapse_annotated_names_for_given_probe(
        DF_FINAL_ANNOT_LIST_DATA,
        list_LIST_DATA__ = list_LIST_DATA_,
        input_id_col__ = input_id_col_,
        str_sep = str_sep_)
    
    DF_FINAL_ANNOT_LIST_DATA[["external_gene_name"]] <- tolower(DF_FINAL_ANNOT_LIST_DATA[["external_gene_name"]]) ### !!! so that all externals were in lower for further comparisons 280421
    
    return(DF_FINAL_ANNOT_LIST_DATA)
  }
  ###################
  ### LEGACY PART ###
  ###################
}

















# If we have one, specific identifer, than we want to apply it to all experiments. If we have a list of identifers, than this means that each of these identifier should be applyied to corresponding experiment
set_identifiers_used_for_annotation_if_not_probeID <- function(str_identifier_type, list_LIST_DATA) 
{
  if(length(str_identifier_type) == 1)
  {
    return( rep(str_identifier_type, length(list_LIST_DATA)) )
  }
  else # This is not good, because if we have single identifier for a list of experiments, but it is not any of abovementioned identifier, than it will fuck up program. the solution should be to ask if the lenght(str_identifier_type) == 1
  {
    return(str_identifier_type)
  }
}




















collapse_annotated_names_for_given_probe <- function(df_annotated, list_LIST_DATA__, input_id_col__, str_sep)
{
  tryCatch({df_annotated <-
             aggregate(as.formula(paste0('external_gene_name ~ ', input_id_col__)),
                       data = df_annotated,
                       FUN = stringr::str_c)},
           error = function(cond) {
             message("there were no rows to aggregate")
             return(df_annotated)
           })
  

  # This part collapses multiple same values to single value and saves it as temp column
  df_annotated$temp_gene_name <-
    as.character(lapply(
      X = df_annotated$external_gene_name,
      FUN = function(x) {
        paste(x, collapse = str_sep)
      }
    ))
  
  df_annotated$external_gene_name <- change_vector_of_mixed_normal_and_c_geneNames_into_unique_geneName(
    df_annotated$temp_gene_name, 
    regex_pattern_to_split_with = str_sep)
  
  df_annotated <- subset(x = df_annotated, select = c(input_id_col__, 'external_gene_name'))
  
  bind_list_LIST_DATA__ <- rlist::list.rbind(list_LIST_DATA__)
  
  merged_df_annotated <-
    merge(x = bind_list_LIST_DATA__,
          y = df_annotated,
          by = input_id_col__,
          all.x = T)
  
  merged_df_annotated <- dplyr::select(merged_df_annotated, colnames(bind_list_LIST_DATA__), dplyr::everything())
  
  merged_df_annotated <- unique(merged_df_annotated)
  
  return(merged_df_annotated)
}


















## !! This is not wrong, but im not sure why is this neccesary
# # Input: Vector of gene names, needs to be a proper char vector, with gene names separ
change_vector_of_mixed_normal_and_c_geneNames_into_unique_geneName <- function(char_vec, regex_pattern_to_split_with)
{
  # This part collapses multiple same values to single value and saves it as temp column
  temp_list <-
    lapply(
      X = char_vec,
      FUN = function(x) {
        temp <- as.character(stringr::str_split(
          string = x,
          pattern = regex_pattern_to_split_with,
          simplify = T
        ))

        temp <- unique(temp)

        return(temp)
      }
    )

  # Check if the value for given row is proper char, or leftover list ' c("blabla") '
  was_the_value_uniqued <-
    as.logical(lapply(
      X = temp_list,
      FUN = function(x) {
        if (length(x) == 1) {
          return(F)
        }
        else{
          return(T)
        }
      }
    ))

  bullshit_list_derived_strings <-
    as.character(
      purrr::map_if(
        .x = temp_list,
        .p = was_the_value_uniqued,
        .f = function(x) {
          paste(x, collapse = regex_pattern_to_split_with)
        },
        .else = function(x) {
          return('')
        }
      )
    )

  proper_strings <-
    as.character(
      purrr::map_if(
        .x = temp_list,
        .p = !was_the_value_uniqued,
        .f = function(x) {
          return(x)
        },
        .else = function(x) {
          return('')
        }
      )
    )

  resulting_char_vector <-
    paste0(bullshit_list_derived_strings,
           proper_strings)

  return(resulting_char_vector)
}

















return_all_usable_id_types <- function(ENSG_, ENST_, Gene_ID, NM_, Accession, Unigene, NR_, XM_, XR_, list_LIST_DATA_)
{
  IDs <- rlist::list.clean(list('ensembl_gene_id' = ENSG_, 'ensembl_transcript_id' = ENST_, 'entrezgene_id' = Gene_ID, 'refseq_mrna' = NM_, 'embl' = Accession, 'mgi_id' = Unigene, 'refseq_ncrna' = NR_, 'refseq_mrna_predicted' = XM_, 'refseq_ncrna_predicted' = XR_))
  
  id_types <- purrr::map2(.x = IDs, .y = names(IDs), .f = function(x,y) 
  {
    data_for_exp <- list('col_name' = x, 'filter_name' = y)
    
    vals_in_this_col <- unique(list_LIST_DATA_[[data_for_exp$col_name]])
    
    vals_in_this_col <- vals_in_this_col[!is.na(vals_in_this_col)]
    
    if (length(vals_in_this_col) == 0) {
      return(NULL)
    } else{
      return(data_for_exp)
    }
  })
  
  id_types <- rlist::list.clean(id_types)

  return(id_types)
}









































prepare_db_for_pseudo_memoization <- function(LIST_DATA___, des_species_vec, ENSG_, ENST_, Gene_ID, NM_, Accession, Unigene, NR_, XM_, XR_, return_qa_of_pseudomemoization_, string_separator__)
{
  if(length(LIST_DATA___) != length(des_species_vec))
  {
    stop('number of experiments in not the same as number of pieces of informaton on species for these experiments')
  }

  species_subseting_vec <- purrr::map(
    .x = unique(des_species_vec), 
    .f = function(x){
      logic <- des_species_vec == x
      species <- x
      
      binded_list <- rlist::list.rbind(LIST_DATA___[logic])
      
      id_type_to_analyze <- return_all_usable_id_types(ENSG_ = ENSG_, ENST_ = ENST_, Gene_ID = Gene_ID, NM_ = NM_, Accession = Accession, Unigene = Unigene, XM_ = XM_, XR_ = XR_, NR_ = NR_, list_LIST_DATA_ = binded_list)

      id_columns <- purrr::map(.x = id_type_to_analyze, .f = function(x)
      {
        id_column <- subset(x = binded_list, subset = !is.na(binded_list[[x$col_name]]), select = (x$col_name)
        )
        
        id_vector <- purrr::map(
          .x = id_column,
          .f = function(x) {
            stringr::str_split(string = x,
                               pattern = string_separator__,
                               simplify = T)})
        
        id_vector <- unique(as.character(rlist::list.cbind(id_vector)))
        
        id_vector <- id_vector[order(id_vector)]
        
        return(list(id_column, id_vector))
      }
      )
      
      return(list(
          'logic' = logic,
          'species' = species,
          'binded_list' = binded_list,
          'id_type_to_analyze' = id_type_to_analyze,
          'id_columns' = id_columns))
    })
  
  
  
  if (return_qa_of_pseudomemoization_ == T) {
    return(species_subseting_vec)
  } else {
    meaty_return <- purrr::map(
      .x = species_subseting_vec,
      .f = function(x) {
        return_ <- list('species' = x$species,
                        'columns' = x$id_columns)
        
        return_$columns <- purrr::map(
          .x = return_$columns,
          .f = function(z) {
            z[[2]]})
        
        return(return_)
      }
    )
    
    return(meaty_return)
  }
}





















pseudo_memoize <- function(pseudo_memoization_db_)
{
  if (is.null(pseudo_memoization_db_)) {
    stop('pseudo_memoization_db_ missing')
  }
  
  memoized <- purrr::map(.x = pseudo_memoization_db_, .f = function(species__)
  {
    mart_for_using <- set_mart_to_be_used(species_ = species__$species, int_loop = 0)
    memoized_id <- purrr::map2(.x = species__$columns, .y = names(species__$columns), .f = function(column__, col_name)
      {
      memoized_data <- biomaRt::getBM(
        attributes = c(col_name, "external_gene_name"),
        filters = col_name,
        values = column__,
        uniqueRows = T,
        mart = mart_for_using)
      
      memoized_data <- subset(x = memoized_data, subset = !is.na(memoized_data[["external_gene_name"]]))
      
      memoized_data <- subset(x = memoized_data, subset = memoized_data[["external_gene_name"]] != '')
      
      ### !!! Added this, so that the ids that are numbers do not act differently further down the road 280521
      memoized_data[[1]] <- as.character(memoized_data[[1]])
      memoized_data[[2]] <- as.character(memoized_data[[2]])
      ### !!! Added this, so that the ids that are numbers do not act differently further down the road 280521
      
      return(memoized_data)
      })
  })
  
  names_for_memo <- purrr::map_chr(.x = pseudo_memoization_db_, .f = function(species__){
    species__$species
  })
  
  names(memoized) <- names_for_memo
  
  return(memoized)
}



















# list_gemma_platforms is list of annotations for platforms in list_LIST_DATA downloaded from gemma: a result from download_platforms_from_gemma function
get_gemma_annotations_for_data <- function(
  descriptions_df,
  des_platform_col,
  named_list_gemma_platforms, #from download_platforms_from_gemma
  gemma_probe_col,
  input_df,
  exp_id_col,
  input_probe_col)
{
  
  ### Make sure input is proper ###
  validate_col_types(df_ = descriptions_df, col_names_list = list(exp_id_col, des_platform_col), col_types_list = list('numeric', 'character'))
  
  validate_col_types(df_ = input_df, col_names_list = list(exp_id_col, input_probe_col), col_types_list = list('numeric', 'character'))
  
  purrr::walk(.x = named_list_gemma_platforms, .f = function(x) {
    if (!is.na(x)) {
      validate_col_types(df_ = x, col_names_list = list(gemma_probe_col), col_types_list = list('character'))
    }
  })

  named_list_gemma_platforms_no_na <- list()
  
  for (list_nb in seq_along(named_list_gemma_platforms)) {
    if (length(named_list_gemma_platforms[[list_nb]] == 1)) {
      if (!is.na(named_list_gemma_platforms[[list_nb]])) {
        named_list_gemma_platforms_no_na[[list_nb]] <- named_list_gemma_platforms[[list_nb]]
        names(named_list_gemma_platforms_no_na)[[list_nb]] <- names(named_list_gemma_platforms)[[list_nb]]
      }
    }
    
  }
  
  composite_list <- list()
  
  for (platform_nb in seq_along(named_list_gemma_platforms_no_na)) {
    
    platform <- names(named_list_gemma_platforms_no_na)[[platform_nb]]
    
    composite_list[[platform]]$platform_data <- named_list_gemma_platforms_no_na[[platform]]
    
    composite_list[[platform]]$descriptions <- subset(
      x = descriptions_df, 
      subset = toupper(descriptions_df[[des_platform_col]]) == platform)
    
    
    composite_list[[platform]]$input_list <- purrr::map(
      .x = unique(composite_list[[platform]]$descriptions[[exp_id_col]]), 
      .f = function(.x) {
        subset(
          x = input_df, 
          subset = input_df[[exp_id_col]] == .x)
      })
    
    composite_list[[platform]]$input <- rlist::list.rbind(composite_list[[platform]]$input_list)
  }
  
  
  return_list <-
    purrr::map(
      .x = composite_list,
      .f = function(x)
      {
        merge(
          x = x$input,
          y = x$platform_data,
          by.x = input_probe_col,
          by.y = gemma_probe_col,
          all.x = T
        ) #%>%
          #dplyr::select(Probe_ID:GeneSymbols)
      }
    )

  return(return_list)
}




















# str_platforms_ids_to_download are in GEO format
download_platforms_from_gemma <- function(str_platforms_ids_to_download)
{
  #devtools::install_github('PavlidisLab/gemmaAPI.R')
  
  str_platforms_ids_to_download <- toupper(str_platforms_ids_to_download)
  
  username_ <- getPass::getPass(msg = "Gimme Your GEMMA username: ")
  password_ <- getPass::getPass(msg = "Gimme Your GEMMA password: ")
  
  gemmaAPI::setGemmaUser(username = username_, password = password_)
  
  temp_list <- list()
  
  temp_list <- lapply(X = str_platforms_ids_to_download, FUN = function(X)
  {
    tryCatch(
      {
        gemmaAPI::platformInfo(platform = X, 
                               request = 'annotations')
      }, error = function(content){
        warning(paste0('Platform ', X, ' not found'))
        return(NA)})
    
    
  })
  
  gemmaAPI::setGemmaUser()
  
  names(temp_list) <- str_platforms_ids_to_download
  
  return(temp_list)
}

















GPL8160_post_gemma_wrapper <- function(annotation_raw, divider_of_values_in_serie_str_, old_divider_to_be_replaced, uniqualize_, Gene_ID_col, gene_description_col, gene_symbol_col, cols_to_nullify)
{
  annotation_raw[[gene_description_col]] <- annotation_raw$GeneNames
  
  annotation_raw$GeneSymbols <- stringr::str_replace_all(string = annotation_raw$GeneSymbols, pattern = old_divider_to_be_replaced, replacement = divider_of_values_in_serie_str_)
  
  annotation_raw$NCBIids <- stringr::str_replace_all(string = annotation_raw$NCBIids, pattern = old_divider_to_be_replaced, replacement = divider_of_values_in_serie_str_)
  
  annotation_raw[[gene_symbol_col]] <- paste_strings_consisting_of_series_of_values(strings_to_add_list = list(annotation_raw[[gene_symbol_col]], annotation_raw$GeneSymbols), divider_of_values_in_serie_str = divider_of_values_in_serie_str_, uniqualize = uniqualize_)
  
  annotation_raw[[Gene_ID_col]] <- paste_strings_consisting_of_series_of_values(strings_to_add_list = list(annotation_raw[[Gene_ID_col]], annotation_raw$NCBIids), divider_of_values_in_serie_str = divider_of_values_in_serie_str_, uniqualize = uniqualize_)
  
  for(col in cols_to_nullify)
  {
    annotation_raw[[col]] <- NULL
  }
  
  return(annotation_raw)
  
}























annotate_identifiers_to_geneID <- function(
  descriptions_df_,
  desc_exp_id_col_,
  desc_species_col_,
  input_df_,
  input_exp_id_col_,
  input_id_col_,
  chr_gene_identifier_type = 'Gene name',
  string_separator_)
{

  input_plus_species <- merge(x = input_df_, y = descriptions_df_, by.x = input_exp_id_col_, by.y = desc_exp_id_col_)

  unique_input_plus_species <- get_vector_of_single_unique_gene_ids_and_species(
    input_df = input_plus_species, 
    identifer_col = input_id_col_, 
    species_col = desc_species_col_, 
    string_separator = string_separator_)
  
  query_for_ncbi <- get_query_for_ncbi_geneID_annotation(
    char_vec_gene_id_to_query_with = unique_input_plus_species[[input_id_col_]], 
    char_vec_organism = unique_input_plus_species[[desc_species_col_]], 
    chr_gene_identifier = chr_gene_identifier_type) ### !!! Tutaj chyba się psuje, jak query_for_ncbi == 0 length df

  strvec_ncbi_query_for_identifers <- search_for_ids_in_ncbi(query_for_ncbi)

  annotated_data <- make_and_write_table_with_original_and_ncbi_ids(
    entrez_gene_search_output = strvec_ncbi_query_for_identifers, 
    df_original_data = input_plus_species)

  return(annotated_data)
}
















# INPUT: two corresponding char vectors, char_vec_gene_id_to_query_with - gene names to query ncbi with, char_vec_organism - species name for every gene name. OUTPUT: List. [[1]] - actual queries to be sent to ncbi; [[2]] - vector added to Probe_ID (gene name). Without it, it is difficult for next function to return probeid/genename to its original form, which is neccesary for merging queries with original dataframe
get_query_for_ncbi_geneID_annotation <- function(
  char_vec_gene_id_to_query_with, 
  char_vec_organism, 
  chr_gene_identifier)# 'Gene name'
{

  normalized_species_names <- normalize_species_names(return_normalized_names = T)

  char_vec_organism <- stringr::str_replace(string = char_vec_organism, pattern = normalized_species_names[['mouse']], replacement = 'mus musculus')
  char_vec_organism <- stringr::str_replace(string = char_vec_organism, pattern = normalized_species_names[['rat']], replacement = 'rattus norvegicus')
  char_vec_organism <- stringr::str_replace(string = char_vec_organism, pattern = normalized_species_names[['human']], replacement = 'homo sapiens')
  char_vec_organism <- stringr::str_replace(string = char_vec_organism, pattern = normalized_species_names[['sheep']], replacement = 'ovis aries')
  char_vec_organism <- stringr::str_replace(string = char_vec_organism, pattern = normalized_species_names[['saimiri']], replacement = 'saimiri')
  
  linker_string_crucial_for_returning_ProbeID_to_proper_form <- paste0('[', chr_gene_identifier, '] AND ') 

  query_vector <- paste0(char_vec_gene_id_to_query_with, linker_string_crucial_for_returning_ProbeID_to_proper_form, char_vec_organism, '[Organism]') 

  return(query_vector)
}




















### Annotate unidentified, !!!but unique!!! gene identifiers to gene-id using ncbi gene database. INPUT: delay_between_queries - 0.5 produces errors even. 3 times per second my asshole. We can query ncbi databases only 3 times per second - that is the reason for sys.sleep time! A
search_for_ids_in_ncbi <- function(query_vector, entrez_db = 'gene', delay_between_queries = 0.5)
{
  ids <- list()
  counter <- 1
  for (id_ in query_vector)
  {
    done <- F
    
    while (done == F) {
      message(id_)
      message(counter)
      
      tryCatch( {
        ids[[counter]] <- rentrez::entrez_search(db = entrez_db, term = id_)
        
        done <- T
        
        counter = counter + 1
        
        Sys.sleep(delay_between_queries)
      },
      error=function(e) {
        Sys.sleep(5)
        warning('Failed to download data, trying again in 5 s')
      })
    }
  }
  return(ids)
}





















#This function should except output of search_for_ids function. BEWARE - this only returns first geneID found INPUT. 
make_and_write_table_with_original_and_ncbi_ids <- function(
  entrez_gene_search_output, 
  df_original_data)
{
  
  extracted_data <- purrr::map(.x = entrez_gene_search_output, .f = function(x){
    if (length(x$ids) != 0) {
      entrez_id <- x$ids[[1]]
    } else { entrez_id <- NA }
    
    
    query <- stringr::str_split(string = x$QueryTranslation, pattern = ' AND ', simplify = T)
    
    input_id <- stringr::str_remove(string = query[1], pattern = '\\[.*')
    
    organism <- stringr::str_remove(string = tolower(query[2]), pattern = '\\[.*')
    organism <- stringr::str_remove_all(string = organism, pattern = '"')
    
    return(c(entrez_id, input_id, organism))
    
  })
  
  extracted_data_df <- as.data.frame(t(as.data.frame(extracted_data)))
  
  colnames(extracted_data_df) <- c('Gene_ID', 'input_id', 'organism')
  rownames(extracted_data_df) <- seq_along(extracted_data_df[[1]])
  
  extracted_data_df$dummy_exp <- as.numeric(as.factor(extracted_data_df$organism))
  
  extracted_data_df$Gene_ID <- as.character(extracted_data_df$Gene_ID)
  extracted_data_df$input_id <- as.character(extracted_data_df$input_id)
  extracted_data_df$organism <- as.character(extracted_data_df$organism)

  extracted_data_df$organism <- normalize_species_names(extracted_data_df$organism)

  return(extracted_data_df)
}

















get_vector_of_single_unique_gene_ids_and_species <- function(
  input_df,
  identifer_col,
  species_col,
  string_separator
)
{
  split_input_df <- split(input_df, f = input_df[[species_col]])
  
  uniquelized_ids <- purrr::map(.x = split_input_df, .f = function(x) {
    temp <- tibble::tibble(
      split_string_by_pattern_and_extract_vec_of_unique_values(chr_vec = x[[identifer_col]], pattern_to_split_individual_strings_with = string_separator),
      x[[species_col]][1])
    
    return(temp)
  })

  uniquelized_ids_df <- rlist::list.rbind(uniquelized_ids)
  
  colnames(uniquelized_ids_df) <- c(identifer_col, species_col)
  
  uniquelized_ids_df <- subset(x = uniquelized_ids_df, subset = !is.na(uniquelized_ids_df[[identifer_col]]))
  
  return(uniquelized_ids_df)
}























add_new_gene_id_col_originating_from_ncbi_annotation <- function(
  descriptions_df,
  desc_exp_id_col,
  organism_col,
  input_df, 
  input_input_id_col,
  input_exp_id_col,
  perform_ncbi_annotation_output, 
  output_gene_id_col = 'Gene_ID', 
  output_input_id_col = 'input_id', 
  output_organism_col = 'organism', 
  new_col_name = 'Gene_ID_from_ncbi',
  str_sep){

  perform_ncbi_annotation_output[[output_organism_col]] <- normalize_species_names(perform_ncbi_annotation_output[[output_organism_col]])

  input_df[[organism_col]] <- recode_values_based_on_key(
    to_recode_chrvec = as.character(input_df[[input_exp_id_col]]), 
    replace_this_chrvec = as.character(descriptions_df[[desc_exp_id_col]]), 
    with_this_chrvec = descriptions_df[[organism_col]])
  
  input_list <- list()
  output_list <- list()
  
  input_list <- purrr::map(.x = unique(input_df[[organism_col]]), .f = function(species){
    input <- subset(x = input_df, subset = input_df[[organism_col]] == species)
    output <- subset(x = perform_ncbi_annotation_output, subset = perform_ncbi_annotation_output[[output_organism_col]] == species)
    
    input[[new_col_name]] <- NA
    
    for (row_nb in seq_along(input[[1]])) {
      
      input[[new_col_name]][[row_nb]] <- split_string_by_pattern_and_replace_values_according_to_key(
        string_to_split = input[[input_input_id_col]][[row_nb]], 
        pattern_to_split_with = str_sep, 
        key_replace_this_chrvec = output[[output_input_id_col]], 
        key_with_this_chrvec = output[[output_gene_id_col]], 
        return_unique_values = T)
    }
    
    return(input)
  })
  
  input_df <- rlist::list.rbind(input_list)
  
  input_df[[organism_col]] <- NULL
  
  return(input_df)
}







assert_desc_and_input_have_the_same_exp_col_wrapper <- function(data_, desc_, exp_col)
{
  data_exps <- subset(x = data_, select = exp_col)
  desc_exps <- subset(x = desc_, select = exp_col)
  
  data_exps <- unique(data_exps)
  desc_exps <- unique(desc_exps)

  are_exps_the_same <- all(data_exps[[1]] %in% desc_exps[[1]])
  
  if (are_exps_the_same != T) {
    warning('experiments in input data are in descriptions')
  }
  
  if (class(data_exps[[exp_col]]) != class(desc_exps[[exp_col]])) {
    warning('experiments columns in input and descriptions are of different class')
  }

  return(are_exps_the_same)
}






















check_nb_of_options_selected <- function(A, B, C, D, E)
{
  options_selected <- as.logical(c(A, B, C, D, E))
  
  nb_T <- options_selected[options_selected == T]
  
  if (length(nb_T) == 1) {
    return(T)
  } else stop('you need to select exactly 1 PERFORM_* argument')
}























get_usable_symbols_subset_of_leftovers <- function(true_leftovers_, input_symbol_col, input_geneid_col = NA, other_input_cols = NA, output_ext_gene_name_col)
{
  identifiers <- as.character(c(input_symbol_col, input_geneid_col, other_input_cols))
  
  identifiers <- identifiers[!is.na(identifiers)]
  
  usable_leftovers <- list()
  leftovers_ <- true_leftovers_

  for (identifier in identifiers) {
    usable_leftovers[[identifier]] <- subset(
      x = leftovers_, 
      subset = !is.na(leftovers_[[identifier]]))
    
    usable_leftovers[[identifier]][[output_ext_gene_name_col]] <- usable_leftovers[[identifier]][[identifier]]
    
    
    if (identifier == input_geneid_col) {
      usable_leftovers[[identifier]][[output_ext_gene_name_col]] <- stringr::str_remove(
        string = usable_leftovers[[identifier]][[output_ext_gene_name_col]], 
        pattern = ',.*')
      
      usable_leftovers[[identifier]][[output_ext_gene_name_col]] <- paste0(
        'LOC', 
        usable_leftovers[[identifier]][[output_ext_gene_name_col]])
      
      } else if (identifier != input_symbol_col) {
        usable_leftovers[[identifier]][[output_ext_gene_name_col]] <- stringr::str_remove(
          string = usable_leftovers[[identifier]][[output_ext_gene_name_col]], 
          pattern = ',.*')
      }
    
    if (length(usable_leftovers[[identifier]][[1]]) == 0) {
      usable_leftovers[[identifier]] <- NULL }

    leftovers_ <- subset(
      x = leftovers_, 
      subset = is.na(leftovers_[[identifier]]))
  }

  usable_leftovers <- rlist::list.rbind(usable_leftovers)

  
  assertthat::assert_that(
    length(true_leftovers_[[1]]) == length(usable_leftovers[[1]]) + length(leftovers_[[1]]) )
  
  return(list('usable' = usable_leftovers, 'unusable' = leftovers_))
}









### !!! added 280421
# INPUT: single string containing multiple gene names separeted by separator
select_best_geneName_wrapper_for_single_string <- function(string_to_be_vectorised, separator = ', ')
{
  string_to_be_vectorised <- stringr::str_remove(string = string_to_be_vectorised, pattern = "^, ") ### These symbols were returned
  
  vectorised_string <- as.character(stringr::str_split(
    string = string_to_be_vectorised,
    pattern = separator,
    simplify = T
  ))
  
  the_best_name_ <- select_best_geneName(char_vec = vectorised_string)
  
  return(the_best_name_)
}



# INPUT: char_vec - includes all gene names returned by ensembl for given gene. Hence, the function should be iterated over list of vectors, each vector/list element for single gene. regex_to_detect_bad_names_with - genes can have 4 numbers: Olr1237 OUTPUT: the_best_name - char_vec of length 1
select_best_geneName <- function(char_vec, regex_to_detect_bad_names_with = '(\\d{5})|(^Gm\\d)', regex_to_detect_less_bad_names_with = '(Rik)|(LOC)|(Gm)', ignore_case_ = T)
{
  log_is_this_name_bad <-
    stringr::str_detect(string = char_vec, pattern = regex_to_detect_bad_names_with)
  
  good_names <- subset(x = char_vec, subset = !log_is_this_name_bad)
  bad_names <- subset(x = char_vec, subset = log_is_this_name_bad)
  
  log_is_this_name_less_bad <-
    stringr::str_detect(
      string = bad_names,
      pattern = stringr::regex(regex_to_detect_less_bad_names_with, ignore_case = ignore_case_)
    )
  
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
### !!! added 280421
