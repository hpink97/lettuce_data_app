
# Load necessary libraries
library(tidyverse)

# Set the path to the Shiny app directory
#app_path <- 'G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/'

# Source the file containing the database connection functions
#source(file.path(app_path, 'sql_db/db_connect.R'))
#source('sql_db/db_connect.R')

# Define vectors of allowed fungi identifiers
allowed_bot <- c('bot', 'botrytis', 'b.cinerea', 'bc', 'bcin')
allowed_sclero <- c('scl', 'sclero', 'sclerotinia', 'ss', 'ssclerotiorum', 'ssclerotiroum')
allowed_fungi <- c(allowed_bot, allowed_sclero)


# Function to get a formatted query for valid genes
get_genes_query <- function(input_genes) {
  # Filter out genes with specific patterns using regex
  valid_genes <- input_genes[grepl('Lsat_1_v5_gn_\\d|AT\\d[g|G]\\d+', input_genes)]
  valid_genes <- unique(valid_genes)
  genes_query <- paste0("('", paste(valid_genes, collapse = "', '"), "')")

  if (length(valid_genes) >= 1) {
    return(genes_query)
  } else {
    return(NULL)
  }
}


# Function to retrieve expression data from the database for a given set of genes and fungi
get_divset_exp <- function(genes = c('Lsat_1_v5_gn_5_78841', 'Lsat_1_v5_gn_9_69201', 'Lsat_1_v5_gn_1_56161'),
                           top_n_by_lesion_cor= 1e6,
                           signif_cor_only = FALSE,
                           DEGs_filtered= TRUE,
                           fungi = c('Bc', 'Ss'),
                           return_long = TRUE
) {
  # Process fungi names and filter allowed ones
  fungi = tolower(fungi) %>% gsub(' ', '', .) %>% gsub('\\.|_', '', .)
  fungi <- fungi[fungi %in% allowed_fungi]


  if(top_n_by_lesion_cor>length(unique(genes))){
    top_n_by_lesion_cor <- NULL
  }



  genes_query <- get_genes_query(genes)
  # Check if there are valid fungi and genes
  if (length(fungi) < 1 | is.null(genes_query)) {
    return(NULL)
  }

  WHERE <- paste0('WHERE GeneID IN ',genes_query)
  if(!is.null(top_n_by_lesion_cor)){
    WHERE <- paste0("WHERE GeneID IN (SELECT GeneID FROM sclero_divset_corr " ,
                    "WHERE ",
                    ifelse(DEGs_filtered,"confounding_DE=0 AND ",""),
                    "GeneID IN ",genes_query,
                    "ORDER BY ABS(sclero_lesion_cor) DESC LIMIT ",top_n_by_lesion_cor, ")")
  }else if(signif_cor_only & !DEGs_filtered ){#& !signif_cor_DEGs_filtered_only & is.null(abs_lesion_cor_thresh)){
    WHERE <- paste0(WHERE,"AND GeneID IN (SELECT GeneID FROM sclero_divset_corr)")
  }else if(signif_cor_only & DEGs_filtered){
    WHERE <- paste0(WHERE,"AND GeneID IN (SELECT GeneID FROM sclero_divset_corr WHERE confounding_DE=0)")
  }else if(!signif_cor_only & DEGs_filtered){
    WHERE <- paste0(WHERE,"AND GeneID NOT IN (SELECT GeneID FROM sclero_divset_corr WHERE confounding_DE=1)")
  }








  # Construct the SQL query based on the number of fungi
  if (length(fungi) == 1) {
    my_query <- paste0("SELECT * FROM ",
                       ifelse(fungi %in% allowed_sclero, 'sclero', 'bot'), "_divset_expr ", WHERE)
  } else {
    my_query <- paste0("SELECT bot.GeneID, bot.* , scl.* ",
                       "FROM (SELECT * FROM bot_divset_expr ",
                       WHERE, ") AS bot ",
                       "LEFT JOIN (SELECT * FROM sclero_divset_expr ",
                       WHERE, ") AS scl ",
                       "ON bot.GeneID = scl.GeneID")
  }

  # Execute the SQL query and return data frame
  genes_exp.wide <- db_query(my_query)
  genes_exp.wide <- dplyr::select(genes_exp.wide,
                                  all_of(unique(colnames(genes_exp.wide))))

  # Convert wide data frame to long format if required
  if (return_long) {
    genes_exp.long <- pivot_longer(genes_exp.wide,
                                   cols = -1,
                                   names_to = c('accession', 'biorep', 'fungi'),
                                   names_sep = '_',
                                   values_to = 'log2Exp')
    return(genes_exp.long)
  } else {
    return(genes_exp.wide)
  }
}


# Function to calculate confidence intervals for a given row of data
calculate_confidence_interval <- function(row) {
  selected_values <- row[!is.na(row)]
  lower_bound <- quantile(selected_values, 0.025)
  upper_bound <- quantile(selected_values, 0.975)
  return(c(lower_bound, upper_bound))
}


##----------------------------------------------------------------------------------------------------------------------------
# This function retrieves the time-series expression data for specified genes and fungi from a database.
#
# Parameters:
# - genes: Vector of gene IDs to be queried. Default values are provided as examples.
# - fungi: Vector of fungi names. Default values are provided as examples.
# - include_raw_points: Boolean flag to indicate whether raw data points should be included.
# - include_mock: Boolean flag to determine if mock treatments should be included in the results.
# - include_conf_intervals: Boolean flag to determine if confidence intervals should be computed for the results.
# - ts_DEGs_only: Boolean flag to filter only the differentially expressed genes (DEGs) in the time series.
# - up_down_only: String indicating which type of genes to retrieve based on their expression direction (upregulated, downregulated, or both).
#
# Returns:
# - A data frame with the expression data for the specified genes and fungi.
#
get_timeseries_exp <- function(genes = c('Lsat_1_v5_gn_3_121961', 'Lsat_1_v5_gn_9_69201', 'invalid_id'),
                               fungi = c('Bc', 'Ss', 'invalid_fungi'),
                               include_raw_points = TRUE,
                               include_mock = TRUE,
                               include_conf_intervals = TRUE,
                               ts_DEGs_only = FALSE,
                               up_down_only = 'both'
) {
  # Process fungi names to standardize them and filter out valid fungi.
  fungi = tolower(fungi) %>% gsub(' ', '', .) %>% gsub('\\.|_', '', .)
  fungi <- ifelse(fungi[fungi %in% allowed_fungi] %in% allowed_sclero, 'Scl', 'Bot')

  # Retrieve a formatted query string for valid genes.
  genes_query = get_genes_query(genes)

  # Check if there are valid fungi and genes. If not, return NULL.
  if (length(fungi) < 1 | is.null(genes_query)) {
    return(NULL)
  }

  ## ensure that up_down_only length ==1
  if(length(up_down_only)>1){
    up_down_only <- up_down_only[1]
  }

  # Construct the SQL query based on input parameters.
  # SELECT clause will determine which columns will be included in the result set.
  SELECT <- paste0("SELECT GeneID, Fungi, Treatment, hpi, ",
                   ifelse(include_raw_points, "rep1_log2 as rep1, rep2_log2 as rep2, rep3_log2 as rep3, ", ""),
                   "(COALESCE(rep1_log2, 0) + COALESCE(rep2_log2, 0) + COALESCE(rep3_log2, 0)) /
                   (CASE WHEN rep1_log2 IS NOT NULL THEN 1 ELSE 0 END +
                   CASE WHEN rep2_log2 IS NOT NULL THEN 1 ELSE 0 END +
                   CASE WHEN rep3_log2 IS NOT NULL THEN 1 ELSE 0 END) AS mean_expr")

  # WHERE clause will filter out the records based on the input parameters.
  WHERE <- paste0(" WHERE GeneID IN ", genes_query)

  #if mock data not needed, query just infected data
  if (!include_mock) {
    WHERE <- paste0(WHERE, " AND Treatment = 'Infected'")
  }

  ## query expression in response to specific fungi selected by user input
  if (length(fungi) == 1) {
    WHERE <- paste0(WHERE, " AND Fungi = '", fungi, "'")
  }

  ## select
  up_down_only <- tolower(trimws(up_down_only))
  if (ts_DEGs_only|up_down_only %in% c('up','down')) {
    if(up_down_only %in% c('up','down')){
      subquery_where <- paste0(" WHERE Bot_Inf_De_Dir = '",up_down_only,"' ")
    }else{
      subquery_where <- ""
    }
    WHERE <- paste0(WHERE, " AND GeneID IN (SELECT GeneID FROM timeseries_overlap_degs",subquery_where,")")
  }

  # Concatenate SELECT and WHERE clauses to form the complete query.
  my_query <- paste0(SELECT, " FROM timeseries", WHERE)

  # Print the query for debugging purposes.
  print(my_query)

  # Execute the SQL query.
  exp <- db_query(my_query)

  # Add confidence intervals if required
  if (include_conf_intervals) {
    confidence_intervals <- t(apply(exp[, c('rep1', 'rep2', 'rep3')], 1, calculate_confidence_interval))
    exp$'LCL' <- confidence_intervals[, 1]
    exp$'UCL' <- confidence_intervals[, 2]
  }

  return(exp)
}


#------------------------------------------------------------------------------------------------------------------------

# This function converts a given Gene Ontology (GO) term description to its corresponding GO ID.
#
# Parameters:
# - description: A string representing the GO term description.
#
# Returns:
# - The GO ID corresponding to the provided description.
#
GO_descrip_2_ID <- function(description) {
  # Convert description to lowercase and trim white spaces.
  description <- tolower(trimws(description))

  # Use AnnotationDbi to query the GO database for the given description and retrieve the corresponding GO ID.
  res <- AnnotationDbi::select(GO.db::GO.db,
                               keys=description,
                               keytype="TERM",
                               columns=c("GOID"))

  # Return the first GO ID from the results.
  return(res$GOID[1])
}



#------------------------------------------------------------------------------------------------------------------

# This function retrieves lettuce genes associated with a given GO term using the TAIR database.
#
# Parameters:
# - go_id: A string representing the GO ID or the GO term description.
# - subset_lettuce_ids: A vector of lettuce gene IDs to limit the results. If NULL, no subsetting is done.
# - ts_degs_only: A boolean flag indicating if only differentially expressed genes should be returned.
# - include_gene_names: A boolean flag indicating if gene names should be included in the result.
#
# Returns:
# - A data frame containing lettuce genes associated with the given GO term.
#
get_GO_genes <- function(go_id = 'GO:0009723',
                         subset_lettuce_ids = NULL,
                         ts_degs_only = FALSE,
                         include_gene_names =TRUE
) {

  # Check if the provided go_id is a GO ID format using regex.
  is_go_id <- grepl("^GO:\\d{7}$", go_id)

  # If the provided input is a GO description, convert it to its corresponding GO ID.
  if (!is_go_id) {
    go_id <- GO_descrip_2_ID(go_id)

    # Handle cases where an invalid GO description is provided (i.e., when the GO ID doesn't match the expected format).
    if (!grepl("^GO:\\d{7}$", go_id)) {
      stop("Invalid GO description provided.")
    }
  }

  # Retrieve Arabidopsis gene IDs (AGI codes) associated with the provided GO term from the TAIR database.
  At_genes_w_GO <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                                         keys = go_id,
                                         columns = c("TAIR"),
                                         keytype = "GO")

  # If fewer than 2 genes are retrieved, return NULL (indicating no associated genes found).
  if (nrow(At_genes_w_GO) < 2) {
    return(NULL)
  }

  # Construct the WHERE clause for the SQL query based on the provided parameters.
  WHERE <- paste('WHERE AGI IN', get_genes_query(At_genes_w_GO$TAIR))

  # If only differentially expressed genes should be included.
  if (ts_degs_only) {
    WHERE <- paste0(WHERE, " AND GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")
  }
  # If a subset of lettuce gene IDs is provided.
  else if (!is.null(subset_lettuce_ids)) {
    let_genes_query <- get_genes_query(subset_lettuce_ids)
    if (!is.null(let_genes_query)) {
      WHERE <- paste0(WHERE, " AND GeneID IN ", let_genes_query)
    }
  }

  # Construct the SQL query based on the constructed WHERE clause and retrieve lettuce genes associated with the provided GO term.
  query <- paste0("SELECT GeneID",
                  ifelse(include_gene_names,
                         paste0(", CASE WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID ",
                                "ELSE 'Ls' || At_ShortName END AS name"),
                         ""),
                  " FROM gene_annotations ", WHERE)

  # Execute the SQL query and store the result in Ls_genes_w_GO.
  Ls_genes_w_GO <- db_query(query)

  # Return the resulting data frame.
  return(Ls_genes_w_GO)
}



#------------------------------------------------------------------------------------------------------------------------

# This function retrieves lettuce orthologs of given Arabidopsis genes based on either AGI code or gene symbol.
#
# Parameters:
# - arabidopsis_genes: A character vector containing Arabidopsis AGI codes or gene symbols.
# - subset_lettuce_ids: A vector of lettuce gene IDs to limit the results. If NULL, no subsetting is done. (unused in this function)
# - ts_degs_only: A boolean flag indicating if only differentially expressed genes should be returned. (unused in this function)
# - include_gene_name: A boolean flag indicating if gene names should be included in the result.
#
# Returns:
# - A data frame containing lettuce genes that are orthologs of the provided Arabidopsis genes.
#
get_orthologs <- function(arabidopsis_genes,
                          subset_lettuce_ids = NULL,
                          ts_degs_only = FALSE,
                          include_gene_name = FALSE
) {

  # Convert the provided Arabidopsis genes to uppercase and trim white spaces.
  arabidopsis_genes <- toupper(trimws(arabidopsis_genes))

  # Check if the provided genes are in the AGI code format.
  is_AGI =  grepl('AT[0-9]G[0-9]+', arabidopsis_genes)

  # Filter out the genes that are in AGI code format.
  AGI <- arabidopsis_genes[is_AGI]

  # If any gene symbols are provided (i.e., not in AGI format), convert them to AGI using the TAIR database.
  if (sum(!is_AGI) >= 1) {
    tair_code <- tryCatch({
      AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                            keys = arabidopsis_genes[!is_AGI],
                            keytype = "SYMBOL",
                            columns = "TAIR")$TAIR
    }, error = function(e) {
      # If an error occurs during the database query, display an error message and return an empty string.
      message("Error occurred while querying the TAIR database.")
      return('')
    })

    # Combine the original AGI codes with the AGI codes retrieved from the database.
    AGI <- c(AGI, tair_code)
  }

  # If gene names should be included in the result, construct the respective SQL query part.
  gene_name_query = ifelse(include_gene_name,
                           ",
  CASE
    WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID
    ELSE 'Ls' || At_ShortName
  END AS name","")

  # Construct the complete SQL query to retrieve lettuce orthologs based on the provided AGI codes.
  query <- paste0("SELECT GeneID",gene_name_query,
                  " FROM gene_annotations WHERE AGI IN ",
                  get_genes_query(AGI),
                  ifelse(include_gene_name, " ORDER BY name",""))

  # Execute the SQL query and store the result.
  Ls_orthologs <- db_query(query)

  # Return the resulting data frame.
  return(Ls_orthologs)
}


#------------------------------------------------------------------------------------------------------------------------

# This function retrieves lettuce genes associated with a specified protein domain.
#
# Parameters:
# - domain_id: The protein domain ID (default is 'PF00067').
# - domain_desc: A description of the protein domain. If provided, it takes precedence over domain_id.
# - subset_lettuce_ids: A vector of lettuce gene IDs to limit the results. If NULL, no subsetting is done.
# - ts_degs_only: A boolean flag indicating if only differentially expressed genes should be returned.
# - include_gene_name: A boolean flag indicating if gene names should be included in the result.
#
# Returns:
# - A data frame containing lettuce genes associated with the specified protein domain.
#
get_genes_w_domain <- function(domain_id = 'PF00067',
                               domain_desc = NULL,
                               subset_lettuce_ids = NULL,
                               ts_degs_only = FALSE,
                               include_gene_name = TRUE
) {
  # Construct the WHERE clause for the SQL query based on input parameters

  # If a domain description is provided, use it in the WHERE clause with a LIKE statement to find matches.
  if (!is.null(domain_desc)) {
    WHERE = paste("WHERE LOWER(da.description) LIKE", sprintf("'%%%s%%'", tolower(domain_desc)))
  } else {
    # Otherwise, use the provided domain ID to filter the results.
    WHERE = paste0("WHERE da.id =", "'", domain_id, "'")
  }

  # If only differentially expressed genes should be included, add this filter to the WHERE clause.
  if (ts_degs_only) {
    WHERE <- paste0(WHERE, " AND da.GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")
  }
  # If a specific subset of lettuce genes is provided, add this filter to the WHERE clause.
  else if (!is.null(subset_lettuce_ids)) {
    let_genes_query <- get_genes_query(subset_lettuce_ids)
    if (!is.null(let_genes_query)) {
      WHERE <- paste0(WHERE, " AND da.GeneID IN ", let_genes_query)
    }
  }

  # Construct the SQL query. If gene names should be included in the result, create a JOIN to the gene_annotations table.
  if (include_gene_name) {
    query <- paste("SELECT DISTINCT da.GeneID,",
                   "CASE WHEN (ga.At_ShortName IS NULL OR ga.At_ShortName LIKE 'AT%G%') THEN da.GeneID",
                   "ELSE 'Ls' || ga.At_ShortName END AS name",
                   "FROM domain_annotations AS da",
                   "LEFT JOIN gene_annotations AS ga ON da.GeneID = ga.GeneID",
                   WHERE)
  } else {
    query <- paste0("SELECT DISTINCT GeneID FROM domain_annotations AS da ", WHERE)
  }

  # Print the generated query (optional)
  print(query)

  # Execute the SQL query and store the result.
  Ls_genes_w_domain <- db_query(query)

  # Return the resulting data frame.
  return(Ls_genes_w_domain)
}




#------------------------------------------------------------------------------------------------------------------------

# This function serves as a centralized wrapper to retrieve lettuce gene IDs based on multiple criteria.
# It prioritizes the criteria in the following order: GeneIDs > At_orthologs > GO_id > protein_domain.
# Only one criterion will be used per function call, based on the given priority.
#
# Parameters:
# - GeneIDs: A vector of provided lettuce gene IDs.
# - At_orthologs: A vector of Arabidopsis gene symbols or AGI codes for which lettuce orthologs are desired.
# - GO_id: A Gene Ontology term (either ID or description) to get associated lettuce genes.
# - protein_domain: A protein domain ID or description to get associated lettuce genes.
#
# Returns:
# - A vector of lettuce gene IDs that match the given criterion.

get_lettuce_genes_from_inputs <- function(GeneIDs=NULL, At_orthologs=NULL, GO_id=NULL, protein_domain=NULL) {

  # If GeneIDs are provided, filter those that match the lettuce gene ID format and return them.
  if (!is.null(GeneIDs)) {
    lettuce_genes <- GeneIDs[grepl('Lsat_1_v5_gn_\\d_\\d+',GeneIDs)]
    return(GeneIDs)
  }

  # If At_orthologs are provided, retrieve the lettuce orthologs of the provided Arabidopsis genes.
  else if (!is.null(At_orthologs)) {
    return(get_orthologs(arabidopsis_genes = At_orthologs)$GeneID)
  }

  # If a GO_id is provided, retrieve lettuce genes associated with the given GO term.
  else if (!is.null(GO_id)) {
    return(get_GO_genes(GO_id, ts_degs_only = FALSE, include_gene_names = FALSE)$GeneID)
  }

  # If a protein_domain is provided:
  else if (!is.null(protein_domain)) {
    # If it matches the standard protein domain ID patterns from Pfam or Panther (PFxxxxx or PTHRxxxxx),
    # retrieve lettuce genes associated with the given protein domain ID.
    if (grepl('PF\\d+|PTHR\\d+', protein_domain)) {
      return(get_genes_w_domain(domain_id = protein_domain,
                                include_gene_name = FALSE, ts_degs_only = FALSE)$GeneID)
    }
    # If the provided protein_domain doesn't match the standard ID patterns,
    # assume it's a description and retrieve genes associated with the given protein domain description.
    else {
      return(get_genes_w_domain(domain_desc = protein_domain,
                                include_gene_name = FALSE, ts_degs_only = FALSE)$GeneID)
    }
  }

  # If none of the above criteria is provided, return NULL.
  else {
    return(NULL)
  }
}



#------------------------------------------------------------------------------------------------------------------------

get_gene_names <- function(geneids, allow_dups = FALSE){
  #if no names provided, automatically generate them from annotations file
  gene_df <- db_query(paste0("SELECT GeneID,
  CASE
    WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID
    ELSE 'Ls' || At_ShortName
  END AS name
  FROM gene_annotations WHERE GeneID IN ",
                             get_genes_query(unique(geneids))))

  if(allow_dups){
    return(gene_df)
  }


  while(any(duplicated(gene_df$name))){
      dup_name <- gene_df$name [match(TRUE,duplicated(gene_df$name) )]
      dups <- filter(gene_df, name == dup_name)
      for( i in 1:nrow(dups)){
        row_i <- gene_df$GeneID == dups$GeneID[i]
        gene_df$name[row_i] <- paste0(gene_df$name[row_i], LETTERS[i])
      }
    }

  return(gene_df)
}

helperFuncsLoaded = TRUE




