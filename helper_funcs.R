
# Load necessary libraries
library(tidyverse)

# Set the path to the Shiny app directory
app_path <- 'G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/'

# Source the file containing the database connection functions
source(file.path(app_path, 'sql_db/db_connect.R'))


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
                           top_n_by_lesion_cor= NULL,
                           signif_cor_only = FALSE,
                           signif_cor_DEGs_filtered_only = FALSE,
                           fungi = c('Bc', 'Ss'),
                           return_long = TRUE
) {
  # Process fungi names and filter allowed ones
  fungi = tolower(fungi) %>% gsub(' ', '', .) %>% gsub('\\.|_', '', .)
  fungi <- fungi[fungi %in% allowed_fungi]


  genes_query <- get_genes_query(genes)
  # Check if there are valid fungi and genes
  if (length(fungi) < 1 | is.null(genes_query)) {
    return(NULL)
  }

  WHERE <- paste0('WHERE GeneID IN ',genes_query)
  if(signif_cor_only ){#& !signif_cor_DEGs_filtered_only & is.null(abs_lesion_cor_thresh)){
    WHERE <- paste0(WHERE,"AND GeneID IN (SELECT GeneID FROM sclero_divset_corr)")
  }else if(top_n_by_lesion_cor){
    WHERE <- paste0("WHERE GeneID IN (SELECT GeneID FROM sclero_divset_corr " ,
                    "WHERE confounding_DE=0 AND GeneID IN ",genes_query,
                    "ORDER BY ABS(sclero_lesion_cor) DESC LIMIT ",top_n_by_lesion_cor, ")")
  }else if(signif_cor_DEGs_filtered_only){
    WHERE <- paste0(WHERE,"AND GeneID IN (SELECT GeneID FROM sclero_divset_corr WHERE confounding_DE=0)")
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


# Function to get time-series expression data for specified genes and fungi
get_timeseries_exp <- function(genes = c('Lsat_1_v5_gn_3_121961', 'Lsat_1_v5_gn_9_69201', 'invalid_id'),
                               fungi = c('Bc', 'Ss', 'invalid_fungi'),
                               include_raw_points = TRUE,
                               include_mock = TRUE,
                               include_conf_intervals = TRUE,
                               ts_DEGs_only = FALSE,
                               up_down_only = 'both'
) {
  # Process fungi names and filter allowed ones
  fungi = tolower(fungi) %>% gsub(' ', '', .) %>% gsub('\\.|_', '', .)
  fungi <- ifelse(fungi[fungi %in% allowed_fungi] %in% allowed_sclero, 'Scl', 'Bot')

  # Get the formatted query for valid genes
  genes_query = get_genes_query(genes)

  # Check if there are valid fungi and genes
  if (length(fungi) < 1 | is.null(genes_query)) {
    return(NULL)
  }

  # Construct the SQL query based on input parameters
  SELECT <- paste0("SELECT GeneID, Fungi, Treatment, hpi, ",
                   ifelse(include_raw_points, "rep1_log2 as rep1, rep2_log2 as rep2, rep3_log2 as rep3, ", ""),
                   "(COALESCE(rep1_log2, 0) + COALESCE(rep2_log2, 0) + COALESCE(rep3_log2, 0)) /
                   (CASE WHEN rep1_log2 IS NOT NULL THEN 1 ELSE 0 END +
                   CASE WHEN rep2_log2 IS NOT NULL THEN 1 ELSE 0 END +
                   CASE WHEN rep3_log2 IS NOT NULL THEN 1 ELSE 0 END) AS mean_expr")

  WHERE <- paste0(" WHERE GeneID IN ", genes_query)
  if (!include_mock) {
    WHERE <- paste0(WHERE, " AND Treatment = 'Infected'")
  }
  if (length(fungi) == 1) {
    WHERE <- paste0(WHERE, " AND Fungi = '", fungi, "'")
  }
  up_down_only <- tolower(trimws(up_down_only))
  if (ts_DEGs_only|up_down_only %in% c('up','down')) {
    if(up_down_only %in% c('up','down')){
      subquery_where <- paste0(" WHERE Bot_Inf_De_Dir = '",up_down_only,"' ")
    }else{
      subquery_where <- ""
    }
    WHERE <- paste0(WHERE, " AND GeneID IN (SELECT GeneID FROM timeseries_overlap_degs",subquery_where,")")
  }

  my_query <- paste0(SELECT, " FROM timeseries", WHERE)

  print(my_query)

  # Execute the SQL query and return data frame
  exp <- db_query(my_query)

  # Add confidence intervals if required
  if (include_conf_intervals) {
    confidence_intervals <- t(apply(exp[, c('rep1', 'rep2', 'rep3')], 1, calculate_confidence_interval))
    exp$'LCL' <- confidence_intervals[, 1]
    exp$'UCL' <- confidence_intervals[, 2]
  }

  return(exp)
}


GO_descrip_2_ID <- function(description) {
  description <- tolower(trimws(description))
  res <- AnnotationDbi::select(GO.db::GO.db,
                               keys=description,
                               keytype="TERM",
                               columns=c("GOID"))
  return(res$GOID[1])
}


# Function to retrieve lettuce genes associated with a given Gene Ontology (GO) term
get_GO_genes <- function(go_id = 'GO:0009723',
                         subset_lettuce_ids = NULL,
                         ts_degs_only = FALSE,
                         include_gene_names =TRUE
) {

  # Check if the input is a GO ID or a description using regex
  is_go_id <- grepl("^GO:\\d{7}$", go_id)

  # If it's a description, get its GO ID
  if (!is_go_id) {
    go_id <- GO_descrip_2_ID(go_id)

    # Handle invalid descriptions (i.e., when the returned GO ID is NULL or has no rows)
    if (!grepl("^GO:\\d{7}$", go_id)) {
      stop("Invalid GO description provided.")
    }
  }



  # Retrieve AGI codes associated with the GO term from the TAIR database
  At_genes_w_GO <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                                         keys = go_id,
                                         columns = c("TAIR"),
                                         keytype = "GO")

  if (nrow(At_genes_w_GO) < 2) {
    return(NULL)
  }

  WHERE <- paste('WHERE AGI IN', get_genes_query(At_genes_w_GO$TAIR))

  if (ts_degs_only) {
    WHERE <- paste0(WHERE, " AND GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")
  } else if (!is.null(subset_lettuce_ids)) {
    let_genes_query <- get_genes_query(subset_lettuce_ids)
    if (!is.null(let_genes_query)) {
      WHERE <- paste0(WHERE, " AND GeneID IN ", let_genes_query)
    }
  }



  # Construct the SQL query and retrieve lettuce genes with the GO term

  query <- paste0("SELECT GeneID",
                  ifelse(include_gene_names,
                         paste0(", CASE WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID ",
                                "ELSE 'Ls' || At_ShortName END AS name"),
                         ""),
                  " FROM gene_annotations ", WHERE)

  Ls_genes_w_GO <- db_query(query)

  return(Ls_genes_w_GO)
}


# Function to retrieve lettuce genes associated with a given protein domain
get_genes_w_domain <- function(domain_id = 'PF00067',
                               domain_desc = NULL,
                               subset_lettuce_ids = NULL,
                               ts_degs_only = FALSE,
                               include_gene_name = TRUE
) {
  # Construct the WHERE clause for the SQL query based on input parameters
  if (!is.null(domain_desc)) {
    WHERE = paste("WHERE LOWER(da.description) LIKE", sprintf("'%%%s%%'", tolower(domain_desc)))
  } else {
    WHERE = paste0("WHERE da.id =", "'", domain_id, "'")
  }

  if (ts_degs_only) {
    WHERE <- paste0(WHERE, " AND da.GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")
  } else if (!is.null(subset_lettuce_ids)) {
    let_genes_query <- get_genes_query(subset_lettuce_ids)
    if (!is.null(let_genes_query)) {
      WHERE <- paste0(WHERE, " AND da.GeneID IN ", let_genes_query)
    }
  }

  # Construct the SQL query based on whether to include gene names or not
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

  # Execute the SQL query and return data frame
  Ls_genes_w_domain <- db_query(query)

  return(Ls_genes_w_domain)
}


# Function to retrieve lettuce genes that are orthologs of a given AGI (Arabidopsis thaliana gene identifier)
get_orthologs <- function(arabidopsis_genes,
                          subset_lettuce_ids = NULL,
                          ts_degs_only = FALSE,
                          include_gene_name = FALSE

) {

  arabidopsis_genes <- toupper(trimws(arabidopsis_genes))

  is_AGI =  grepl('AT[0-9]G[0-9]+', arabidopsis_genes)
  AGI <- arabidopsis_genes[is_AGI]



  # If shortname is provided, try to convert it to AGI using TAIR database
  if (sum(!is_AGI)>=1) {
    tair_code <- tryCatch({
      AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                            keys = arabidopsis_genes[!is_AGI],
                            keytype = "SYMBOL",
                            columns = "TAIR")$TAIR
    }, error = function(e) {
      # If an error occurs, return an empty data frame
      message("Error occurred while querying the TAIR database.")
      return('')
    })

    AGI <- c(AGI, tair_code)
  }

  gene_name_query = ifelse(include_gene_name,
         ",
  CASE
    WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID
    ELSE 'Ls' || At_ShortName
  END AS name","")


  # Construct the SQL query and retrieve lettuce genes with orthologous relationship
  query <- paste0("SELECT GeneID",gene_name_query,
                  " FROM gene_annotations WHERE AGI IN ",
                  get_genes_query(AGI),
                  ifelse(include_gene_name, " ORDER BY name",""))

  Ls_orthologs <- db_query(query)

  return(Ls_orthologs)
}

get_lettuce_genes_from_inputs <- function(GeneIDs, At_orthologs, GO_id, protein_domain) {
  if (!is.null(GeneIDs)) {
    lettuce_genes <- GeneIDs[grepl('Lsat_1_v5_gn_\\d_\\d+',GeneIDs)]
    return(GeneIDs)
  } else if (!is.null(At_orthologs)) {
    return(get_orthologs(arabidopsis_genes = At_orthologs)$GeneID)
  } else if (!is.null(GO_id)) {
    return(get_GO_genes(GO_id, ts_degs_only = FALSE, include_gene_names = FALSE)$GeneID)
  } else if (!is.null(protein_domain)) {
    if (grepl('PF\\d+|PTHR\\d+', protein_domain)) {
      return(get_genes_w_domain(domain_id = protein_domain,
                                include_gene_name = FALSE, ts_degs_only = FALSE)$GeneID)
    } else {
      return(get_genes_w_domain(domain_desc = protein_domain,
                                include_gene_name = FALSE, ts_degs_only = FALSE)$GeneID)
    }
  } else {
    return(NULL)
  }
}

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

# get_GO_orthologs(go_id = 'GO:0009723', ts_degs_only=TRUE)
#test <- get_genes_w_domain(domain_desc = 'AP2', ts_degs_only=TRUE)
#get_orthologs(shortname = 'LBD29')


