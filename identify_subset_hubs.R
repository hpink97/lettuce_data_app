library(tidyverse)


if(!'helperFuncsLoaded' %in% ls()){
  funcs_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/plotting_funcs/helper_funcs.R"
  source(funcs_path)
}


get_pred_regulators<- function(GeneIDs = NULL,
                               At_orthologs = NULL,
                               GO_id = NULL,
                               protein_domain = NULL,
                               n_TFs = 5,
                               min_subset_targets = 5,
                               order_by_column = "FoldEnrichment"){ # Or any other column name you'd like){
  inputs <- c(!is.null(GeneIDs),
              !is.null(At_orthologs),
              !is.null(GO_id),
              !is.null(protein_domain))

  if(sum(inputs)!=1){return(NULL)}


  if(inputs[1]){
    target_query = get_genes_query(GeneIDs)
  }else if(inputs[2]){
    lettuce_orthologs <- get_orthologs(arabidopsis_genes = At_orthologs)
    target_query = get_genes_query(lettuce_orthologs$GeneID)
  }else if(inputs[3]){
    go_genes <- get_GO_genes(GO_id,ts_degs_only = TRUE,include_gene_names = FALSE)
    target_query <- get_genes_query(go_genes$GeneID)
  }else if(inputs[4]){
    if(grepl('PF\\d+|PTHR\\d+', protein_domain)){
      domain_genes = get_genes_w_domain(domain_id = protein_domain,include_gene_name = FALSE,ts_degs_only = TRUE)
    }else{
      domain_genes = get_genes_w_domain(domain_desc = protein_domain,include_gene_name = FALSE,ts_degs_only = TRUE)
    }

    target_query = get_genes_query(domain_genes$GeneID)
  }else{
    return(data.frame())
  }


  if(is.null(target_query)){
    return(data.frame())
  }






  query <- paste0(
    # Create a CTE (Common Table Expression) to get all unique targets from grn_edges
    "WITH AllUniqueTargets AS (",
    "SELECT DISTINCT Target FROM grn_edges),",

    # Create a CTE to get unique targets from grn_edges which are in the target_query
    "SubsetUniqueTargets AS (",
    "SELECT DISTINCT Target FROM grn_edges WHERE Target IN ", target_query, "),",

    # Compute TF statistics; count out-degrees in the subset and the total out-degrees
    "TFStats AS (",
    "SELECT TF, ",
    "COUNT(DISTINCT CASE WHEN Target IN ", target_query, " THEN Target END) AS SubsetOutdegrees,",
    "COUNT(*) AS TotalOutdegrees ",
    "FROM grn_edges GROUP BY TF),",

    # Calculate the expected out-degrees in the subset based on proportions and total out-degrees
    "ExpectedSubset AS (",
    "SELECT t.TF, ",
    "ROUND(CAST((SELECT COUNT(*) FROM SubsetUniqueTargets) AS FLOAT) / CAST((SELECT COUNT(*) FROM AllUniqueTargets) AS FLOAT) * CAST(s.TotalOutdegrees AS FLOAT), 2) AS ExpectedValue ",
    "FROM TFStats t JOIN TFStats s ON t.TF = s.TF),",

    # Extract TF names
    "TFNames AS (",
    "SELECT GeneID, ",
    "CASE ",
    "WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID ",  # if no At_Shortname use GeneID
    "ELSE 'Ls' || At_ShortName ",  # otherwise use Atshortname with prefix with 'Ls'
    "END AS TFName ",
    "FROM gene_annotations WHERE GeneID IN (SELECT DISTINCT TF FROM grn_edges)) ",  # only TFs in grn_edges are considered

    # Main query to join all CTEs and calculate metrics
    "SELECT s.TF, n.TFName, s.TotalOutdegrees, s.SubsetOutdegrees, ",
    "e.ExpectedValue AS ExpectedSubsetOutdegrees, ", # Renamed ExpectedValue for clarity
    # Calculate proportion of targets in subset over all unique targets in subset
    "ROUND(CAST(s.SubsetOutdegrees AS FLOAT)/CAST((SELECT COUNT(*) FROM SubsetUniqueTargets) AS FLOAT), 2) AS SubsetProportion,",
    # Calculate proportion of targets in subset over total out-degrees for a TF
    "ROUND(CAST(s.SubsetOutdegrees AS FLOAT)/CAST(s.TotalOutdegrees AS FLOAT), 2) AS TFProportion, ",
    # Sum and average importance, respectively
    "ROUND(CAST(SUM(Importance) AS FLOAT), 2) AS SumImportance, ",
    "ROUND(CAST(AVG(Importance) AS FLOAT), 2) AS AvgImportance, ",
    # Calculate fold enrichment as the ratio of actual to expected out-degrees in the subset
    "ROUND(CAST(s.SubsetOutdegrees AS FLOAT) / e.ExpectedValue, 2) AS FoldEnrichment ",
    "FROM grn_edges ",
    # Join with other CTEs based on TF
    "JOIN TFStats s ON grn_edges.TF = s.TF ",
    "JOIN ExpectedSubset e ON grn_edges.TF = e.TF ",
    "JOIN TFNames n ON s.TF = n.GeneID ",
    # Filter results based on target_query and minimum subset targets threshold
    "WHERE Target IN ", target_query,
    "AND s.SubsetOutdegrees >= ", min_subset_targets, " ",
    # Group by necessary columns to get aggregated metrics
    "GROUP BY s.TF, n.TFName, s.TotalOutdegrees, s.SubsetOutdegrees, e.ExpectedValue ",
    # Sort results based on a specified order column (dynamically provided)
    "ORDER BY ", order_by_column, " DESC LIMIT ", n_TFs, ";"
  )

  df <- db_query(query)


  return(df)


}

#
# get_pred_regulators(At_orthologs = c('NSL1','NSL2','RBOHD'),
#                     min_subset_targets = 2)



