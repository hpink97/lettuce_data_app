library(tidyverse)

app_path <- 'G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/'
source(file.path(app_path,'sql_db/db_connect.R'))



allowed_bot <- c('bot','botrytis','b.cinerea','bc','bcin')
allowed_sclero <- c('scl','sclero','sclerotinia','ss','ssclerotiorum','ssclerotiroum')
allowed_fungi <- c(allowed_bot,allowed_sclero)


get_genes_query <- function(input_genes){
  valid_genes <- input_genes[grepl('Lsat_1_v5_gn_\\d|AT\\d[g|G]\\d+',input_genes)]
  valid_genes <- unique(valid_genes)
  genes_query = paste0("('", paste(valid_genes, collapse = "', '"), "')")

  if(length(valid_genes)>=1){
    return(genes_query)
  }else{
    return(NULL)
  }

}


get_divset_exp <- function(genes=c('Lsat_1_v5_gn_5_78841','Lsat_1_v5_gn_9_69201','Lsat_1_v5_gn_1_56161'),
                           fungi=c('Bc','Ss'),
                           return_long = TRUE
){



  fungi = tolower(fungi) %>% gsub(' ','',.) %>% gsub('\\.|_','',.)
  fungi <- fungi[fungi %in% allowed_fungi]
  genes_query = get_genes_query(genes)


  if(length(fungi) <1 | is.null(genes_query)){
    return(NULL)
  }




  if(length(fungi) ==1){
    my_query<- paste0("SELECT * FROM ",
                    ifelse(fungi%in% allowed_sclero,'sclero','bot'),"_divset_expr WHERE GeneID IN ",genes_query)

  }else{
    my_query<- paste0("SELECT bot.GeneID, bot.* , scl.* ",
                    "FROM (SELECT * FROM bot_divset_expr ",
                    "WHERE GeneID IN ", genes_query, ") AS bot ",
                    "LEFT JOIN (SELECT * FROM sclero_divset_expr ",
                    "WHERE GeneID IN ", genes_query, ") AS scl ",
                    "ON bot.GeneID = scl.GeneID")
  }


  genes_exp.wide <- db_query( my_query)

  if(return_long){
    genes_exp.long <- pivot_longer(genes_exp.wide,
                              cols=-1,
                              names_to = c('accession','biorep','fungi'),
                              names_sep = '_',
                              values_to = 'log2Exp')
    return(genes_exp.long)
  }else{
    return(genes_exp.wide)
  }


}


calculate_confidence_interval <- function(row) {
  selected_values <- row[!is.na(row)]
  lower_bound <- quantile(selected_values, 0.025)
  upper_bound <- quantile(selected_values, 0.975)
  return(c(lower_bound, upper_bound))
}

get_timeseries_exp <- function(genes=c('Lsat_1_v5_gn_3_121961','Lsat_1_v5_gn_9_69201','invalid_id'),
                               fungi=c('Bc','Ss','invalid_fungi'),
                               include_raw_points = TRUE,
                               include_mock = TRUE,
                               include_conf_intervals = TRUE,
                               ts_DEGs_only =FALSE){

  fungi = tolower(fungi) %>% gsub(' ','',.) %>% gsub('\\.|_','',.)
  fungi <- ifelse(fungi[fungi %in% allowed_fungi] %in% allowed_sclero,
                  'Scl','Bot')
  genes_query = get_genes_query(genes)


  if(length(fungi) <1 | is.null(genes_query)){
    return(NULL)
  }


  SELECT <- paste0("SELECT GeneID, Fungi, Treatment, hpi, ",
                  ifelse(include_raw_points,
                   "rep1_log2 as rep1, rep2_log2 as rep2, rep3_log2 as rep3, ",""),
                  "(rep1_log2 + rep2_log2 + rep3_log2) / 3 AS mean_expr")

  WHERE <-paste0("WHERE GeneID IN ",genes_query)
  if(!include_mock){
    WHERE <- paste0(WHERE," AND Treatment = 'Infected'")
  }
  if(length(fungi) ==1){
    WHERE <- paste0(WHERE," AND Fungi = '",fungi, "'")
  }
  if(ts_DEGs_only){
    WHERE <- paste0(WHERE," AND GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")
  }

  my_query = paste0(SELECT, " FROM timeseries ",WHERE)

  exp <- db_query(my_query)


  ##add confidence intervals
  if(include_conf_intervals){
    confidence_intervals <- t(apply(exp[,c('rep1','rep2','rep3')], 1, calculate_confidence_interval))
    exp$'LCL' <- confidence_intervals[,1]
    exp$'UCL' <- confidence_intervals[,2]



  }

  return(exp)


}


get_timeseries_exp(ts_DEGs_only = TRUE)


get_GO_genes <- function(go_id='GO:0009723',
                             subset_lettuce_ids=NULL,ts_degs_only = FALSE){




  At_genes_w_GO <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                        keys=go_id,
                        columns=c("TAIR"),
                        keytype="GO")

  if(nrow(At_genes_w_GO)<2){return(NULL)}


  WHERE <- paste('WHERE AGI IN', get_genes_query(At_genes_w_GO$TAIR))

  if(ts_degs_only){
    WHERE <- paste0(WHERE, " AND GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")

  }else if(!is.null(subset_lettuce_ids)){
    let_genes_query <- get_genes_query(subset_lettuce_ids)
    if(!is.null(let_genes_query)){
      WHERE <- paste0(WHERE, " AND GeneID IN ", let_genes_query)
    }
  }

  query <- paste0("SELECT GeneID,",
                  "CASE WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN GeneID ",
                  "ELSE 'Ls' || At_ShortName END AS name ",
                  "FROM gene_annotations ", WHERE)


  Ls_genes_w_GO <- db_query(query)




  return(Ls_genes_w_GO)



}


get_genes_w_domain <- function(domain_id='PF00067',
                               domain_desc=NULL,
                               subset_lettuce_ids=NULL,
                               ts_degs_only = FALSE,
                               include_gene_name = TRUE){


  if( !is.null(domain_desc)){
    WHERE= paste("WHERE LOWER(da.description) LIKE",
                  sprintf("'%%%s%%'", tolower(domain_desc)))
  }else{
      WHERE = paste0("WHERE da.id =", "'",domain_id,"'")

  }

  if(ts_degs_only){
    WHERE <- paste0(WHERE, " AND da.GeneID IN (SELECT GeneID FROM timeseries_overlap_degs)")
  }else if(!is.null(subset_lettuce_ids)){
    let_genes_query <- get_genes_query(subset_lettuce_ids)
    if(!is.null(let_genes_query)){
      WHERE <- paste0(WHERE, " AND da.GeneID IN ", let_genes_query)

    }
  }

  if(include_gene_name){
    query <-paste("SELECT DISTINCT da.GeneID,",
           "CASE WHEN (ga.At_ShortName IS NULL OR ga.At_ShortName LIKE 'AT%G%') THEN da.GeneID",
           "ELSE 'Ls' || ga.At_ShortName END AS name",
           "FROM domain_annotations AS da",
           "LEFT JOIN gene_annotations AS ga ON da.GeneID = ga.GeneID",
           WHERE)

  }else{
    query <- paste0("SELECT DISTINCT GeneID FROM domain_annotations AS da ", WHERE)
    }

  print(query)




  Ls_genes_w_domain <- db_query(query)


  return(Ls_genes_w_domain)


}

get_orthologs <- function(AGI= 'AT3G10500', shortname=NULL,
                          subset_lettuce_ids=NULL,ts_degs_only = FALSE){



  if(!is.null(shortname)){

    AGI <- tair_code <- tryCatch({
      AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                                       keys = shortname,
                                       keytype = "SYMBOL",
                                       columns = "TAIR")$TAIR
    }, error = function(e) {
      # If an error occurs, return an empty data frame
      message("Error occurred while querying the TAIR database.")
      return('')
    })
  }
  AGI = toupper(AGI)
  if(!grepl('AT[0-9]G[0-9]+', AGI)){
    print('invalid AGI')
    return(NULL)
  }

  query <- paste0("SELECT GeneID FROM gene_annotations WHERE AGI ='",
                  AGI,"'")

  Ls_genes_w_domain <- db_query(query)




  return(Ls_genes_w_domain)




}


# get_GO_orthologs(go_id = 'GO:0009723', ts_degs_only=TRUE)
#test <- get_genes_w_domain(domain_desc = 'AP2', ts_degs_only=TRUE)
#get_orthologs(shortname = 'LBD29')


