library(tidyverse)

app_path <- 'G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/'
source(file.path(app_path,'sql_db/db_connect.R'))



allowed_bot <- c('bot','botrytis','b.cinerea','bc','bcin')
allowed_sclero <- c('scl','sclero','sclerotinia','ss','ssclerotiorum','ssclerotiroum')
allowed_fungi <- c(allowed_bot,allowed_sclero)


get_genes_query <- function(input_genes){
  valid_genes <- input_genes[grepl('Lsat_1_v5_gn_\\d',input_genes)]
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

get_timeseries_exp <- function(genes=c('Lsat_1_v5_gn_5_78841','Lsat_1_v5_gn_9_69201','Lsat_1_v5_gn_1_56161'),
                               fungi=c('Bc','Ss','lkjhjk'),
                               include_raw_points = TRUE,
                               include_mock = TRUE){

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

  my_query = paste0(SELECT, " FROM timeseries ",WHERE)

  exp <- db_query(my_query)


  ##add confidence intervals





}



get_GO_genes <- function(go_id='GO:0009723',
                             subset_lettuce_ids=NULL,ts_degs_only = FALSE){

  conn <- db_connect()



  At_genes_w_GO <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
                        keys=go_id,
                        columns=c("TAIR"),
                        keytype="GO")

  if(ts_degs_only){
    ts_degs <- dbGetQuery(conn= conn,statement = 'SELECT GeneID FROM timeseries_overlap_degs' )
    query = paste0("SELECT GeneID FROM gene_annotations WHERE AGI IN ",
                   paste0("('", paste(unique(At_genes_w_GO$TAIR), collapse = "', '"), "')"),
                   'AND GeneID IN', paste0("('", paste(ts_degs$GeneID, collapse = "', '"), "')"))


  }else if(is.null(subset_lettuce_ids)){
    query = paste0("SELECT GeneID FROM gene_annotations WHERE AGI IN ",
                   paste0("('", paste(unique(At_genes_w_GO$TAIR), collapse = "', '"), "')"))
  }else{

    query = paste0("SELECT GeneID FROM gene_annotations WHERE AGI IN ",
                   paste0("('", paste(unique(At_genes_w_GO$TAIR), collapse = "', '"), "')"),
                   'AND GeneID IN', paste0("('", paste(subset_lettuce_ids, collapse = "', '"), "')"))

  }
  Ls_genes_w_GO <- dbGetQuery(conn= conn,statement = query )

  db_disconnect(conn)



  return(Ls_genes_w_GO)



}


get_genes_w_domain <- function(domain_id='PF00067', domain_desc=NULL,
                               subset_lettuce_ids=NULL,ts_degs_only = FALSE){

  conn <- db_connect()

  if( !is.null(domain_desc)){
    domain_WHERE= paste(" LOWER(description) LIKE",
                  sprintf("'%%%s%%'", tolower(domain_desc)))
  }else{
    domain_WHERE = paste0(" id =", "'",domain_id,"'")

  }

  if(ts_degs_only){
    ts_degs <- dbGetQuery(conn= conn,statement = 'SELECT GeneID FROM timeseries_overlap_degs' )
    geneID_WHERE = paste0(' AND GeneID IN', paste0("('", paste(ts_degs$GeneID, collapse = "', '"), "')"))
  }else if(!is.null(subset_lettuce_ids)){
    paste0(' AND GeneID IN', paste0("('", paste0(subset_lettuce_ids, collapse = "', '"), "')"))

  }else{
    geneID_WHERE = ""
  }

  query <- paste0("SELECT DISTINCT GeneID FROM domain_annotations WHERE", domain_WHERE, geneID_WHERE)

  #print(query)



  Ls_genes_w_domain <- dbGetQuery(conn= conn,statement = query )

  db_disconnect(conn)

  return(Ls_genes_w_domain)


}

get_orthologs <- function(AGI= 'AT3G10500', shortname=NULL,
                          subset_lettuce_ids=NULL,ts_degs_only = FALSE){


  conn <- db_connect()

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

  Ls_genes_w_domain <- dbGetQuery(conn= conn,statement = query )




  db_disconnect(conn)
  return(Ls_genes_w_domain)




}


# get_GO_orthologs(go_id = 'GO:0009723', ts_degs_only=TRUE)
# test <- get_genes_w_domain(domain_desc = 'ethylene', ts_degs_only=FALSE)
#get_orthologs(shortname = 'LBD29')


