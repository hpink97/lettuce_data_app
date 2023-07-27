library(DBI)
library(RSQLite)

sqlite_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/sql_db/lettuce_data.sqlite"


db_connect <- function(db_path = sqlite_path){
  db_connection <- dbConnect(RSQLite::SQLite(), db_path)

  return(db_connection)
}

db_disconnect <- function(db_connection){
  dbDisconnect(db_connection)
}

db_query <- function(query=''){
  conn <- db_connect()
  df <- dbGetQuery(conn, query)
  db_disconnect(conn)

  return(df)
}

db_print_head <- function(table,topn=5){
  query <- paste('SELECT * FROM',table, 'LIMIT',topn)
  df <- db_query(query)
  return(df)
}

db_print_tables <- function(db_connection){
  conn <- db_connect()
  print(dbListTables(conn))
  db_disconnect(conn)
}



