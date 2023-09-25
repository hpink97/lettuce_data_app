# Load necessary libraries
library(DBI)
library(RSQLite)

# Set the path to the SQLite database file
sqlite_path <- "/mnt/shiny/lettuce_transcriptomics/lettuce_data.sqlite"

# Function to establish a connection to the SQLite database
db_connect <- function(db_path = sqlite_path){
  db_connection <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  return(db_connection)
}

# Function to disconnect from the SQLite database
db_disconnect <- function(db_connection){
  DBI::dbDisconnect(db_connection)
}

# Function to execute a custom SQL query and return results as a data frame
db_query <- function(query=''){
  conn <- db_connect()
  df <- dbGetQuery(conn, query)
  db_disconnect(conn)

  return(df)
}

# Function to print the first few rows of a specified table in the database
db_print_head <- function(table, topn=5){
  query <- paste('SELECT * FROM', table, 'LIMIT', topn)
  df <- db_query(query)
  return(df)
}

# Function to print the names of all tables in the database
db_print_tables <- function(db_connection){
  conn <- db_connect()
  print(dbListTables(conn))
  db_disconnect(conn)
}
