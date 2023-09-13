library(shiny)
library(shinythemes)
library(shinycssloaders)
library(DT)

rm(list=ls())


source('sql_db/db_connect.R')
source('plotting_funcs/helper_funcs.R')
source('plotting_funcs/plot_divset_cor.R')
source('plotting_funcs/identify_subset_hubs.R')
source('plotting_funcs/plot_timeseries_expr.R')
source('ui.R')
source('server.R')




shinyApp(ui = ui, server = server)
