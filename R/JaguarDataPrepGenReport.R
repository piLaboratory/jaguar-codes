#'
#' #          **Jaguar Data Preparation Generator Report**
#' 
#' #### *Alan E. de Barros, Bernardo Niebuhr, Vanesa Bejarano, Julia Oshima,Claudia Kanda, Milton Ribeiro, Ronaldo Morato,Paulo Prado*
#' date: "March, 08 2019"
#' 
#' #### *Preamble*
#' For a fresh start, clean everything in working memory
rm(list= ls())                                                 
#' 
#' Install and load packages
if(!require(install.load)) install.packages('install.load'); library(install.load)
#'
#' #### Produce an Rmd and html outputs using ezknitr
install.load::install_load("knitr", "ezknitr") # To render documents 
ezspin(file = "JaguarDataPrep.R", out_dir = "reports",
       params = list("DATASET_NAME" = "jaguar.dat"), 
       keep_html = TRUE, keep_rmd = TRUE)

open_output_dir()

