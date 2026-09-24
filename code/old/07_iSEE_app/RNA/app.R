library(SpatialExperiment)
library(iSEE)
library(shiny)
library(scuttle)
library(qs2)
# library(here)

#   For interactive testing at JHPCE
# setwd(here('code', '07_iSEE_app', 'RNA'))

source("initial.R")

sce = qs_read('sce_RNA_iSEE.qs2')

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce)
# 2.32 GB

iSEE(
    sce,
    appTitle = "Habenula Atlas Multiome RNA Data",
    initial = initial
)
