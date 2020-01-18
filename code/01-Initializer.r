#This script takes a user specified Yaml file and creates a markdown file from it
#To run: Rscript Run.r <Configuration YAML FILE> <

library(yaml)
args <- commandArgs(trailingOnly=TRUE)
setwd( "~/AD_TargetRank/")

#READ IN Config
#config<-'configs/TestConfig.yaml'
configuration <- args[1]
config <- read_yaml(configuration)

#CREATE the RMD to RUN
system( paste0( "sed 's/YAML_FILE/", sub("/", "\\\\/",config$name), "/g' ~/AD_TargetRank/code/02-Master.Rmd > ~/AD_TargetRank/code/03-Run.Rmd"))

setwd("~/AD_TargetRank/")
source("~/AD_TargetRank/utilityFunctions/knitfile2synapseClient.R")
source("~/AD_TargetRank/utilityFunctions/hook_synapseMdSyntax_plot.R")
createAndKnitToFolderEntityClient(file = "~/AD_TargetRank/code/03-Run.Rmd",
                                  parentId =config$parentID,
                                  folderName = config$runID)
