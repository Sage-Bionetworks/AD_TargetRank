#This script takes a user specified Yaml file and creates a markdown file from it
#To run: Rscript code/01-Initializer.r <Configuration YAML FILE>
#eg. Rscript code/01-Initializer.r configs/TestConfig.yaml
#eg. Rscript code/01-Initializer.r configs/M62.yaml
#eg. Rscript code/01-Initializer.r configs/pQTL.yaml
#eg. Rscript code/01-Initializer.r configs/Retromer.yaml
#eg. Rscript code/01-Initializer.r configs/Moesin.yaml
#eg. Rscript code/01-Initializer.r configs/EmoryTargets.yaml

library(yaml)
args <- commandArgs(trailingOnly=TRUE)

#READ IN Config
#configuration<-'configs/EmoryTargets.yaml'
configuration <- args[1]
config <- read_yaml(configuration)

setwd( config$filedir )
system( paste0('mkdir -p runs/', config$runname, '/figures' ) )
system( paste0('mkdir -p runs/', config$runname, '/tables' ) )

#CREATE the RMD to RUN
system( paste0( "sed 's/YAML_FILE/", sub("/", "\\\\/",config$name), "/g' ~/AD_TargetRank/code/02-Master.Rmd > ~/AD_TargetRank/code/03-Run.Rmd"))

reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

setwd("~/AD_TargetRank/")
source("~/AD_TargetRank/utilityFunctions/knitfile2synapseClient.R")
source("~/AD_TargetRank/utilityFunctions/hook_synapseMdSyntax_plot.R")
createAndKnitToFolderEntityClient(file = "~/AD_TargetRank/code/03-Run.Rmd",
                                  parentId =config$parentID,
                                  folderName = config$runID)
