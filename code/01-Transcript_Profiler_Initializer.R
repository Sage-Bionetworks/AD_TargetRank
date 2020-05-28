#This script takes a user specified Yaml file and creates a markdown file from it
#To run: Rscript code/01-Transcript_Profiler_Initializer.R <Configuration YAML FILE>
#eg. Rscript code/01-Initializer.r configs/Initial_49_Targets.yaml

library(yaml)
args <- commandArgs(trailingOnly=TRUE)

#READ IN Config
#configuration<-'configs/Initial_49_Targets.yaml'
configuration <- args[1]
config <- read_yaml(configuration)

setwd( config$filedir )

#Load Synapse and our DE Data:
reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

#DE Results
MSBB_DE <- data.table::fread( syn_temp$get('syn22092611')$path )
colnames(Mayo_DE)[ colnames(Mayo_DE) == 'Tissue.ref' ] <- 'Tissue'
Mayo_DE <- data.table::fread( syn_temp$get('syn22093747')$path )
Rosmap_DE <- data.table::fread( syn_temp$get('syn22097286')$path )
colnames(Rosmap_DE)[ colnames(Rosmap_DE) == 'Region' ] <- 'Tissue'

#Normalized Expression
MSBB_EXP <- as.data.frame(data.table::fread( syn_temp$get('syn22092607')$path ))
row.names(MSBB_EXP) <- MSBB_EXP$ensembl_gene_id
MSBB_EXP <- MSBB_EXP[,colnames(MSBB_EXP) != 'ensembl_gene_id']

Mayo_EXP <- as.data.frame( data.table::fread( syn_temp$get('syn22093746')$path ) )
row.names(Mayo_EXP) <- Mayo_EXP$ensembl_gene_id
Mayo_EXP <- Mayo_EXP[,colnames(Mayo_EXP) != 'ensembl_gene_id']

Rosmap_EXP <- as.data.frame( data.table::fread( syn_temp$get('syn22097285')$path ) )
row.names(Rosmap_EXP) <- Rosmap_EXP$ensembl_gene_id
Rosmap_EXP <- Rosmap_EXP[,colnames(Rosmap_EXP) != 'ensembl_gene_id']

#Covariates
MSBB_Cov <- data.table::fread( syn_temp$get('syn22092576')$path )
row.names(MSBB_Cov) <- MSBB_Cov$SampleID
Mayo_Cov <- data.table::fread( syn_temp$get('syn22093739')$path )
row.names(Mayo_Cov) <- Mayo_Cov$SampleID
Rosmap_Cov <- data.table::fread( syn_temp$get('syn22097275')$path )
row.names(Rosmap_Cov) <- Rosmap_Cov$SampleID

#Get Transcripts for Gene List
Genes <- as.vector( read.table(config$genelistfile)$V1 )

library("biomaRt")
gene_id <- 'ensembl_transcript_id'

if( config$genelisttype == "GNAME" ){
  filters <- 'hgnc_symbol'
  attrs <- c(filters, "ensembl_transcript_id", "ensembl_gene_id")
}else{
  if( config$genelisttype == "ENSG" ){
    filters <- 'ensembl_gene_id'
    attrs <- c(filters, "ensembl_transcript_id", "hgnc_symbol")
  }else{
    warning( "GENE NAME TYPE NOT SPECIFIED CORRECTLY Please Use ENSG or GNAME")
  }
}

host <- 'ensembl.org'
organism <- 'hsa'
attrs <- c(filters, "ensembl_transcript_id", "ensembl_gene_id")
ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)

ds <- biomaRt::listDatasets(ensembl)[, "dataset"]
ds <- grep(paste0("^", organism), ds, value = TRUE)
if (length(ds) == 0) {
  stop(paste("Mart not found for:", organism))
} else if (length(ds) > 1) {
  message("Found several marts")
  sapply(ds, function(d) message(paste(which(ds == d), d, sep = ": ")))
  n <- readline(paste0("Choose mart (1-", length(ds), ") : "))
  ds <- ds[as.integer(n)]
}
ensembl <- biomaRt::useDataset(ds, mart = ensembl)

Trans <- getBM(filters = 'hgnc_symbol',
               attributes = attrs,
               values = Genes,
               mart = ensembl )


#Filter the Larget tables and write out the tiny tidy data
##Make into temp files and store in temp file directory
features <- Trans$ensembl_transcript_id

#DE Results - store
MSBB_DE <- MSBB_DE[ MSBB_DE$ensembl_gene_id %in% features, ] 
write.table( as.data.frame(MSBB_DE), file = 'Transcript_Profiles/Temp_Files/MSBB_DE.tsv', col.names=T, sep ='\t', row.names = F, quote = F )
Mayo_DE <- Mayo_DE[ Mayo_DE$ensembl_gene_id %in% features, ] 
write.table( as.data.frame(Mayo_DE), file = 'Transcript_Profiles/Temp_Files/Mayo_DE.tsv', col.names=T, sep ='\t', row.names = F, quote = F )
Rosmap_DE <- Rosmap_DE[ Rosmap_DE$ensembl_gene_id %in% features, ] 
write.table( as.data.frame(Rosmap_DE), file = 'Transcript_Profiles/Temp_Files/RosMap_DE.tsv', col.names=T, sep ='\t', row.names = F, quote = F )

#Normalized Expression - store
MSBB_EXP <- MSBB_EXP[ row.names(MSBB_EXP)[ row.names(MSBB_EXP) %in% features ], ]
write.table( as.data.frame(MSBB_EXP), file = 'Transcript_Profiles/Temp_Files/MSBB_Exp.tsv', col.names=T, sep ='\t', row.names = T, quote = F )
Mayo_EXP <- Mayo_EXP[ row.names(Mayo_EXP)[ row.names(Mayo_EXP) %in% features ], ]
write.table( as.data.frame(Mayo_EXP), file = 'Transcript_Profiles/Temp_Files/Mayo_Exp.tsv', col.names=T, sep ='\t', row.names = T, quote = F )
Rosmap_EXP <- Rosmap_EXP[ row.names(Rosmap_EXP)[ row.names(Rosmap_EXP) %in% features ], ]
write.table( as.data.frame(Rosmap_EXP), file = 'Transcript_Profiles/Temp_Files/RosMap_Exp.tsv', col.names=T, sep ='\t', row.names = T, quote = F )

#Write Trans to temp file
write.table( as.data.frame(Trans), file = 'Transcript_Profiles/Temp_Files/Translate.tsv', col.names=T, sep ='\t', row.names = F, quote = F )

## Create a Gene Specific Markdown file for each gene and execute it. Then push to synapse
#CREATE the RMD to RUN
genes <- names( table(Trans$hgnc_symbol) )
#name <- 'ACE'

if( dir.exists( paste0( '~/AD_TargetRank/Transcript_Profiles/', config$runname )) ){
}else{
  system(paste0( 'mkdir ~/AD_TargetRank/Transcript_Profiles/', config$runname ))
}

parentId = config$parentID;
activityName = 'Transcript Analysis';
activityDescription = 'A transcript analysis of all transcripts across an AMP-AD gene';
thisFileName <- '02-Transcript_Master.Rmd'
thisFileName2 <- '01-Transcript_Profiler_Initializer.R'
# Github link
thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink( repository = thisRepo, repositoryPath=paste0( 'code/',thisFileName ))
thisFile2 <- githubr::getPermlink( repository = thisRepo, repositoryPath=paste0( 'code/',thisFileName2 ))

Syns_Used <- c("syn22092611", "syn22093747", "syn22097286", "syn22092607", "syn22093746", "syn22097285", "syn22092576", "syn22093739", "syn22097275")
CODE <- syn_temp$store(synapseclient$Folder(name = config$runname, parentId = parentId))

for( name in genes ){
  system( paste0( "sed 's/GENEX/", sub("/", "\\\\/", name), "/g' ", config$filedir, "/code/02-Transcript_Master.Rmd > ", config$filedir, "/code/03-", name,".Rmd"))
  eval(parse(text= paste0( 'rmarkdown::render(\'~/AD_TargetRank/code/03-', 
                           name,
                           '.Rmd\', \'pdf_document\', output_dir = \'~/AD_TargetRank/Transcript_Profiles/',
                           config$runname, '\' )' ) ))
  
  #Load to Synapse and scrub from files
  system( paste0( 'rm ~/AD_TargetRank/code/03-', name, '.Rmd' ) )
  system( paste0( 'rm -r ~/AD_TargetRank/Transcript_Profiles/Initial_49_Targets/03-', name, '_files/' ) )
  
  reticulate::use_python("/usr/bin/python", required = TRUE)
  synapseclient <- reticulate::import("synapseclient")
  syn_temp <- synapseclient$Synapse()
  syn_temp$login()
  
  ENRICH_OBJ <- syn_temp$store( synapseclient$File( path= paste0( '~/AD_TargetRank/Transcript_Profiles/', config$runname,'/03-', name,'.pdf'), name = paste0( 'Transcript analysis for ', name ), 
                                                     parentId=CODE$properties$id ), 
                                                     used = Syns_Used,
                                                     activityName = activityName, 
                                                     executed = c(thisFile,thisFile2),
                                                     activityDescription = activityDescription)
}


