library(parallel)
library(doParallel)
library(data.table)
reticulate::use_python("/usr/bin/python", required = TRUE)
#reticulate::use_python("../miniconda2/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

setwd("~/AD_TargetRank/")
source("~/AD_TargetRank/utilityFunctions/knitfile2synapseClient.R")
source("~/AD_TargetRank/utilityFunctions/hook_synapseMdSyntax_plot.R")
parentId = 'syn21532474'

ExpressionDS <- c('syn21291908','syn21292041','syn21285564','syn21285564','syn21285564','syn21285564','syn21291908')
names( ExpressionDS ) <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

#Study ID Translator
Study <- c( 'RosMap', 'Mayo', 'Mayo', 'MSBB', 'MSBB', 'MSBB', 'MSBB')
names(Study) <- c('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

for( Tissue in c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX') ){
  #Tissue<-'DLPFC'
  Syns_Used <- NULL
  
  #Load expression for tissue
  exp <- read.table(syn_temp$get(as.character(ExpressionDS[Tissue]))$path, header =T, sep ='\t', row.names=1)
  colnames( exp ) <- gsub( "X", "", colnames( exp ) )
  Syns_Used <- c(Syns_Used, as.character(ExpressionDS[Tissue] ) )
  #Seperate exp by tissue
  if( Tissue == 'DLPFC'){
  }else{
    if( Tissue == 'TCX' | Tissue == 'CBE' ){
      if( Tissue == 'CBE' ){
        slec <- 'CER'
      }else{ slec <- 'TCX' }
      colnames( exp ) <- gsub( "TC", "TCX", colnames( exp ) )
      exp <- exp[ , grepl( slec, colnames(exp)) ]
    }else{
      if( Tissue %in% c('FP', 'IFG', 'PHG', 'STG') ){
        Syns_Used <- c( Syns_Used, 'syn21285520' )
        Meta <- read.table( syn_temp$get('syn21285520')$path, header =T, sep ='\t', row.names=1 )
        exp <- exp[ colnames(exp) %in% row.names(Meta[grepl( Tissue, Meta$Tissue.Diagnosis),])]
      }else{
        stop(paste0("ERROR: SOURCE=Config.yaml Issue=Tissue: ", Tissue," is improper must be one of: CBE, DLPFC, FP, IFG, PHG, STG, TCX"))
      }
    }
  }
  
  #Impute svalues for given gene-patient NA values
  foo <- bcv::impute.svd( t(exp) )
  Exp <- foo$x
  row.names(Exp) <- row.names(t(exp)) 
  colnames(Exp) <- colnames(t(exp)) 
  
  x <- Exp
  
  RUNNEr <- function(n, x){
    test <- cor.test( x[ ,n[1] ], x[ ,n[2] ], method = "spearman" )
    return( c( colnames(x)[n[1]],
               colnames(x)[n[2]], 
               as.numeric(test$estimate), 
               as.numeric(test$p.value), 
               p.adjust( as.numeric(test$p.value), method = "fdr", n=dim(x)[2]) )
          )
  }
  
  #mark <- Sys.time()
  LIST <- combn(1:length(colnames(x)),2)
  #Sys.time()-mark
  
  rm(foo)
  rm(Exp)
  rm(exp)
  
  cl <- makeCluster(detectCores()-12)
  registerDoParallel(cl)
  #mark <- Sys.time()
  foo <- t(parApply(cl,LIST,2,RUNNEr,x))
  #Sys.time()-mark
  stopCluster(cl)
  
  rm(LIST)
  
  #Write A to File
  write.table(foo, "Attemp.txt", quote=F, row.names=F, col.names=F )
  rm(foo)
  rm(x)
  
  system( paste0('bash code/prepscripts/SpearmanAnnotator.sh ', Tissue, ' ', Study[Tissue] ) )
  
  #Push files to synapse:
  thisRepo <- githubr::getRepo(repository = "jgockley62/AD_TargetRank", ref="branch", refName='master' )
  thisFile  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/prepscripts/07_1-SpearmanCor.R' )
  thisFileAlso  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/prepscripts/SpearmanAnnotator.sh' )
  
  CODE <- syn_temp$store(synapseclient$Folder(name = 'Reference Data', parentId = 'syn21532474') )
  
  # Set annotations for synapse objects ( figures and tables )
  all.annotations = list(
    dataType = 'Coexpression',
    dataSubType = 'spearman correlations',
    summaryLevel = 'gene',
    assay	 = 'RNAseq',
    tissueTypeAbrv	= Tissue, 
    study = c('Mayo', "MSBB", 'RosMap'), 
    organism = 'HomoSapiens',
    consortium	= 'AMPAD',
    normalizationStatus	= TRUE,
    normalizationType	= 'CQN',
    rnaquantification = 'RSEM',
    genomeAssemblyID = 'GRCh38'
  )
  
  activityName = paste0( Tissue, 'Spearman Correlations')
  activityDescription = paste0( Tissue, 'Spearman Correlations and pValues')
  
  CODE <- syn_temp$store(synapseclient$Folder(name = 'Reference Data', parentId = parentId))
  
  ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0(Tissue,'_SpearmanCor.txt'), 
                                                     name = paste0(Tissue,' Spearman Correlations'), 
                                                     parentId=CODE$properties$id ), 
                                                     used = Syns_Used,
                                                     activityName = activityName, 
                                                     executed = list(thisFile, thisFileAlso),
                                                     activityDescription = activityDescription)
  
  all.annotations$dataSubType = 'Spearman Correlations'
  syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)
  
  #Remove Latent Files
  file.remove('Attemp.txt')
  file.remove(paste0(Tissue, '_SpearmanCor.txt'))
}


