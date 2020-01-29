library(parallel)
library(doParallel)
reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

setwd("~/AD_TargetRank/")
source("~/AD_TargetRank/utilityFunctions/knitfile2synapseClient.R")
source("~/AD_TargetRank/utilityFunctions/hook_synapseMdSyntax_plot.R")

ExpressionDS <- c('syn21291908','syn21292041','syn21285564','syn21285564','syn21285564','syn21285564','syn21291908')
names( ExpressionDS ) <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

#Study ID Translator
Study <- c( 'RosMap', 'Mayo', 'Mayo', 'MSBB', 'MSBB', 'MSBB', 'MSBB')
names(Study) <- c('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

for( Tissue in c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX') ){
  
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
  
  mark <- Sys.time()
  LIST <- combn(1:length(colnames(x)),2)
  Sys.time()-mark
  
  rm(foo)
  rm(Exp)
  rm(exp)
  
  mark <- Sys.time()
  foo <- t(parApply(cl,LIST,2,RUNNEr,x))
  Sys.time()-mark
  
  rm(LIST)
  
  colnames(FOO) <- c("Target_Gene", "Hit_Gene", "Spearman_Correlation", "Spearman_pVal", "FDR_CorPval")
  FOO$Spearman_Correlation <- as.numeric(as.character( FOO$Spearman_Correlation ))
  FOO$Spearman_pVal <- as.numeric(as.character( FOO$Spearman_pVal ))
  FOO$FDR_CorPval <- as.numeric(as.character( FOO$FDR_CorPval ))
  FOO$Hit_Tissue <- rep( Tissue,dim(FOO)[1] )
  FOO$Hit_Data_Set <- rep( as.character(Study[Tissue]),dim(FOO)[1] )
  #Write A to File
  
  #Write B to File - No header
  #_#foo[,c(2,1,3,4,5)])
  
  #System Call Concatenate A and B
  
  #Erase A and B
  
  #Clear R Memory...
  rm(foo)
  
  #Push files to synapse:
  thisRepo <- githubr::getRepo(repository = "jgockley62/AD_TargetRank", ref="branch", refName='master' )
  thisFile  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/prepscripts/07_1-SpearmanCor.R' )
  
  CODE <- syn_temp$store(synapseclient$Folder(name = 'Reference Data', parentId = 'syn21532474') )
  
  # Set annotations for synapse objects ( figures and tables )
  all.annotations = list(
    dataType = 'Coexpression',
    dataSubType = 'spearman correlations',
    summaryLevel = 'gene',
    assay	 = 'RNAseq',
    tissueTypeAbrv	= config$tissue, 
    study = c('Mayo', "MSBB", 'RosMap'), 
    organism = 'HomoSapiens',
    consortium	= 'AMPAD',
    normalizationStatus	= TRUE,
    normalizationType	= 'CQN',
    rnaquantification = 'RSEM',
    genomeAssemblyID = 'GRCh38'
  )
  
  #syns_Used <- eval(parse(text=paste0( 'c(', paste0('Syn$', config$tissue, collapse =", " ), ')')))
  
  activityName = paste0( Tissue, 'Coexpression heatmap')
  activityDescription = paste0( Tissue, 'Coexpression heatmap for genelist: ', config$genelistfile)
  
  CODE <- syn_temp$store(synapseclient$Folder(name = 'data', parentId = RunParent))
  
  write.table(OUT_Cor, file = paste0('./', Tissue, '_OUT_Cor.tsv'), row.names = T, col.names = T, quote=F, sep='\t' )
  write.table(OUT_pVal, file = paste0('./', Tissue, '_OUT_pval.tsv'), row.names = T, col.names = T, quote=F, sep='\t')
  
  ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0('runs/',config$runID,'PartialCorrelationPerms.tsv'), 
                                                     name = paste0('PartialCorrelationPerms.tsv'), 
                                                     parentId=CODE$properties$id ), 
                                 used = Syns_Used,
                                 activityName = activityName, 
                                 executed = list(thisFile),
                                 activityDescription = activityDescription)
  
  all.annotations$dataSubType = 'Spearman Correlations'
  syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)
  
}

stopCluster(cl)