##Coexpression Analysis module:
#devtools::install_github("blogsdon/spike/spike/")#, repos = NULL, type = "source")
source("~/AD_TargetRank/utilityFunctions/Parallel_vbsrBootstrap.R")
# Load Libraries and se working directory 
library(pheatmap)
library(parallel)
library(doParallel)
library(spike)
library(reshape2)
setwd( config$filedir )

# Set Github Provenance links
thisRepo <- githubr::getRepo(repository = "jgockley62/AD_TargetRank", ref="branch", refName='master')
GL_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath=config$genelistfile)
CF_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath=config$name)
I_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/01-Initializer.r')
M_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/02-Master.Rmd')

CODE <- syn_temp$store(synapseclient$Folder(name = 'figures', parentId = RunParent))

#SynIDs of expression data to pull from
ExpressionDS <- c('syn21291908','syn21292041','syn21285564','syn21285564','syn21285564','syn21285564','syn21291908')
names( ExpressionDS ) <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

#Study ID Translator
Study <- c( 'RosMap', 'Mayo', 'Mayo', 'MSBB', 'MSBB', 'MSBB', 'MSBB')
names(Study) <- c('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

Syn <- list('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

#Useable Chrs for genelist ENSG  translations and wiki analysis preamble
row.names(Trans) <- Trans$ensembl_gene_id
Trans <- Trans[ Trans$chromosome_name %in% c(1:22,"X","Y"), ]

writeLines(paste0('For each tissue specified in ', 
                  paste( config$name, collapse = ', '), 
                  " the coexpression of genes specified in ", 
                  config$genelistfile, 
                  " is given. This is measured as a percentage of 1000 pemutations in which the gene was a significant predictor of expression of its comparison gene by a linear model."
)
)
# Make objects for list of Genes Not Expressed In Tissue and dataframes to plot
Missing <- list()
Plots <- list()

FinalTis <- data.frame()

library(parallel)
library(doParallel)
cores <- detectCores()-2 
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)


for( Tissue in config$tissue ){
  
  parentId = 'syn21532474';
  activityName = paste0( Tissue,' Partial Correlations');
  activityDescription = paste0( Tissue,'Pairwise Partial Correlation Calculation');
  thisFileName <- 'PartialCorelations.R'
  # Github link
  thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
  thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/prepscripts/',thisFileName))
  
  activityName = 'Covariate and Diagnosis Regression';
  activityDescription = 'Covariate analysis and Regression of aligned effective counts with GRCh38 with CQN normalisation (IFG, STG, FP, PHG)';
  
  CODE <- syn_temp$store(synapseclient$Folder(name = "Reference Data", parentId = parentId))
  
  #Set Used SynIDs For Provenance
  # Set annotations
  all.annotations = list(
    dataType = 'RNA-Seq',
    dataSubType = 'Partial Correlations',
    summaryLevel = 'gene',
    assay	 = 'RNAseq',
    tissueTypeAbrv	= Tissue, 
    study = paste0( Study[Tissue] ), 
    organism = 'HomoSapiens',
    consortium	= 'AMPAD',
    normalizationStatus	= TRUE,
    normalizationType	= 'CQN',
    rnaquantification = 'RSEM',
    genomeAssemblyID = 'GRCh37'
  )
  
  
  #Tissue<-'CBE'
  Syns_Used <- as.character(ExpressionDS[Tissue])
  #Tissue <- 'STG'
  
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
  exp <- exp[rowSums(is.na(exp))<ncol(exp), ]
  
  foo <- bcv::impute.svd( t(exp) )
  Exp <- foo$x
  row.names(Exp) <- row.names(t(exp)) 
  colnames(Exp) <- colnames(t(exp)) 
  
  #Record Genes missing in tissue of interest
  if( as.numeric(table(row.names(Trans) %in% colnames(x) )["TRUE"] ) == dim(Trans)[1] ){
    Missing[[ Tissue ]] <- paste0("All Listed Genes Expressed in ", Tissue)
  }else{
    Missing[[ Tissue ]] <- paste0( "Queried Genes missing from ", Tissue, " : ",
                                   paste(c(Trans[ (row.names(Trans) %in% colnames(x) ) == F, ]$hgnc_symbol), collapse = ', ')
    )
  }
  #Build output matrix for partial correlation
  
  #Pairwise spearman correlation of genes
  XCor <- cor(x, method = "spearman")
  
  #Prep for partial correlation detection
  Final <- data.frame()
  
  source("~/AD_TargetRank/utilityFunctions/Parallel_vbsrBootstrap.R")
  #Run Partial correlations for each gene
  
  RUNNe <- function( i=i, x=x ){ 
    source("~/AD_TargetRank/utilityFunctions/Parallel_vbsrBootstrap.R")
    OBS <- i
    #OBS <- 'ENSG00000227232'
    y <- as.matrix(x[,OBS]) 
    colnames(y) <- i
    X <- x[,(colnames(x) %in% OBS) == F ]
    att <- pvbsrBootstrap( y=y, x=X, nsamp=100, cores=1 )
    names( att ) <- gsub( "intercept", OBS, names(att) )
    att <- c( 'SeedGene' = OBS, att[ colnames(x) ]) 
    return( att[c( 'SeedGene', colnames(x) )] )
    #return(eval(parse( text = paste( c( OBS,att[colnames(x)]), sep='\t') )))
  }
  
  LIST <- as.matrix(colnames(x))
  
  mark <- Sys.time()
  foo <- t( parApply(cl, as.matrix( LIST[1:dim(LIST)[1],] ), 1, RUNNe, x) )
  Sys.time()-mark
  
  row.names( foo ) <-  foo[,1]
  foo <-  foo[, colnames(foo)[ (colnames(foo) %in% "SeedGene") == F ]]
  write.table( foo, file=paste0(Tissue,"_PartialCorMatrix.tsv"), row.names=T, col.names = T, quote=F, sep='\t' )
  DT <- melt(t(foo))
  colnames(DT) <- c( "Hit Gene", "Target Gene", "PartialCor") 
  write.table(DT, file=paste0(Tissue,"_PartialCorTable.tsv"), row.names=T, col.names = T, quote=F, sep='\t' )

  
  #Push to Synapse
  
  # Store Matrix
  ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0(Tissue,"_PartialCorMatrix.tsv"), name = paste0(Tissue, ' Pairwise Partial correlation Matrix'), parentId=CODE$properties$id ), used = Syns_Used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
  all.annotations$dataSubType = 'PartialCorrelationsForADTargetRank'
  syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)
  # Store Tables
  ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0(Tissue,"_PartialCorTable.tsv"), name = paste0(Tissue, ' Pairwise Partial correlation Table'), parentId=CODE$properties$id ), used = Syns_Used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
  all.annotations$dataSubType = 'PartialCorrelationsForADTargetRank'
  syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)
} 

stopCluster(cl)