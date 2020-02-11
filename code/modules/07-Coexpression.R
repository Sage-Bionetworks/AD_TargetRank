##Coexpression Analysis module:
#install.packages("~/Desktop/Programs/spike/spike/", repos = NULL, type = "source")
source("~/AD_TargetRank/utilityFunctions/Parallel_vbsrBootstrap.R")
# Load Libraries and se working directory 
library(pheatmap)
library(parallel)
library(doParallel)

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

for( Tissue in config$tissue ){
  
  Syns_Used <- NULL
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
  
  #x <- Exp[ , colnames(Exp) %in% row.names(Trans)]
  x <- Exp
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
  
  #mark<-Sys.time()
  #FOO <- psych::corr.test(x, method = "spearman")
  #Sys.time()-mark
  
  #Prep for partial correlation detection
  Final <- data.frame()
  core <- parallel::detectCores()-2 
  source("~/AD_TargetRank/utilityFunctions/Parallel_vbsrBootstrap.R")
  #Run Partial correlations for each gene
  RUNNe <- function( i=i, x=x){
    OBS <- i
    y <- as.matrix(x[,OBS]) 
    colnames(y) <- i
    X <- x[,(colnames(x) %in% OBS) == F ]
    att <- pvbsrBootstrap( y=y, x=X, nsamp=10, cores=1 )
    names( att["intercept"] ) <- eval( parse(text = OBS))
    return(att)
  }
  mark <- Sys.time()
  
  Sys.time()-mark
  
  for( i in row.names(Trans) ){
    OBS <- i
    #OBS <- 'ENSG00000084234'
    y <- as.matrix(x[,OBS]) 
    colnames(y) <- i
    X <- x[,(colnames(x) %in% OBS) == F ]
    
    
    mark <- Sys.time()
    att <- pvbsrBootstrap( y=y, x=X, nsamp=10, cores=1 )
    Sys.time()-mark
    
    mark <- Sys.time()
    att <- pvbsrBootstrap( y=y, x=X, nsamp=100, cores=1 )
    Sys.time()-mark
    
    mark <- Sys.time()
    att <- pvbsrBootstrap( y=y, x=X, nsamp=200, cores=1 )
    Sys.time()-mark
    
    mark <- Sys.time()
    att <- pvbsrBootstrap( y=y, x=X, nsamp=500, cores=1 )
    Sys.time()-mark
    
    mark <- Sys.time()
    att <- pvbsrBootstrap( y=y, x=X, nsamp=1000, cores=1 )
    Sys.time()-mark
    
    
    
    att <- pvbsrBootstrap( y=y, x=X, nsamp=100, cores=core )
    
    #Build Iteration based DF
    temp <- as.data.frame( cbind( Target_Gene = rep(i,length(att)), 
                Seed_Set = rep(config$runID,length(att)),
                Hit_Tissue = rep(Tissue,length(att)),
                Hit_Data_Set = rep( as.character(Study[Tissue]),length(att) ),
                Hit_Gene = names(att),
                Hit_Probability= as.numeric(att),
                Spearman_Correlation= rep( NA,length(att) )
              ))
    
    temp$Hit_Gene <- as.character(temp$Hit_Gene)
    temp[ temp$Hit_Gene == 'intercept',]$Hit_Gene <- i
    temp$Spearman_Correlation <- XCor[ i, temp$Hit_Gene ]
    
    #Add to Genelist based DF
    Final <- as.data.frame(rbind(temp,Final))
  }
  Sys.time() - mark
  
  SwapName <- GeneData( colnames(x) , 'ENSG' )
  #row.names(SwapName) <- SwapName$ensembl_gene_id
  # Replace blank gene names
  SwapName[SwapName$hgnc_symbol=="",]$hgnc_symbol <- SwapName[SwapName$hgnc_symbol=="",]$ensembl_gene_id
  SwapName <- SwapName[ !duplicated(SwapName$ensembl_gene_id), ]
  row.names(SwapName) <- SwapName$ensembl_gene_id
  
  SwapName$ensembl_gene_id <- as.character(SwapName$ensembl_gene_id)
  SwapName$hgnc_symbol <- as.character(SwapName$hgnc_symbol)
  SwapName$chromosome_name <- as.character(SwapName$chromosome_name)
  #SwapName$start_position <- as.numeric(as.character(SwapName$hgnc_symbol))
  #SwapName$end_position <- as.numeric(as.character(SwapName$hgnc_symbol))
  
  addon <- as.data.frame( matrix( NA, length(colnames(x)[ (colnames(x) %in% row.names(SwapName)) == F ]), 5) )
  addon[,1] <- colnames(x)[ (colnames(x) %in% row.names(SwapName)) == F ]
  addon[,2] <- colnames(x)[ (colnames(x) %in% row.names(SwapName)) == F ]
  row.names(addon) <- addon[,1]
  colnames(addon) <- colnames(SwapName)
  
  addon$chromosome_name <- as.character(addon$chromosome_name)
  addon$start_position <- as.integer(addon$start_position)
  addon$end_position <- as.integer(addon$end_position)
  
  SwapName <- as.data.frame( rbind(addon,SwapName) )
  #Rename the output partial correlation matrix:
  
  #row.names(Final) <- Trans[ row.names(Final), ]$hgnc_symbol
  Final$Target_GeneName <- SwapName[ as.character(Final$Target_Gene),]$hgnc_symbol
  Final$Hit_GeneName <- SwapName[ as.character(Final$Hit_Gene),]$hgnc_symbol
  #colnames(Final) <- Trans[ colnames(Final), ]$hgnc_symbol
  
  ORD <- c("Target_Gene", "Target_GeneName", 
    "Seed_Set", "Hit_Tissue", "Hit_Data_Set",
    "Hit_Gene", "Hit_GeneName",
    "Hit_Probability", "Spearman_Correlation"
  )
  
  Final <- Final[,ORD]
  FinalTis <- do.call( rbind, list(FinalTis, Final) )
  
  #Write Plot info to list, plot to file, and plot to wiki
  Ps <- data.frame( matrix( 0, length(table(Final$Target_Gene)), as.numeric(table(Final$Target_Gene)[1]) ))
  row.names( Ps ) <- names(table(Final$Target_GeneName))
  colnames( Ps ) <- names(table(Final$Hit_Gene))
  
  for( nom in row.names( Ps ) ){
    small <- Final[ Final$Target_GeneName == nom, ]
    row.names(small) <- as.character(small$Hit_Gene)
    Ps[ nom, ] <- small[ colnames(Ps), ]$Hit_Probability
  }
  
  #Exclude Hit Targets with zero probability:
  #Ps <- as.numeric(Ps)
  at <- as.data.frame(sapply(Ps, as.numeric))
  row.names(at) <- row.names(Ps)
  
  Ps = at[,colSums(at) > 0.4  ]
  
  
  Plots[[ Tissue ]] <- list(ParCor = t(Ps), main=Tissue)
  #pheatmap(t(Ps), main=Tissue, show_rownames=F)
  
  eval( parse( text= paste0( 'Syn$', Tissue, '<- Syns_Used' ) ))
}

for( Tissue in config$tissue ){
  eval( parse( text = paste0( 'pheatmap(Plots$',
                              Tissue,
                              '$ParCor, main=Plots$',
                              Tissue,
                              '$main, show_rownames=F)'
                         )))
} 
for( Tissue in config$tissue ){
  #pheatmap(Final, main=paste0(Tissue))
  writeLines(paste0(Missing[[ Tissue ]]))
  
  dev.off()
  plot.new()
  pdf( file = paste0(plotdir,'/', Tissue,'_Coexpression.pdf') )
    #pheatmap(Final, main=Tissue)
    eval( parse( text = paste0( 'pheatmap(Plots$',
                              Tissue,
                              '$ParCor, main=Plots$',
                              Tissue,
                              '$main, show_rownames=F)'  
                            )))
  Sys.sleep(2)
  dev.off()
}
dev.off()
#Push Plots to synapse
for( Tissue in config$tissue ){
  
  # Set annotations for synapse objects ( figures and tables )
  all.annotations = list(
    dataType = 'Coexpression',
    dataSubType = 'pheatmap PDF',
    summaryLevel = 'gene',
    assay	 = 'RNAseq',
    tissueTypeAbrv	= Tissue, 
    study = Study[Tissue], 
    organism = 'HomoSapiens',
    consortium	= 'AMPAD',
    normalizationStatus	= TRUE,
    normalizationType	= 'CQN',
    rnaquantification = 'RSEM',
    genomeAssemblyID = 'GRCh38'
  )
  
  syns_Used <- eval( parse( text= paste0( 'Syn$', Tissue ) ))
  #Tissue <- 'STG'
  activityName = paste0( Tissue, 'Coexpression heatmap')
  activityDescription = paste0( Tissue, 'Coexpression heatmap for genelist: ', config$genelistfile);
  
  #Push To Synapse:
  ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0(plotdir,'/', Tissue,'_Coexpression.pdf'), 
                               name = paste0(Tissue,' Coexpression.pdf'), 
                               parentId=CODE$properties$id ), 
                               used = syns_Used,
                               activityName = activityName, 
                               executed = list(GL_File, CF_File, I_File, M_File),
                               activityDescription = activityDescription)

  all.annotations$dataSubType = 'pheatmap PDF'
  syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)
}

#Save the Partial Cor table/List Object
write.table( FinalTis, file = paste0('runs/',config$runID,'/PartialCorrelationPerms.tsv'), quote = F, row.names = F, col.names = T,  )

# Set annotations for synapse objects ( figures and tables )
all.annotations = list(
  dataType = 'Coexpression',
  dataSubType = 'partial correlations',
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

syns_Used <- eval(parse(text=paste0( 'c(', paste0('Syn$', config$tissue, collapse =", " ), ')')))

activityName = paste0( Tissue, 'Coexpression heatmap')
activityDescription = paste0( Tissue, 'Coexpression heatmap for genelist: ', config$genelistfile)

CODE <- syn_temp$store(synapseclient$Folder(name = 'data', parentId = RunParent))

ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0('runs/',config$runID,'/PartialCorrelationPerms.tsv'), 
                                                   name = paste0('PartialCorrelationPerms.tsv'), 
                                                   parentId=CODE$properties$id ), 
                               used = syns_Used,
                               activityName = activityName, 
                               executed = list(GL_File, CF_File, I_File, M_File),
                               activityDescription = activityDescription)

all.annotations$dataSubType = 'Partial Correlations'
syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)

