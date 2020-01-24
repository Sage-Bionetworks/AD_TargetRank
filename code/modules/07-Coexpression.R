##Coexpression Analysis module:
#install.packages("~/Desktop/Programs/spike/spike/", repos = NULL, type = "source")

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
  
  Syns_Used <- NULL
  #Tissue <- 'STG'
  activityName = paste0( Tissue, 'Coexpression heatmap')
  activityDescription = paste0( Tissue, 'Coexpression heatmap for genelist: ', config$genelistfile);
  
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
  
  #Prep for partial correlation detection
  Final <- data.frame()
  core <- parallel::detectCores()-2 
  source("~/AD_TargetRank/utilityFunctions/Parallel_vbsrBootstrap.R")
  #Run Partial correlations for each gene
  mark <- Sys.time()
  for( i in row.names(Trans) ){
    OBS <- i
    y <- as.matrix(x[,OBS]) 
    colnames(y) <- i
    X <- x[,(colnames(x) %in% OBS) == F ]
    
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
  
  #Replace ENSGs with gene names for plotting
  
  
  
  #row.names(Final) <- Trans[ row.names(Final), ]$hgnc_symbol
  #colnames(Final) <- Trans[ colnames(Final), ]$hgnc_symbol
  
  #Write Plot info to list, plot to file, and plot to wiki
  Plots[[ Tissue ]] <- list(ParCor = Final, main=Tissue)
}

for( Tissue in config$tissue ){
  eval( parse( text = paste0( 'pheatmap(Plots$',
                              Tissue,
                              '$ParCor, main=Plots$',
                              Tissue,
                              '$main)' 
                         )))
} 
for( Tissue in config$tissue ){
  #pheatmap(Final, main=paste0(Tissue))
  writeLines(paste0(Missing[[ Tissue ]]))
  
  #dev.off()
  #plot.new()
  pdf( file = paste0(plotdir,'/', Tissue,'_Coexpression.pdf') )
    #pheatmap(Final, main=Tissue)
    eval( parse( text = paste0( 'pheatmap(Plots$',
                              Tissue,
                              '$ParCor, main=Plots$',
                              Tissue,
                              '$main)'  
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
  
  #Syns_Used <- NULL
  #Tissue <- 'STG'
  activityName = paste0( Tissue, 'Coexpression heatmap')
  activityDescription = paste0( Tissue, 'Coexpression heatmap for genelist: ', config$genelistfile);
  
  #Push To Synapse:
  ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path=paste0(plotdir,'/', Tissue,'_Coexpression.pdf'), 
                               name = paste0(Tissue,' Coexpression.pdf'), 
                               parentId=CODE$properties$id ), 
                               used = Syns_Used,
                               activityName = activityName, 
                               executed = list(GL_File, CF_File, I_File, M_File),
                               activityDescription = activityDescription)

  all.annotations$dataSubType = 'pheatmap PDF'
  syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)
}
