##Coexpression Analysis module:
#install.packages("~/Desktop/Programs/spike/spike/", repos = NULL, type = "source")

#config$tissue <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

# Github link
thisRepo <- githubr::getRepo(repository = "jgockley62/AD_TargetRank", ref="branch", refName='master')

setwd( config$filedir )

GL_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath=config$genelistfile)
CF_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath=config$name)
I_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/01-Initializer.r')
M_File  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/02-Master.Rmd')

CODE <- syn_temp$store(synapseclient$Folder(name = 'figures', parentId = RunParent))

#SynIDs of expression data to pull from
ExpressionDS <- c('syn21291908','syn21292041','syn21285564','syn21285564','syn21285564','syn21285564','syn21291908')
names( ExpressionDS ) <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

Study <- c( 'RosMap', 'Mayo', 'Mayo', 'MSBB', 'MSBB', 'MSBB', 'MSBB')
names(Study) <- c('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

row.names(Trans) <- Trans$ensembl_gene_id
Trans <- Trans[ Trans$chromosome_name %in% c(1:22,"X","Y"), ]

writeLines(paste0('For each tissue specified in ', 
                  paste( config$name, collapse = ', '), 
                  " the coexpression of genes specified in ", 
                  config$genelistfile, 
                  " is given. This is measured as a percentage of 1000 pemutations in which the gene was a significant predictor of expression of its comparison gene by a linear model."
                  )
           )
#List of Genes Not Expressed In Tissue
#dev.off()
Missing <- list()
Plots <- list()

for( Tissue in config$tissue ){
  
  # Set annotations
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
  
  #Load expression
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
  
  foo <- bcv::impute.svd( t(exp) )
  
  Exp <- foo$x
  row.names(Exp) <- row.names(t(exp)) 
  colnames(Exp) <- colnames(t(exp)) 
  
  x <- Exp[ , colnames(Exp) %in% row.names(Trans)]
  
  #Record Genes missing in tissue of interest
  if( as.numeric(table(row.names(Trans) %in% colnames(x) )["TRUE"] ) == dim(Trans)[1] ){
    Missing[[ Tissue ]] <- paste0("All Listed Genes Expressed in ", Tissue)
  }else{
    Missing[[ Tissue ]] <- paste0( "Queried Genes missing from ", Tissue, " : ",
                                  paste(c(Trans[ (row.names(Trans) %in% colnames(x) ) == F, ]$hgnc_symbol), collapse = ', ')
                                 )
  }
  
  if(length(Missing))
  Final <- matrix( NA, length(colnames(x)), length(colnames(x)) )
  row.names(Final) <- colnames(x)
  colnames(Final) <- colnames(x)
  
  for( i in colnames(x) ){
    OBS <- i
    y <- as.matrix(x[,OBS]) 
    colnames(y) <- i
    X <- x[,(colnames(x) %in% OBS) == F ]
    att <- spike::vbsrBootstrap( y=y,x=X,nsamp=100,cores=4)
    
    eval(parse(text = paste0("Final[ i,] <- c( att[ (names(att) %in% 'intercept')==F ], ", i, " = 1)[colnames(Final)]") ))
  }
  row.names(Final) <- Trans[ row.names(Final), ]$hgnc_symbol
  colnames(Final) <- Trans[ colnames(Final), ]$hgnc_symbol
  
  #Plots[[ Tissue ]] <- pheatmap::pheatmap(Final, main=Tissue)
  Plots[[ Tissue ]] <- list(Final, main=Tissue)
  #for( Tissue in config$tissue ){
  do.call( "pheatmap", Plots[[ Tissue ]])
  writeLines(Missing[[ Tissue ]])
  #}
  
  
  #plot.new()
  pdf( file = paste0(plotdir,'/', Tissue,'_Coexpression.pdf') )
    pheatmap::pheatmap(Final, main=Tissue)
  dev.off()
  #dev.off()
  
  #Push To Synapse:
  # Store SMR
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

