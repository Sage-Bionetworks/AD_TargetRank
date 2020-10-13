
reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

library(Matrix)
library( DescTools)
library(factoextra)

GeneTrain <-   read.csv(file = syn_temp$get( 'syn23268088' )$path )
counts <- readMM( file = syn_temp$get('syn22975647')$path )
mathys_colnames <- read.table( file=syn_temp$get('syn22975687')$path )
mathys_rownames <- read.table( file=syn_temp$get('syn22965976')$path )
meta <- read.csv( file=syn_temp$get('syn22975699')$path, row.names = 1)

rownames(counts) <- mathys_rownames[,1]
colnames(counts) <- mathys_colnames[,1]

uncounts <- 2^counts
uncounts <- uncounts-1

Log_Counts <- log2(uncounts)


#Wizorize scale and realign the gene expression
Log_Counts <- as.data.frame(as.matrix(Log_Counts))
Log_Counts[ Log_Counts == -Inf] <- NA
Wizor_Log_Cnts <- apply( Log_Counts, 1, Winsorize, na.rm = T )
Scale_Win_Log_Cnts <- scale(Wizor_Log_Cnts)
Scale_Win_Log_Cnts <- as.data.frame( t(Scale_Win_Log_Cnts) )

write.csv( Scale_Win_Log_Cnts, 'Scaled_Winzorized_Log2Counts.csv' )
# Scale_Win_Log_Cnts <- as.data.frame( data.table::fread( 'Scaled_Winzorized_Log2Counts.csv', sep=',', check.names=FALSE ) )
# row.names( Scale_Win_Log_Cnts ) <- Scale_Win_Log_Cnts$V1
# Scale_Win_Log_Cnts <- Scale_Win_Log_Cnts[ , (colnames(Scale_Win_Log_Cnts) %in% 'V1') == F]

#Count Zeros By Wide Cell Type and Case Control: -- Start Here on Monday
Diagnosis_CT <- table( meta$broad.cell.type, meta$simpleDiagnosis )
Iter <- apply(expand.grid( colnames(Diagnosis_CT), row.names(Diagnosis_CT) ), 1, paste, collapse=".")

meta$broad.cell.type <- as.character( meta$broad.cell.type )
meta$Diagnosis <- as.character( meta$Diagnosis )

Cell_Div <- c("Cont","AD")
names(Cell_Div) <- c("Control","AD")

Summer <- function( inpt ){
  #return( sum(is.na(inpt)) )
  return(length( inpt[is.na(inpt)==F]  )/length(inpt))
}

MED <- function( inpt ){
  #return( sum(is.na(inpt)) )
  return(median( inpt[is.na(inpt)==F] ))
}

MEAN <- function( inpt ){
  #return( sum(is.na(inpt)) )
  return(mean( inpt[is.na(inpt)==F]  ))
}

Upper <- function( inpt ){
  #return( sum(is.na(inpt)) )
  return( as.numeric(quantile( inpt[is.na(inpt)==F], 0.975) ))
}
Lower <- function( inpt ){
  #return( sum(is.na(inpt)) )
  return( as.numeric(quantile( inpt[is.na(inpt)==F], 0.025) ))
}

Zeros <- list()
MEan <- list()
MEDIan <- list()
CI.L <- list()
CI.H <- list()
  
# 5% of cells
for( i in 1:length(Iter) ){
  
  Val <- strsplit( Iter[i] , '[.]')[[1]]
  Cells <- as.character( meta[ meta$simpleDiagnosis %in% names(Cell_Div[ Val[1] ]) & meta$broad.cell.type %in% Val[2], ]$TAG )
  CT_Cells <- as.character( meta[ meta$simpleDiagnosis %in% names(Cell_Div[ Val[1] ])  & meta$broad.cell.type %in% Val[2], ]$Subcluster )
  
  temp <- Scale_Win_Log_Cnts[ , paste0( CT_Cells, '_', as.character(Cell_Div[ Val[1] ]), '_', Cells) ]
  
  eval( parse( text=paste0( 'Zeros$', Iter[i], ' <- apply( temp , 1, Summer )' ) ) )
  eval( parse( text=paste0( 'MEan$', Iter[i], ' <- apply( temp , 1, MEAN )' ) ) )
  eval( parse( text=paste0( 'MEDIan$', Iter[i], ' <- apply( temp , 1, MED )' ) ) )
  eval( parse( text=paste0( 'CI.H$', Iter[i], ' <- apply( temp , 1, Upper )' ) ) )
  eval( parse( text=paste0( 'CI.L$', Iter[i], ' <- apply( temp , 1, Lower )' ) ) )
  
}

Template <- as.data.frame(matrix(NA, length(MEan$AD.Ast),8))
colnames( Template ) <-  c('GName', 
                           'ESNG',
                           'CellType',
                           'Diagnosis',
                           'Mean',
                           'Median',
                           'CI.L',
                           'CI.H'
                           )
Template$GName <- names( Zeros$AD.Ast )


#Remove Dpubles of Pseduo Genes:
#GeneTrain[ GeneTrain$gene_short_name %in% c( 'CCDC7', 'CRHR1', 'PAK6', 'PIK3R3', 'TBC1D26', 'TIMM10B', 'TMEM256-PLSCR3' ),]
GeneTrain <- GeneTrain[ (GeneTrain$ENSG_name %in% c( 'ENSG00000278139', 'ENSG00000265264', 
                                              'ENSG00000150076', 'ENSG00000259288',
                                              'ENSG00000187838', 'ENSG00000214946',
                                              'ENSG00000263715') == F
            ),]
row.names(GeneTrain) <- GeneTrain$gene_short_name
Template$ESNG <- GeneTrain[ Template$GName, ]$ENSG_name

write.csv( Template,'celltype_Gene_Expression_template.csv')


row.names( Template ) <- Template$GName
Filled_Templates <- list()

for( Cond in c( 'AD','Control') ){
  for( CT in row.names(Diagnosis_CT) ){
    Filled <- Template
    Filled$CellType <- CT
    Filled$Diagnosis <- Cond
    
    for( gene in Filled$GName ){
      if( eval(parse( text=paste0( 'Zeros$', Cond, '.', CT, '[\'', gene, '\'] > 0.05 ' ))) ){
        eval(parse( text=paste0( 'Filled[\'', gene, '\',]$Mean <- as.numeric( MEan$', Cond, '.', CT, '[\'', gene, '\'])' )))
        eval(parse( text=paste0( 'Filled[\'', gene, '\',]$Median <- as.numeric( MEDIan$', Cond, '.', CT, '[\'', gene, '\'])' )))
        eval(parse( text=paste0( 'Filled[\'', gene, '\',]$CI.L <- as.numeric( CI.L$', Cond, '.', CT, '[\'', gene, '\'])' )))
        eval(parse( text=paste0( 'Filled[\'', gene, '\',]$CI.H <- as.numeric( CI.H$', Cond, '.', CT, '[\'', gene, '\'])' )))
        
      }else{
        
      }
    }
    #Add to list object
    eval(parse( text=paste0( 'Filled_Templates$', Cond, '.', CT, ' <- Filled' ) ))
  }
}

Final_Prelim_sc <- do.call( rbind, Filled_Templates )
write.csv( Final_Prelim_sc, file= 'preliminary_scRNAExpression.csv',  row.names = FALSE )

parentId <- 'syn21534582'
activityName = 'Template For Cell Type Expression';
activityDescription = 'Template For Cell Type Expression';
CODE <- syn_temp$store(synapseclient$Folder(name = "Cell Type Specificity", parentId = parentId))

thisFileName <- 'Mathys_SingleCell_Processer.R'

# Github link
thisRepo <- githubr::getRepo(repository = "jgockley62/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))

#Set Used SynIDs For Provenance
Syns_Used <- c( 'syn22686575', 'syn22414716' )
# Set annotations
all.annotations = list(
  dataType = 'Network',
  summaryLevel = 'gene',
  assay	 = 'RNAseq',
  tissueTypeAbrv	= c('IFG', 'STG', 'FP', 'PHG', 'TCX', 'DLFPC'), 
  study = c( 'MSBB', 'ROSMAP', 'Mayo' ), 
  organism = 'HomoSapiens',
  consortium	= 'TreatAD',
  genomeAssemblyID = 'GRCh38'
)


ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='/home/jgockley/celltype_Gene_Expression_template.csv', name = 'GeneTemplate', parentId=CODE$properties$id ), executed = thisFile, used = Syns_Used, activityName = activityName, activityDescription = activityDescription)

ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='/home/jgockley/celltype_Gene_Expression_template.csv', name = 'GeneTemplate', parentId=CODE$properties$id ), executed = thisFile, used = Syns_Used, activityName = activityName, activityDescription = activityDescription)
