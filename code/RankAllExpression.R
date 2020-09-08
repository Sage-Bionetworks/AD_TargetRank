## This script Creates the RNA-Expression Raw Ranks for every gene Across all Data Sets
## Non-significant DE between Case and Control are pushed to Zero

library(yaml)

library(readr)
library(biomaRt)
library(knitr)
library(kableExtra)
library(tidyverse)
library(grid)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(yaml)

setwd('~/AD_TargetRank')

reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

print(getwd())
#Assign your YAML formated config file to a varable
#configuration <- "configs/EmoryTargets.yaml"
#config <- read_yaml(file=paste0('./',configuration))

setwd(config$filedir)
plotdir=paste0( 'runs/totalRank/figures' )
tabledir=paste0( 'runs/totalRank/tables' )

Tissues <- config$tissue

#Create Synapse Run Folder and retrieve parentID
CODE <- syn_temp$store(synapseclient$Folder(name = 'totalRank', parentId = 'syn21532475'))

foo <- syn_temp$getChildren('syn21532475')
Temp <- reticulate::iterate(foo, simplify = F)
for( i in 1:length(Temp) ){
  if( Temp[[i]]$name == 'totalRank' ){
    RunParent <- Temp[[i]]$id
  }else{
  }
}


#Load in the differential expression data
DE <- read.table( syn_temp$get('syn21534583')$path, header=T, sep ="\t")
#Genes <- Tab$ENSG
Genes <- as.character(DE$ensembl_gene_id)
Genes <- Genes[!duplicated(Genes)]

#Construct the rank model
ModelConstruct <- function( de, Comp, Mod, Sex ){
  #C#reates a model and plots the target genes
  #'@de the differential  expression Construct (DE)
  #'@Comp The comparison ie NA or "AD-CONTROL" 
  #'@Mod The model is NA or 'Diagnosis'
  #'@Sex The Sex Comparison ie NA or 'MALE' or 'FEMALE'
  
  #de<- DE
  #Comp <- 'AD-CONTROL'
  #Mod <- 'Diagnosis'
  #Sex <- 'ALL' 
  
  #Find the gene expression for the value of the genes where it exists
  if( isTRUE(is.na(Comp)) ){
  }else{
    de <- de[ de$Comparison == Comp ,  ] 
  }
  if( isTRUE(is.na(Mod)) ){
  }else{
    de <- de[ de$Model == Mod ,   ] 
  }
  if( isTRUE(is.na(Sex)) ){
  }else{
    de <- de[ de$Sex == Sex ,  ] 
  }
  de <- de[ , c('hgnc_symbol', 'ensembl_gene_id','Tissue','logFC','adj.P.Val') ]
  
  Models <- list()
  Plots <- list()
  
  #Loop through induvidual tissues
  for( Tissue in c('CBE', 'DLPFC', 'FP', "IFG", "PHG", "STG", "TCX" )){
    Work <- de[ de$Tissue == Tissue, ]
    Work$y <- rank(abs( Work$logFC)) / max(rank(abs( Work$logFC)))
    Work$Type <- "Actual Rank"
    #plot( log(abs( Work$logFC)), Work$y )
    
    #Fit Logistic Model:
    #mylogit <- glm(Work$y  ~ abs( Work$logFC), data = FOO, family = "binomial")
    mylogit <- glm(Work$y  ~ abs( Work$logFC), data = Work, family = "binomial")
    
    Work$log_abs_PVal <- abs( Work$logFC)
    
    Work2 <- Work  
    Work2$Type <- "Predicted Weight"
    #Work2$y <- y <- 1/(1+exp( -( mylogit$coefficients[1]+mylogit$coefficients[2]* abs( Work$logFC) ) ))
    Work2$y <- 1/(1+exp( -( mylogit$coefficients[1]+mylogit$coefficients[2]* abs( Work$logFC) ) ))
    
    tWork <- as.data.frame(rbind( Work ,Work2 ))
    tWork$Sig <- ifelse(tWork$adj.P.Val < 0.05, "YES", "NO" )
    eval(parse( text=paste0( 'Models$', Tissue, ' <- as.data.frame( tWork )' ) ))
    
    #Colors for the sig and non-sig names
    COLS <- c( "blueviolet","grey29")
    names(COLS) <- c("YES","NO")
    
    #Plot the model
  
  }
  return( list( Plots=Plots, Models=Models ) )
}

#Plot/Run/Build the DE Models
ALL_Mod <- ModelConstruct( DE, 'AD-CONTROL', 'Diagnosis', 'ALL' )

#Pull All DE Data and push non-significant values to zero:

RNA <- list()
for( i in names(ALL_Mod$Models) ){
  print(i)
  eval(parse(text = paste0( 'RNA$', i, ' <- ALL_Mod$Models$',i,'[ALL_Mod$Models$',i,'$Type==\'Predicted Weight\', ]' )))
  eval(parse(text = paste0( 'RNA$', i, '<- RNA$', i, '[ is.na(RNA$', i, '$Sig)==F,]' )))
  eval(parse(text = paste0( 'RNA$',i,'[ RNA$',i,'$Sig == \'NO\',]$y <- 0' )))
  eval(parse(text = paste0( 'RNA$',i,' <- RNA$',i,'[  !duplicated(RNA$',i,'$ensembl_gene_id), ] ' )))
  eval(parse(text = paste0( 'row.names( RNA$',i,' ) <-RNA$',i,'$ensembl_gene_id ' )))
}

Transctipt_Weight <- as.data.frame( matrix( NA, length(Genes), 9 ) )
colnames(Transctipt_Weight) <- c( 'ENSG', 'GeneName', "DLPFC", "FP", "IFG", "PHG", "STG", "TCX", "CBE")
row.names(Transctipt_Weight) <- Genes
Transctipt_Weight$ENSG <- Genes

for( i in names(ALL_Mod$Models) ){
  print(i)
  eval(parse(text = paste0( 'Transctipt_Weight[ row.names(RNA$', i,') ,]$', i, '<- RNA$', i, '$y' )))
  eval(parse(text = paste0( 'RNA$', i, '$hgnc_symbol <- as.character(RNA$', i, '$hgnc_symbol)' )))
  eval(parse(text = paste0( 'Transctipt_Weight[ row.names(RNA$', i,') ,]$GeneName <- RNA$', i, '$hgnc_symbol' )))
}

Transctipt_Weight$Cortex6 <- NA
Transctipt_Weight$All7 <- NA

Meaner6 <- function(boo){
  mean( boo[,c('DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX' )][ is.na(boo[,c('DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX' )])==F ])
}
Meaner7 <- function(boo){
  mean( boo[,c('DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX','CBE' )][ is.na(boo[,c('DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX', 'CBE' )])==F ])
}


for(i in 1:dim(Transctipt_Weight)[1] ){
  Transctipt_Weight[ i, ]$Cortex6 <- Meaner6( Transctipt_Weight[ i, ] )
}
for(i in 1:dim(Transctipt_Weight)[1] ){
  Transctipt_Weight[ i, ]$All7 <- Meaner7( Transctipt_Weight[ i, ] )
}

write.csv(Transctipt_Weight,'runs/totalRank/RNA_Seq_GeneWeights.csv')

####Proteomics Score:
Genes <- as.character(Prot$ENSG)[ !duplicated(as.character(Prot$ENSG)) ]
  
Prot <- read.csv( file = syn_temp$get('syn21534585')$path, header = T)
Prot <- Prot[ complete.cases(Prot[,c("Log2_FC", "CI_Upr", "CI_Lwr", "PVal", "Cor_PVal")]), ]
Prot <- Prot[complete.cases(Prot),]

#TargetP <- Prot[ Prot$ENSG %in% Genes,]
TargetP <- Prot
  
Names <- paste0(TargetP$UniqID, '|', TargetP$ENSG)
Names <- Names[!duplicated(Names)]
PTab <- do.call( rbind, strsplit( Names, '\\|' ))

if( table(Genes %in% PTab[,3])['TRUE'] == length(Genes) ){
  Res <- PTab
}else{
  Adds <- Genes[ (Genes %in% PTab[,3]) == F ]
  Res <- rbind( PTab, cbind( Tab[Adds,]$GeneName, NA, Adds ))
}

row.names(Res) <- paste0( Res[,1], '|', Res[,2])
Res <- as.data.frame(Res)
colnames( Res ) <- c( "GeneName", "UniProt", "ENSG")

Tiss <- as.character(TargetP$Tissue)[ !duplicated(as.character(TargetP$Tissue)) ]
Adds <- matrix(0, dim(Res)[1], 4)
row.names(Adds) <- row.names(Res)
colnames(Adds) <- c("AntPFC", "DLPFC", "MFG", "TCX" )

TotP <- as.data.frame(cbind(Res, Adds))
TotPval <- as.data.frame(cbind(Res, Adds))

row.names(TotP) <- gsub('\\.', '|', row.names(TotP))
row.names(TotPval) <- gsub('\\.', '|', row.names(TotPval))

for( row in row.names(TotP)){
  for( col in colnames(TotP)){
    if( isTRUE(row %in%TargetP$UniqID) ){
      TP <- TargetP[ TargetP$UniqID == row,]
      if( isTRUE(col %in% TP$Tissue) ){
        TotP[row,col] <- signif( TP[ as.character(TP$Tissue)==col ,]$Log2_FC, 3 )
      }else{}
    }else{
      
    }
  }
}

for( row in row.names(TotPval)){
  for( col in colnames(TotPval)){
    if( isTRUE(row %in%TargetP$UniqID) ){
      TP <- TargetP[ TargetP$UniqID == row,]
      if( isTRUE(col %in% TP$Tissue) ){
        TotPval[row,col] <- TP[ as.character(TP$Tissue)==col ,]$Cor_PVal
      }else{}
    }else{
      
    }
  }
}


#ProtPlots
ModelConstructProt <- function( de ){
  #C#reates a model and plots the target genes
  #'@de the differential  expression Construct (DE)
  #de <- DE
  #de <- Prot
  #Mod <- 'Diagnosis'
  #Sex <- 'ALL'
  
  
  #de <- de[ , c('hgnc_symbol', 'ensembl_gene_id','Tissue','logFC','adj.P.Val') ]
  de <- de[ , c("GeneName", "ENSG", "Tissue", 'Log2_FC', 'Cor_PVal') ]
  colnames(de) <- c('hgnc_symbol', 'ensembl_gene_id','Tissue','logFC','adj.P.Val')
  
  Modelz <- list()
  Plotz <- list()
  
  #Loop through induvidual tissues
  for( Tissue in c("AntPFC", "DLPFC", "MFG", "TCX" )){
    Work <- de[ de$Tissue == Tissue, ]
    Work$y <- rank(abs( Work$logFC)) / max(rank(abs( Work$logFC)))
    Work$Type <- "Actual Rank"
    #plot( log(abs( Work$logFC)), Work$y )
    
    #Fit Logistic Model:
    #mylogit <- glm(Work$y  ~ abs( Work$logFC), data = FOO, family = "binomial")
    mylogit <- glm(Work$y  ~ abs( Work$logFC), data = Work, family = "binomial")
    
    Work$log_abs_PVal <- abs( Work$logFC)
    
    Work2 <- Work  
    Work2$Type <- "Predicted Weight"
    #Work2$y <- y <- 1/(1+exp( -( mylogit$coefficients[1]+mylogit$coefficients[2]* abs( Work$logFC) ) ))
    Work2$y <- 1/(1+exp( -( mylogit$coefficients[1]+mylogit$coefficients[2]* abs( Work$logFC) ) ))
    
    tWork <- as.data.frame(rbind( Work ,Work2 ))
    tWork$Sig <- ifelse(tWork$adj.P.Val < 0.05, "YES", "NO" )
    eval(parse( text=paste0( 'Modelz$', Tissue, ' <- as.data.frame( tWork )' ) ))
    
    #Colors for the sig and non-sig names
    COLS <- c( "blueviolet","grey29")
    names(COLS) <- c("YES","NO")
    
  }
  return( list( Plotz=Plotz, Modelz=Modelz ) )
}

Pro_ALL_Mod <- ModelConstructProt( Prot )

PROT <- list()
for( i in names(Pro_ALL_Mod$Modelz) ){
  print(i)
  eval(parse(text = paste0( 'PROT$', i, ' <- Pro_ALL_Mod$Modelz$',i,'[Pro_ALL_Mod$Modelz$',i,'$Type==\'Predicted Weight\', ]' )))
  eval(parse(text = paste0( 'PROT$', i, '<- PROT$', i, '[ is.na(PROT$', i, '$Sig)==F,]' )))
  eval(parse(text = paste0( 'PROT$',i,'[ PROT$',i,'$Sig == \'NO\',]$y <- 0' )))
  #eval(parse(text = paste0( 'PROT$',i,' <- PROT$',i,'[  !duplicated(PROT$',i,'$ensembl_gene_id), ] ' )))
  #eval(parse(text = paste0( 'row.names( PROT$',i,' ) <-PROT$',i,'$ensembl_gene_id ' )))
}


#Take care of ENSGs with multiple perpresentations
#Average the predicted values of Significant genes
#Ignore zero vals as they could be a different isoform or bad peptide feature

#table( PROT$DLPFC[ PROT$DLPFC$Sig == 'YES', ]$ensembl_gene_id )[ table( PROT$DLPFC[ PROT$DLPFC$Sig == 'YES', ]$ensembl_gene_id ) > 1 ]

Prot <- list()
for( DAT in PROT ){
  Reps <- names( table( DAT$ensembl_gene_id )[ table( DAT$ensembl_gene_id ) > 1 ] )
  DAT$Keeps <- 'Y'
  DAT$ensembl_gene_id <- as.character(DAT$ensembl_gene_id)
  for( i in Reps){
    if( 'YES' %in% DAT[ DAT$ensembl_gene_id  == i ,]$Sig){
      if( 'NO' %in% DAT[ DAT$ensembl_gene_id  == i ,]$Sig){
        DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'NO',]$Keeps <- "N"
        DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$y <- mean( DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$y )
        len <- length( DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$Keeps )
        if(len > 1){
          DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$Keeps[2:len] <- 'N'
        }else{}
      }else{
        DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$y <- mean( DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$y )
        len <- length( DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$Keeps )
        DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'YES',]$Keeps[2:len] <- 'N'
      }
    }else{
      len <- length( DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'NO',]$Keeps )
      DAT[ DAT$ensembl_gene_id  == i & DAT$Sig == 'NO',]$Keeps[2:len] <- 'N'
    }
  }
  DAT$Keeps <- as.character(DAT$Keeps)
  eval(parse( text= paste0( 'Prot$', as.character(DAT$Tissue)[1], ' <- DAT[ DAT$Keeps == \'Y\', ]' )))
  eval(parse( text= paste0( 'row.names(Prot$', as.character(DAT$Tissue)[1], ') <- as.character(Prot$', as.character(DAT$Tissue)[1], '$ensembl_gene_id)') ))
}

PGnames <- c( row.names( Prot$DLPFC ),
              row.names( Prot$AntPFC ),
              row.names( Prot$MFG ),
              row.names( Prot$TCX ) 
)
PGnames <- PGnames[!duplicated(PGnames)]

Prot_Weight <- as.data.frame( matrix( NA, length(PGnames), 6 ) )
colnames(Prot_Weight) <- c( 'ENSG', 'GeneName', "DLPFC", "AntPFC", "MFG", "TCX" )
row.names(Prot_Weight) <- PGnames
Prot_Weight$ENSG <- PGnames

for( i in names(Prot) ){
  print(i)
  eval(parse(text = paste0( 'Prot_Weight[ row.names(Prot$', i,') ,]$', i, '<- Prot$', i, '$y' )))
  eval(parse(text = paste0( 'Prot$', i, '$hgnc_symbol <- as.character(Prot$', i, '$hgnc_symbol)' )))
  eval(parse(text = paste0( 'Prot_Weight[ row.names(Prot$', i,') ,]$GeneName <- Prot$', i, '$hgnc_symbol' )))
}

Prot_Weight$Avg4 <- NA

Meaner4 <- function(boo){
  mean( boo[,c( "AntPFC", "DLPFC", "MFG", "TCX" )][ is.na(boo[,c("AntPFC", "DLPFC", "MFG", "TCX" )])==F ])
}


for(i in 1:dim(Prot_Weight)[1] ){
  Prot_Weight[ i, ]$Avg4 <- Meaner4( Prot_Weight[ i, ] )
}

write.csv(Prot_Weight, 'runs/totalRank/Protien_Seq_GeneWeights.csv')

##############################################################################################################################
##############################################################################################################################
#Meta Analysis Weights.
