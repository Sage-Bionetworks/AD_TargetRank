# Code to Combine RNA and Proteomics Weights and Apply Scoring Harness
library(biomaRt)

reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

Prot <- read.table( syn_temp$get('syn22686575')$path, header=T, sep='\t', stringsAsFactors=F)
RNA  <- read.table( syn_temp$get('syn22414716')$path, header=T, sep='\t', stringsAsFactors=F)

ensgs <- c( Prot$ENSG, RNA$ENSG) 
ensgs <- ensgs[ !duplicated(ensgs) ]

#Build Combined File Template
Comb <- as.data.frame( matrix(NA, length(ensgs), 12) )
rownames(Comb) <- ensgs
colnames( Comb ) <- c( 'ENSG', 'GName',
                        'RNA_TE', 'RNA_fdr_CorPVal', 'RNA_Sig', 'RNA_Direction', 'RNA_Weight',
                        'Pro_TE', 'Pro_fdr_CorPVal', 'Pro_Sig', 'Pro_Direction', 'Pro_Weight'
                        )

#Seed the Dataframe
Comb$ENSG <- row.names( Comb )

prot <- Prot[ Prot$Type == 'Predicted Weight' , ]
row.names(prot) <- prot$ENSG

rna <- RNA[ RNA$Type == 'Predicted Weight' , ]
row.names(rna) <- rna$ENSG

#ProtVals
Comb[ row.names(prot), ]$Pro_Direction <- prot$Direction
Comb[ row.names(prot), ]$Pro_fdr_CorPVal <- prot$fdr.random
Comb[ row.names(prot), ]$Pro_Sig <- prot$Sig
Comb[ row.names(prot), ]$Pro_TE <- prot$TE.random
Comb[ row.names(prot), ]$GName <- prot$GName
Comb[ row.names(prot), ]$Pro_Weight <- prot$y

#RNA Vals
Comb[ row.names(rna), ]$RNA_fdr_CorPVal <- rna$fdr.random
Comb[ row.names(rna), ]$RNA_TE <- rna$TE.random
Comb[ row.names(rna), ]$GName <- rna$Gene
Comb[ row.names(rna), ]$RNA_Weight <- rna$PreAdjWeight

for( i in 1:dim(Comb)[1] ){
  if( is_NA(Comb[i,]$RNA_TE) ){
  }else{
    if( Comb[i,]$RNA_fdr_CorPVal < 0.05 ){
      Comb[i,]$RNA_Sig <- 'YES'
      if( Comb[i,]$RNA_TE < 0 ){
        Comb[i,]$RNA_Direction <- 'DOWN'
      }else{
        Comb[i,]$RNA_Direction <- 'UP'
      }
    }else{
      Comb[i,]$RNA_Sig <- 'NO'
      Comb[i,]$RNA_Direction <- 'NONE'
    }
  }
}


#Coherenane...
table( Comb$Pro_Sig, Comb$RNA_Sig )

temp <- Comb[ Comb$RNA_Sig == 'YES' &  Comb$Pro_Sig == 'YES', ]
temp<-temp[complete.cases(temp),]
table( temp$RNA_Direction, temp$Pro_Direction )

#Apply Harness:
Comb <- Comb[ order(-Comb$RNA_TE ),]
Comb$Harness <- NA
poo <- NULL
for( i in 1:dim(Comb)[1] ){
  
  #NA for either Gene
  if( is_NA(Comb[i,]$RNA_TE) || is_NA(Comb[i,]$Pro_TE) ){
    if( is_NA(Comb[i,]$RNA_TE)){
      
      if( Comb[i,]$Pro_fdr_CorPVal > 0.05 ){
        Comb[i,]$Harness <- ( 1.45 * Comb[i,]$Pro_Weight )
        #Comb[i,]$Harness <- .6*( 2 - abs(Comb[i,]$Pro_Weight ))
      }else{
        Comb[i,]$Harness <- ( 1.175 * Comb[i,]$Pro_Weight )
        #Comb[i,]$Harness <- ( 1.0 - abs(Comb[i,]$Pro_Weight ))
        
      }
      
    }else{
      if( Comb[i,]$RNA_fdr_CorPVal > 0.05 ){
        Comb[i,]$Harness <- ( 1.45 * Comb[i,]$RNA_Weight )
        #Comb[i,]$Harness <- .6 * ( 2 - abs(Comb[i,]$RNA_Weight ))
        
      }else{
        Comb[i,]$Harness <- ( 0.775 * Comb[i,]$RNA_Weight )
      }
    }
  }else{
    
    if( Comb[i,]$RNA_fdr_CorPVal > 0.05 && Comb[i,]$Pro_fdr_CorPVal > 0.05 ){
      Comb[i,]$Harness <- 1.45 * ( mean( Comb[i,]$RNA_Weight + Comb[i,]$Pro_Weight ) )
      #Comb[i,]$Harness <- .6* (2 - ( mean( abs(Comb[i,]$RNA_Weight) + abs(Comb[i,]$Pro_Weight )) ))
      poo <- c( poo, mean( abs(Comb[i,]$RNA_Weight) + abs(Comb[i,]$Pro_Weight )) )
    }else{
      if( Comb[i,]$RNA_fdr_CorPVal < 0.05 && Comb[i,]$Pro_fdr_CorPVal > 0.05 ){
        Comb[i,]$Harness <- .75 * ( Comb[i,]$RNA_Weight ) #+ .1 * ( Comb[i,]$Pro_Weight )
      }else{
        if( Comb[i,]$RNA_fdr_CorPVal > 0.05 && Comb[i,]$Pro_fdr_CorPVal < 0.05 ){
          Comb[i,]$Harness <- 1.15 * ( Comb[i,]$Pro_Weight ) #+ .1 * ( Comb[i,]$RNA_Weight )
        }else{
          Comb[i,]$Harness <- .95 * mean( .75 * Comb[i,]$RNA_Weight + 1.15 * Comb[i,]$Pro_Weight )
        }
      }
    }
  }
}

d <- density( Comb$Harness ) 
plot(d, las =1, bty='n')


#Harness....
Comb <- Comb[ order( -Comb$Harness ), ]
Comb$Rank <-  dim(Comb)[1]:1 / dim(Comb)[1]
Comb$Rank_Adj <- Comb$Rank  

mylogit <- glm(Comb$Rank_Adj  ~  poly( log(Comb$Harness),2)   , data = Comb, family = "binomial")
preds = predict(mylogit, newdata = list(age = Comb$Harnes), se = TRUE, type="response")


Comb$test <- preds$fit

Comb$TYPE <- NA
Comb[ Comb$RNA_Sig %in% 'YES', ]$TYPE <- 'Only_RNA'
Comb[ Comb$Pro_Sig %in% 'YES', ]$TYPE <- 'Only_Pro'
Comb[ Comb$Pro_Sig %in% 'YES' & Comb$RNA_Sig %in% 'YES', ]$TYPE <- 'Both_Sig'
Comb[ (Comb$Pro_Sig %in% 'NO' | is.na(Comb$Pro_Sig)) & (Comb$RNA_Sig %in% 'NO' | is.na(Comb$RNA_Sig)), ]$TYPE <- 'None_Sig'
Comb$TYPE<-as.factor(Comb$TYPE)

ggplot( data=Comb ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.2, width=0.05) 

ggplot( data=Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ] ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05) 


Comb$TYPE<-as.character(Comb$TYPE)

Comb$TYPE <- factor(Comb$TYPE, levels = c( "None_Sig", "Only_RNA", 'Only_Pro', "Both_Sig"),ordered = TRUE)

ggplot() + 
  geom_violin( data=Comb, aes(y = test, x=TYPE, fill=TYPE) )
  
Comb <- Comb[order(-Comb$test),]

head(Comb[, c('GName','test' )],n=20)
head(Comb[, c('ENSG', 'GName', 'RNA_TE', 'RNA_fdr_CorPVal', 'Pro_TE', 'Pro_fdr_CorPVal','TYPE', 'test' )],n=50)


Comb$Final_Weight <- Comb$test
Comb[ (Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES') == F , ]$Final_Weight <- 0
Comb$Final_Weight<- Comb$Final_Weight + 1-max(Comb$Final_Weight)
Comb[ (Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES') == F , ]$Final_Weight <- 0

summary(Comb$Final_Weight)

summary(Comb[ Comb$TYPE == 'Both_Sig', ]$Final_Weight)
summary(Comb[ Comb$TYPE == 'Only_Pro', ]$Final_Weight)
summary(Comb[ Comb$TYPE == 'Only_RNA', ]$Final_Weight)

Comb <- Comb[ order(-Comb$Final_Weight), ]
#head( Comb )

#### Write to Disk:

pdf('~/AD_TargetRank/Target_Type_Filt.pdf')
ggplot( data=Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ] ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05) + 
  ylab("Weight")
dev.off()

pdf('~/AD_TargetRank/Target_Type_Filt_Box.pdf')
ggplot() + 
  geom_violin( data=Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ], aes(y = test, x=TYPE, fill=TYPE) ) + 
  ylab("Weight")
dev.off()

pdf('~/AD_TargetRank/Target_Type_Tot.pdf')
ggplot( data=Comb ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05)  + 
  ylab("Weight")
dev.off()

pdf('~/AD_TargetRank/Target_Type_Tot_Box.pdf')
ggplot() + 
  geom_violin( data=Comb, aes(y = test, x=TYPE, fill=TYPE) ) + 
  ylab("Weight")
dev.off()

Comb <- Comb[ , c('ENSG', 'GName', 'RNA_TE', 'RNA_fdr_CorPVal', 'RNA_Sig', 'RNA_Direction', 'RNA_Weight', 
                  'Pro_TE', 'Pro_fdr_CorPVal', 'Pro_Sig', 'Pro_Direction', 'Pro_Weight', 'Harness', 'Rank',
                  'test', 'TYPE', 'Final_Weight')]
colnames(Comb)[ colnames(Comb) == 'test'] <- 'Predicted_Weight'
Comb$OmicsScore <- 2 * Comb$Final_Weight

write.csv( Comb, '~/AD_TargetRank/OmicsScores.csv', 
           row.names = F, quote = F)

parentId = 'syn22414618';
activityName = 'Second Generation of Omics Scoring Harness';
activityDescription = 'Omics Scoring harnesss based on RNA-Seq and Proteomics meta analysis';
thisFileName <- 'Omics_Harness.R'

# Github link
library(githubr)
thisRepo <- getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))

CODE <- syn_temp$store(synapseclient$Folder(name = "Plots", parentId = parentId))
Syns_Used <- c( 'syn22686575', 'syn22414716' )


ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='~/AD_TargetRank/Target_Type_Filt.pdf', 
                                                   name = 'Distribution of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='~/AD_TargetRank/Target_Type_Filt_Box.pdf', 
                                                   name = 'Box Plots of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='~/AD_TargetRank/Target_Type_Tot.pdf', 
                                                   name = 'Distribution of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='~/AD_TargetRank/Target_Type_Tot_Box.pdf', 
                                                   name = 'Box Plots of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

CODE <- syn_temp$store(synapseclient$Folder(name = "MetaAnalysis", parentId = 'syn22351719'))
Syns_Used <- c( 'syn22686575', 'syn22414716' )

ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='~/AD_TargetRank/OmicsScores.csv', 
                                                   name = 'Omics Scores Version 2.0', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

