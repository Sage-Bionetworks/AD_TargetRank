# Code to Combine RNA and Proteomics Weights and Apply Scoring Harness
library(biomaRt)

#reticulate::use_python("/usr/bin/python", required = TRUE)
#synapseclient <- reticulate::import("synapseclient")
#syn_temp <- synapseclient$Synapse()
#syn_temp$login()

synapser::synLogin()

Prot <- read.table( synapser::synGet('syn22686575')$path, header=T, sep='\t', stringsAsFactors=F)
RNA  <- read.table( synapser::synGet('syn22414716')$path, header=T, sep='\t', stringsAsFactors=F)

Prot[ Prot$peptide_id =='KIAA1107|Q9UPP5-2', ]$EnsemblGene <- 'ENSG00000189195'
Prot[ Prot$peptide_id =='KIAA1107|Q9UPP5-2', ]$SYMBOL <- 'BTBD8'

# Remove the one ID with no ENSG: 0|P04434
Prot <- Prot[grepl('ENSG',Prot$EnsemblGene),]

Prot$ENSG <- Prot$EnsemblGene
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

#Reduce Duplicated Gene Features to favor isoforms with greatest weight
prot$keep <- 'No'
ensgs <- prot$ENSG[!duplicated(prot$ENSG)]
for(gene in ensgs){
  if(as.numeric(table(gene == prot$ENSG)['TRUE']) == 1){
    prot[ prot$ENSG %in% gene, ]$keep <- 'Yes'
  }else{
    #If none are significant
    if((as.numeric(table('YES' == prot[ prot$ENSG %in% gene, ]$Sig)['FALSE']) == dim(prot[ prot$ENSG %in% gene, ])[1]) | is.na(as.numeric(table('YES' == prot[ prot$ENSG %in% gene, ]$Sig)['FALSE']))){
      weights <- prot[ prot$ENSG %in% gene, ]$y
      prot[ prot$ENSG %in% gene , ][1,]$keep <- 'Yes'
    }else{
      # If theres only one significant isoform
      if(as.numeric(table('YES' == prot[ prot$ENSG %in% gene, ]$Sig)['TRUE']) == 1){
        prot[ prot$ENSG %in% gene & prot$Sig == 'YES', ]$keep <- 'Yes'
      }else{
        # If multiple isoforms are significant
        if(as.numeric(table('YES' == prot[ prot$ENSG %in% gene, ]$Sig)['TRUE']) > 1){
          weights <- prot[ prot$ENSG %in% gene & prot$Sig == 'YES', ]$y_Final
          prot[ prot$ENSG %in% gene & prot$Sig == 'YES' & prot$y_Final == max(weights), ]$keep <- 'Yes'
        }
      }
    }
  }
}
prot <- prot[ prot$keep == 'Yes',]
row.names(prot) <- prot$ENSG

rna <- RNA[ RNA$Type == 'Predicted Weight' , ]
row.names(rna) <- rna$ENSG

#ProtVals
Comb[ row.names(prot), ]$Pro_Direction <- prot$Direction
Comb[ row.names(prot), ]$Pro_fdr_CorPVal <- prot$fdr.random
Comb[ row.names(prot), ]$Pro_Sig <- prot$Sig
Comb[ row.names(prot), ]$Pro_TE <- prot$TE.random
Comb[ row.names(prot), ]$GName <- prot$SYMBOL
Comb[ row.names(prot), ]$Pro_Weight <- prot$y

#RNA Vals
Comb[ row.names(rna), ]$RNA_fdr_CorPVal <- rna$fdr.random
Comb[ row.names(rna), ]$RNA_TE <- rna$TE.random
Comb[ row.names(rna), ]$GName <- rna$Gene
Comb[ row.names(rna), ]$RNA_Weight <- rna$PreAdjWeight

## is.na vs is_NA
for( i in 1:dim(Comb)[1] ){
  if( is.na(Comb[i,]$RNA_TE) ){
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
  if( is.na(Comb[i,]$RNA_TE) || is.na(Comb[i,]$Pro_TE) ){
    if( is.na(Comb[i,]$RNA_TE)){
      
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

p <- ggplot2::ggplot(Comb, ggplot2::aes(x=Harness)) + 
  ggplot2::geom_density()

#d <- density( Comb$Harness ) 
#plot(d, las =1, bty='n')


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

Comb$plotname <- NA
#label_name <- c('CPEB1', 'SIRT2', 'DLG5', 'MTCH1', 'MTOR', 'ADAM17', 'APPL1', 'CNTNAP2', 'NEFL')
label_name <- c("MTCH1", "MTOR", "APPL1", 'BAX','STMN2', 'SLC9A9','TOMM40L','GABRA2', 'MTCH2')
Comb[ Comb$GName %in% label_name, ]$plotname <- Comb[ Comb$GName %in% label_name, ]$GName

ggplot( data=Comb ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.2, width=0.05) 

ggplot( data=Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ] ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05) +
  ggrepel::geom_label_repel(data = Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ] , 
                            ggplot2::aes(x=Harness, y=test, label = plotname, size = NULL, color = NULL),
                            nudge_x = .45,
                            nudge_y = -.15,
                            segment.size  = 0.2,
                            segment.color = "grey50",
                            direction     = "x"
  )


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
plots <- '~/Desktop/Projects/TREAT_AD/Figures/'
pdf(paste0(plots,'Target_Type_Filt.pdf'))
ggplot( data=Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ] ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05) + 
  ylab("Weight") +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05) +
  ggrepel::geom_label_repel(data = Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ] , 
                            ggplot2::aes(x=Harness, y=test, label = plotname, size = NULL, color = NULL),
                            nudge_x = .45,
                            nudge_y = -.15,
                            segment.size  = 0.2,
                            segment.color = "grey50",
                            direction     = "x",
  )

dev.off()

pdf(paste0(plots,'Target_Type_Filt_Box.pdf'))
ggplot() + 
  geom_violin( data=Comb[ Comb$RNA_Sig %in% 'YES' | Comb$Pro_Sig %in% 'YES', ], aes(y = test, x=TYPE, fill=TYPE) ) + 
  ylab("Weight")
dev.off()

pdf(paste0(plots,'Target_Type_Tot.pdf'))
ggplot( data=Comb ) + 
  geom_line( aes( y = test, x=Harness  )) +
  geom_jitter( aes( y = test, x=Harness, col= TYPE),alpha = 0.4, shape = 16, size = .5, width=0.05)  + 
  ylab("Weight")
dev.off()

pdf(paste0(plots,'Target_Type_Tot_Box.pdf'))
ggplot() + 
  geom_violin( data=Comb, aes(y = test, x=TYPE, fill=TYPE) ) + 
  ylab("Weight")
dev.off()

Comb <- Comb[ , c('ENSG', 'GName', 'RNA_TE', 'RNA_fdr_CorPVal', 'RNA_Sig', 'RNA_Direction', 'RNA_Weight', 
                  'Pro_TE', 'Pro_fdr_CorPVal', 'Pro_Sig', 'Pro_Direction', 'Pro_Weight', 'Harness', 'Rank',
                  'test', 'TYPE', 'Final_Weight')]
colnames(Comb)[ colnames(Comb) == 'test'] <- 'Predicted_Weight'
Comb$OmicsScore <- 2 * Comb$Final_Weight

write.csv( Comb, '~/Desktop/Projects/TREAT_AD/OmicsScores.csv', 
           row.names = F, quote = F)

parentId = 'syn22414618';
activityName = 'Second Generation of Omics Scoring Harness';
activityDescription = 'Omics Scoring harnesss based on RNA-Seq and Proteomics meta analysis';
thisFileName <- 'Omics_Harness.R'

# Github link
library(githubr)
thisRepo <- getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))


CODE <- synapser::synStore(synapser::Folder(name = "Plots", parentId = parentId))
Syns_Used <- c( 'syn22686575', 'syn22414716' )

ENRICH_OBJ <- synapser::synStore( synapser::File( path=paste0(plots,'Target_Type_Filt.pdf'), 
                                                   name = 'Distribution of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

ENRICH_OBJ <- synapser::synStore( synapser::File( path=paste0(plots,'Target_Type_Filt_Box.pdf'), 
                                                   name = 'Box Plots of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

ENRICH_OBJ <- synapser::synStore( synapser::File( path=paste0(plots,'Target_Type_Tot.pdf'), 
                                                   name = 'Distribution of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

ENRICH_OBJ <- synapser::synStore( synapser::File( path=paste0(plots,'Target_Type_Tot_Box.pdf'), 
                                                   name = 'Box Plots of weights by Omics Modality', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)

CODE <- synapser::synStore(synapser::Folder(name = "MetaAnalysis", parentId = 'syn22351719'))
Syns_Used <- c( 'syn22686575', 'syn22414716' )

ENRICH_OBJ <- synapser::synStore( synapser::File( path='~/Desktop/Projects/TREAT_AD/OmicsScores.csv', 
                                                   name = 'Omics Scores Version 2.0', 
                                                   parentId=CODE$properties$id ), used = Syns_Used, 
                               activityName = activityName, 
                               executed = thisFile, 
                               activityDescription = activityDescription)
#Push Table
Comb$RNA_TE <- as.double(Comb$RNA_TE) 
Comb$RNA_fdr_CorPVal <- as.double(Comb$RNA_fdr_CorPVal) 
Comb$RNA_Weight <- as.double(Comb$RNA_Weight) 

Comb$Pro_TE <- as.double(Comb$Pro_TE) 
Comb$Pro_fdr_CorPVal <- as.double(Comb$Pro_fdr_CorPVal) 
Comb$Pro_Weight <- as.double(Comb$Pro_Weight) 
Comb$Harness <- as.double(Comb$Harness) 
Comb$Rank <- as.double(Comb$Rank) 

Comb$Predicted_Weight <- as.double(Comb$Predicted_Weight) 
Comb$Final_Weight <- as.double(Comb$Final_Weight) 
Comb$OmicsScore <- as.double(Comb$OmicsScore) 

Comb$Pro_fdr_CorPVal <- signif(Comb$Pro_fdr_CorPVal , digits = 6)
Comb$RNA_fdr_CorPVal <- signif(Comb$RNA_fdr_CorPVal , digits = 6)
Comb$Harness <- signif(Comb$Harness , digits = 6)
Comb$Predicted_Weight <- signif(Comb$Predicted_Weight , digits = 6)
Comb$RNA_TE <- signif(Comb$RNA_TE , digits = 6)
Comb$Pro_TE   <- signif(Comb$Pro_TE   , digits = 6)

cols = list( synapser::Column( name='ENSG', columnType='STRING', maximumSize=20),
             synapser::Column( name='GName', columnType='STRING', maximumSize=40),
             synapser::Column( name='RNA_TE', columnType='STRING', maximumSize=20),
             synapser::Column( name='RNA_fdr_CorPVal', columnType='STRING', maximumSize=20),
             synapser::Column( name='RNA_Sig', columnType='STRING', maximumSize=20),
             synapser::Column( name='RNA_Direction', columnType='STRING', maximumSize=20),
             synapser::Column( name='RNA_Weight', columnType='STRING', maximumSize=20),
             synapser::Column( name='Pro_TE', columnType='STRING', maximumSize=20),
             synapser::Column( name='Pro_fdr_CorPVal', columnType='STRING', maximumSize=20),
             synapser::Column( name='Pro_Sig', columnType='STRING', maximumSize=20),
             synapser::Column( name='Pro_Direction', columnType='STRING', maximumSize=20),
             synapser::Column( name='Pro_Weight', columnType='STRING', maximumSize=20),
             synapser::Column( name='Harness', columnType='STRING', maximumSize=20),
             synapser::Column( name='Rank', columnType='STRING', maximumSize=20),
             synapser::Column( name='Predicted_Weight', columnType='STRING', maximumSize=20),
             synapser::Column( name='TYPE', columnType='STRING', maximumSize=20),
             synapser::Column( name='Final_Weight', columnType='STRING', maximumSize=20),
             synapser::Column( name='OmicsScore', columnType='STRING', maximumSize=20)
           )

write.csv( Comb, 'OmicsScores.csv', 
           row.names = F, quote = F)

# Pull Project
#project_treatad <- synapser::synGetEntity('syn21532474')

# Pull Omics Score Table
omics_table <- synapser::synTableQuery(sprintf("select * from %s ", 'syn22758536'))
#omics_table <- read.csv(omics_table$filepath)

# Create a snapshot of the Table Entity
body_json <- rjson::toJSON(c(snapshotComment = "Added TMT and SageSeqR final Processed", snapshotLabel = "v5"))
snapshot <-  synapser::synRestPOST(paste0("/entity/", omics_table$tableId, "/table/snapshot"), body = body_json)

# Pull Version 3 of the Omics scores:
omics_table <- synapser::synTableQuery(sprintf("select * from %s ", 'syn22758536'))
scores <- Comb #read.csv(synapser::synGet('syn22758171', version=2)$path)

# Change the entire Table Entity
deleted <- synapser::synDelete(omics_table)
scores <- scores[, c('ENSG', 'GName', 'RNA_TE', 'RNA_fdr_CorPVal', 'RNA_Sig', 'RNA_Direction',
                     'RNA_Weight',	'Pro_TE',	'Pro_fdr_CorPVal',	'Pro_Sig',	'Pro_Direction',
                     'Pro_Weight',	'Harness',	'Rank',	'Predicted_Weight',	'TYPE',	'Final_Weight',
                     'OmicsScore')]
scores$RNA_TE <- signif(scores$RNA_TE,5)
scores$Pro_TE <- signif(scores$Pro_TE,5)

synapser::synStore(synapser::Table(omics_table$tableId, scores))

  
