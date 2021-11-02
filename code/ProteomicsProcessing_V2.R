library(dplyr)
library(doParallel)

#################################################################################
# Load LFQ Data
load(synapser::synGet('syn24828686')$path )
load(synapser::synGet('syn24828683')$path )
load(synapser::synGet('syn24828684')$path )
load(synapser::synGet('syn24828685')$path )

banner_fixmeta <- read.csv(synapser::synGet('syn26403186')$path) 
table(Imputed_Banner$NuMeta$Sample %in% banner_fixmeta$LFQ.ID.Banner)

### Fix the CERAD Score:
temp <- banner_fixmeta[grepl('_',banner_fixmeta$LFQ.ID.Banner),]
row.names(temp) <- temp$LFQ.ID.Banner

Imputed_Banner$NuMeta[row.names(temp) ,]$CERAD <- temp[row.names(Imputed_Banner$NuMeta[row.names(temp) ,]),]$CERAD
# Manual from Eric
Imputed_Banner$NuMeta['b3_041_03' ,]$CERAD <- 2
Imputed_Banner$NuMeta['b4_007_23' ,]$CERAD <- 2
Imputed_Banner$NuMeta['b4_134_04' ,]$CERAD <- 0
rm(temp)
### 

dim(Imputed_BLSA$ScalWins)
dim(Imputed_Mayo$ScalWins)
dim(Imputed_MSBB$ScalWins)
dim(Imputed_Banner$ScalWins)

#################################################################################
#################################################################################
  ### TMT Data
    ### Expression
Proteomics <- c('syn25006659', 'syn25006639', 'syn25006631', 'syn25006652', 'syn25006866', 'syn25006867')
names(Proteomics) <- c(
  'TMT_rosmap_banner_BA9','TMT_emory_BA24', 'TMT_emory_BA9',
  'TMT_Sinai_BA36', 'TMT_rosmap_BA6','TMT_rosmap_BA37')

counts_load_v1 <- function(synid) {
  exp <- read.csv(synapser::synGet(synid)$path, 
                  header = T,
                  row.names = 1)
}
counts_load_v2 <- function(synid) {
  exp <- data.table::fread(synapser::synGet(synid)$path, 
                  header = T,
                  sep = '\t') %>%
    tibble::column_to_rownames(var = "V1")
  }
#exp <- data.table::fread(synapser::synGet(synid)$path, sep = '\t', header = T)

exp <- list(
  'TMT_rosmap_banner_BA9' = NULL, 
  'TMT_emory_BA24' = NULL,
  'TMT_emory_BA9'= NULL ,
  'TMT_Sinai_BA36' = NULL,
  'TMT_rosmap_BA6_37' = NULL
)

exp[[names(Proteomics)[1]]] <- counts_load_v1( as.character(Proteomics[1]) )
exp[[names(Proteomics)[2]]] <- counts_load_v2( as.character(Proteomics[2]) )
exp[[names(Proteomics)[3]]] <- counts_load_v2( as.character(Proteomics[3]) )
exp[[names(Proteomics)[4]]] <- counts_load_v2( as.character(Proteomics[4]) )
exp[[names(Proteomics)[5]]] <- counts_load_v2( as.character(Proteomics[5]) )

# Scale and Wisorize TMT
#DescTools::Winsorize
parallelThreads=7
clusterLocal <- parallel::makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
doParallel::registerDoParallel(clusterLocal)

exp_scaled_wins <- list(
  'TMT_rosmap_banner_BA9' = NULL, 
  'TMT_emory_BA24' = NULL,
  'TMT_emory_BA9'= NULL ,
  'TMT_Sinai_BA36' = NULL,
  'TMT_rosmap_BA6_37' = NULL
)

for( nam in names(exp)){
  exp_scaled_wins[[nam]] <- foreach(i=1:nrow(exp[[nam]]), 
                                    .combine=rbind) %dopar% { scale(DescTools::Winsorize(as.numeric(exp[[nam]][i,]),na.rm = TRUE))[,1] }
}
for( nam in names(exp)){
  row.names(exp_scaled_wins[[nam]]) <- row.names(exp[[nam]])
  colnames(exp_scaled_wins[[nam]]) <- colnames(exp[[nam]])
}

### Brak up the rosmap/banner data
exp_scaled_wins$TMT_banner_BA9 <- exp_scaled_wins$TMT_rosmap_banner_BA9[ ,
  colnames(exp_scaled_wins$TMT_rosmap_banner_BA9 )[grepl( 'banner', colnames(exp_scaled_wins$TMT_rosmap_banner_BA9))] ]
exp_scaled_wins$TMT_rosmap_BA9 <- exp_scaled_wins$TMT_rosmap_banner_BA9[ ,
  colnames(exp_scaled_wins$TMT_rosmap_banner_BA9 )[grepl( 'rosmap', colnames(exp_scaled_wins$TMT_rosmap_banner_BA9))] ]

    ### MetaData
##### BA-9 - Frontal Cortex (FP?)
tmt_ba9_meta <- readxl::read_excel( synapser::synGet('syn25006658')$path, sheet = 'File0.FullTraits', col_names = T) %>%
  #tibble::column_to_rownames(var = "...1")
  tibble::rownames_to_column() %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  tibble::column_to_rownames(var = "1")

excludes <- tmt_ba9_meta[ tmt_ba9_meta$JohnsonDiagnosis ==  'Exclude', ]$SpecimenID

#tmt_ba9_meta[ tmt_ba9_meta$SpecimenID %in% paste0('rosmap.', c("b01.128N", "b02.127C",  "b04.130C",  "b05.130N",  "b11.127C", 
#              "b22.130N", "b26.129N", "b29.130N", "b36.127N", "b38.128N",
#              "b43.127C", "b45.127C")), ]

## Remove outliers and GIS. Filter for Present profiles
table(tmt_ba9_meta$JohnsonDiagnosis)
tmt_ba9_meta <- tmt_ba9_meta[!(tmt_ba9_meta$JohnsonDiagnosis %in% c( 'Exclude','GIS')),]
tmt_ba9_meta <- tmt_ba9_meta[ !(tmt_ba9_meta$Outlier==TRUE),]
tmt_ba9_meta <- tmt_ba9_meta[tmt_ba9_meta$SpecimenID %in% colnames(exp_scaled_wins$TMT_rosmap_banner_BA9),]

############ isolate the BANNER Stats:
tmt_ba9_banner_meta <- Imputed_Banner$NuMeta[ tmt_ba9_meta$LFQ.ID.Banner[tmt_ba9_meta$LFQ.ID.Banner != 'NA'], ]
tmt_ba9_banner_meta$Sample <- tmt_ba9_meta[  tmt_ba9_meta$LFQ.ID.Banner != 'NA', ]$SpecimenID
row.names(tmt_ba9_banner_meta) <- tmt_ba9_banner_meta$Sample


# Add to the meta Data Object
tmt_meta_data <- list( TMT_banner_BA9 = tmt_ba9_banner_meta )

############ isolate the RosMap Stats:
TMT_Express_Load <- function( SynID, names ){
  #'@param  SynID a synapse ID of a proteomics csv matrix eg. syn21266454 or syn21266454
  #'@param  names number corresponding to the row that are the rownames (Use 0 if none apply )
  if( !is.character(SynID) ){
    return("Invalid Synapse ID: Input SynID is not a character string")
  }
  if( !grepl('syn', SynID) | !((nchar(SynID) == 11) | (nchar(SynID) == 10 ))  ){
    return("Invalid Synapse ID: Input SynID is not a valid synapse ID")
  }
  #if(  grepl('Benefactor not found for', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
  #  return("Syn ID does not exist")
  #}
  if(  'You do not have READ permission for the requested entity' %in% as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) ) ){
    return("User does not have READ access")
  }
  if( !grepl('.csv$', synapser::synGet(SynID)$path ) ){
    return("File to import must be a CSV File")
  }
  if( !is.numeric(names) ){
    return("names is not a number")
  }
  
  if( names == 0 ){
    import <- data.frame( data.table::fread( synapser::synGet(SynID)$path, header=T, sep=',' ))
  }else{
    import <- data.frame( data.table::fread( synapser::synGet(SynID)$path, header=T, sep=',' ), row.names = 1)
  }
  return( import )
}


Meta <- TMT_Express_Load('syn21323404', 1)
# - Public Facing BioSpecimin Data: syn21323366
# - Staged BioSpecimin Data: syn23583548
BioSpecimin <- TMT_Express_Load('syn21323366', 0)
if( "assay" %in% colnames(BioSpecimin) ){
}else{
  BioSpecimin <- TMT_Express_Load('syn23583548', 0)
}
BioSpecimin <- BioSpecimin[ BioSpecimin$assay == 'TMT quantitation', ]
Meta$specimenID <- row.names(Meta)
Meta <- dplyr::left_join(Meta, BioSpecimin, by = 'specimenID' )
Clinical <- TMT_Express_Load( 'syn3191087',0 )
Meta <- dplyr::left_join(Meta, Clinical, by = 'individualID' )
Meta <- Meta[ ,colSums(is.na(Meta))<nrow(Meta) ] 
row.names( Meta ) <- Meta$batchChannel
Meta <- Meta[ colnames(Log2_Normalized), ]
Meta <- Meta[ ,colnames(Meta)[ (colnames(Meta) %in%'controlType' )==F] ]
# Harmonize case-control status
Meta$braaksc <- as.numeric(Meta$braaksc)
Meta$ceradsc <- as.numeric(Meta$ceradsc) 
Meta$cogdx <- as.numeric(Meta$cogdx)
# Harmonize case-control status
Meta$diagnosis <- "other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc >= 3] <- "control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc <= 2] <- "AD"
kableExtra::kable( table(Meta$diagnosis) )
## Add Ages over 90 for modeling
Mast <- TMT_Express_Load('syn23573928', 0)
Mast<- Mast[ !duplicated(Mast$individualID), ]
Meta <- dplyr::left_join(Meta[ , colnames(Meta)[ (colnames(Meta) %in% c('age_at_visit_max', 'age_first_ad_dx', 'age_death') )==F ] ], Mast[, c('individualID', 'age_at_visit_max', 'age_first_ad_dx', 'age_death')],by = 'individualID')
#Convert APOE
Meta$apoe_genotype <- as.numeric( Meta$apoe_genotype )
APOS <- as.numeric(names( table(Meta$apoe_genotype) ))
names(APOS) <- APOS
APOS[names( table(Meta$apoe_genotype) )] <- 0
APOS[grepl("4",names(APOS))] <- 1
APOS[grepl("44",names(APOS))] <- 2
Meta$APOE <- as.numeric( APOS[ as.character(Meta$apoe_genotype) ] )

row.names(Meta) <- Meta$batchChannel
#Code Sex and Diagnosis as Factors
Meta$msex <- as.factor(Meta$msex)
Meta$APOE <- as.factor(Meta$APOE)
####Meta for Diagnosis
#Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
Meta_D <- Meta
Meta_D$diagnosis <- as.factor(Meta_D$diagnosis)

# Sample
Meta_D$Sample <- Meta_D$SampleID
# Gender
Meta_D$Gender <- as.numeric(as.character(Meta_D$msex))
Meta_D[ Meta_D$Gender == 0, ]$Gender <- 2
# Age
Meta_D$Age <- Meta_D$age_death
# PMI
Meta_D$PMI <- Meta_D$pmi
# ApoE
Meta_D$ApoE <- Meta_D$APOE
# CERAD
Meta_D$CERAD <- Meta_D$ceradsc
# cogdx
Meta_D$CogDx <- Meta_D$cogdx 
# Braak
Meta_D$Braak <-  Meta_D$braaksc

# Diagnosis
Meta_D$SampleForCombat <- as.character(Meta_D$diagnosis)
Meta_D[Meta_D$SampleForCombat=='other',]$SampleForCombat <- 'Other'
Meta_D[Meta_D$SampleForCombat=='control',]$SampleForCombat <- 'CT'

Meta_D$AD <- 0
Meta_D[ Meta_D$SampleForCombat=='AD', ]$AD <- 1
Meta_D$Control <- 0
Meta_D[ Meta_D$SampleForCombat=='CT', ]$Control <- 1
Meta_D$Other <- 0
Meta_D[ Meta_D$SampleForCombat=='Other', ]$Other <- 1

# Batch 
Meta_D$Batch <- as.character(Meta_D$batch)

Meta_D  <- Meta_D[ !(row.names(Meta_D) %in% gsub('rosmap.', '', excludes)),]

excludes_ros <- excludes[grepl('rosmap.',excludes)]
Meta_D[ !(row.names(Meta_D) %in% gsub('rosmap.','',row.names(excludes_ros))), ]
tmt_meta_data$TMT_rosmap_BA9 <- Meta_D[,c('Sample', 'Gender', 'Age', 'Batch', 
                                          'PMI', 'ApoE', 'CERAD', 'CogDx',
                                          'Braak', 'AD', 'Control', 'Other', 
                                          'SampleForCombat' )]

#### Load the RosMap Expression: deprecated
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)

foo <- matrix( NA, dim(Log2_Normalized)[1], dim(Log2_Normalized)[2] )
colnames(foo) <- colnames(Log2_Normalized)
row.names(foo)<- row.names(Log2_Normalized)

for( i in 1:dim(Log2_Normalized)[1] ){
  foo[i,] <- scale( DescTools::Winsorize( as.numeric(Log2_Normalized[i,]), na.rm = TRUE ) )
}

exp_scaled_wins$TMT_rosmap_BA9 <- foo[,row.names(tmt_meta_data$TMT_rosmap_BA9)]

########  TMT_emory_BA24
TMT_emory_BA24_met <- readxl::read_excel( synapser::synGet('syn25006635')$path, sheet = 'BA24.ADPD', col_names = T) %>%
  tibble::column_to_rownames(var = "batch.channel")

#match to expression
TMT_emory_BA24_met <- TMT_emory_BA24_met[ colnames(exp_scaled_wins$TMT_emory_BA24), ]

# Sample
TMT_emory_BA24_met$Sample <- row.names(TMT_emory_BA24_met)
# Gender
TMT_emory_BA24_met$Gender <- as.numeric(TMT_emory_BA24_met$MSEX)
TMT_emory_BA24_met[ TMT_emory_BA24_met$Gender ==0,]$Gender <- 2

#Age
TMT_emory_BA24_met$Age <- TMT_emory_BA24_met$AGE 
# Batch, PMI             
# ApoE
TMT_emory_BA24_met$ApoE <- TMT_emory_BA24_met$APOE.Risk
TMT_emory_BA24_met[TMT_emory_BA24_met$ApoE==-1,]$ApoE <- 0 

# SampleForCombat
TMT_emory_BA24_met$SampleForCombat <- TMT_emory_BA24_met$Group
TMT_emory_BA24_met[TMT_emory_BA24_met$SampleForCombat=='CTL',]$SampleForCombat <- 'CT'

# CERAD                     
# Braak
TMT_emory_BA24_met$Braak <- TMT_emory_BA24_met$BRAAK

TMT_emory_BA24_met$AD <- 0
TMT_emory_BA24_met[TMT_emory_BA24_met$SampleForCombat=='AD',]$AD <- 1

TMT_emory_BA24_met$ADPD <- 0
TMT_emory_BA24_met[TMT_emory_BA24_met$SampleForCombat=='ADPD',]$ADPD <- 1

TMT_emory_BA24_met$Control <- 0
TMT_emory_BA24_met[TMT_emory_BA24_met$SampleForCombat=='CT',]$Control <- 1

TMT_emory_BA24_met$PD <- 0
TMT_emory_BA24_met[TMT_emory_BA24_met$SampleForCombat=='PD',]$PD <- 1

TMT_emory_BA24_met$CERAD <- as.numeric(TMT_emory_BA24_met$CERAD) + 1
tmt_meta_data$TMT_emory_BA24 <- TMT_emory_BA24_met[,c('Sample', 'Gender', 'Age', 'Batch', 
                                          'PMI', 'ApoE', 'CERAD',
                                          'Braak', 'AD', 'ADPD', 'Control', 'PD', 
                                          'SampleForCombat' )]
ba24_sync <- TMT_emory_BA24_met[,c('Sample', 'Gender', 'Age', 'Batch', 
                      'PMI', 'ApoE', 'CERAD',
                      'Braak', 'AD', 'ADPD', 'Control', 'PD', 
                      'SampleForCombat', 'Emory Case.ID' )]
foo <- tmt_meta_data$TMT_emory_BA24

########  TMT_emory_BA 9
TMT_emory_BA9_met <- readxl::read_excel( synapser::synGet('syn25006627')$path, sheet = 'BA9.ADPD', col_names = T) %>%
  tibble::column_to_rownames(var = "batch.channel")

#match to expression
TMT_emory_BA9_met <- TMT_emory_BA9_met[ colnames(exp_scaled_wins$TMT_emory_BA9), ]

# Sample
TMT_emory_BA9_met$Sample <- row.names(TMT_emory_BA9_met)
# Gender
TMT_emory_BA9_met$Gender <- as.numeric(TMT_emory_BA9_met$MSEX)
TMT_emory_BA9_met[ TMT_emory_BA9_met$Gender ==0,]$Gender <- 2

#Age
TMT_emory_BA9_met$Age <- TMT_emory_BA9_met$AGE 
# Batch, PMI             
# ApoE
TMT_emory_BA9_met$ApoE <- 0
TMT_emory_BA9_met[TMT_emory_BA9_met$APOE.Genotype %in% c('e23', 'e33'),]$ApoE <- 0
TMT_emory_BA9_met[TMT_emory_BA9_met$APOE.Genotype %in% c('e34'),]$ApoE <- 1
TMT_emory_BA9_met[TMT_emory_BA9_met$APOE.Genotype %in% c('e44'),]$ApoE <- 2

#TMT_emory_BA9_met[TMT_emory_BA9_met$ApoE==-1,]$ApoE <- 0 

# SampleForCombat
TMT_emory_BA9_met$SampleForCombat <- TMT_emory_BA9_met$Group
TMT_emory_BA9_met[TMT_emory_BA9_met$SampleForCombat=='CTL',]$SampleForCombat <- 'CT'

# CERAD                     
# Braak
TMT_emory_BA9_met$Braak <- TMT_emory_BA9_met$BRAAK

TMT_emory_BA9_met$AD <- 0
TMT_emory_BA9_met[TMT_emory_BA9_met$SampleForCombat=='AD',]$AD <- 1

TMT_emory_BA9_met$ADPD <- 0
TMT_emory_BA9_met[TMT_emory_BA9_met$SampleForCombat=='ADPD',]$ADPD <- 1

TMT_emory_BA9_met$Control <- 0
TMT_emory_BA9_met[TMT_emory_BA9_met$SampleForCombat=='CT',]$Control <- 1

TMT_emory_BA9_met$PD <- 0
TMT_emory_BA9_met[TMT_emory_BA9_met$SampleForCombat=='PD',]$PD <- 1

TMT_emory_BA9_met$CERAD <- as.numeric(TMT_emory_BA9_met$CERAD) + 1
tmt_meta_data$TMT_emory_BA9 <- TMT_emory_BA9_met[,c('Sample', 'Gender', 'Age', 'Batch', 
                                                      'PMI', 'ApoE', 'CERAD',
                                                      'Braak', 'AD', 'ADPD', 'Control', 'PD', 
                                                      'SampleForCombat' )]
ba9_sync <- TMT_emory_BA9_met[,c('Sample', 'Gender', 'Age', 'Batch', 
                                   'PMI', 'ApoE', 'CERAD',
                                   'Braak', 'AD', 'ADPD', 'Control', 'PD', 
                                   'SampleForCombat', 'Case ID' )]
########
# Sync all Banner Samples:
colnames(ba24_sync)[colnames(ba24_sync)=='Emory Case.ID'] <- 'Case ID'
table(table( c(ba9_sync$`Case ID` ,  ba24_sync$`Case ID`)))
 

table(
  paste0(ba9_sync$`Case ID`, '_', ba9_sync$Gender)  
  %in% 
    paste0( ba24_sync$`Case ID`, '_', ba24_sync$Gender)
)

total_emory <- rbind(ba9_sync, ba24_sync)
total_emory <- total_emory[ ,colnames(total_emory)[!colnames(total_emory) %in% c('Sample','Batch')]]

total_emory <- total_emory[!duplicated(total_emory),]

test <- total_emory[ !duplicated(total_emory$`Case ID`),]
total_emory[ c(total_emory$`Case ID`=='E04-152'),]
########  TMT_sinai_BA 36
TMT_emory_BA36_met <- readxl::read_excel( synapser::synGet('syn25006648')$path, sheet = 'BA36.MSBB', col_names = T) %>%
  tibble::column_to_rownames(var = "batch.channel")

#remove - GIS 
TMT_emory_BA36_met <- TMT_emory_BA36_met[ !grepl('GIS', TMT_emory_BA36_met$SampleID), ]
TMT_emory_BA36_met <- TMT_emory_BA36_met[ !grepl('Exclude', TMT_emory_BA36_met$JohnsonDx), ]
TMT_emory_BA36_met <- TMT_emory_BA36_met[ !(TMT_emory_BA36_met$Outlier==TRUE),]


TMT_emory_BA36_met <- TMT_emory_BA36_met[ colnames(exp_scaled_wins$TMT_Sinai_BA36),]

table(TMT_emory_BA36_met$SampleID %in% Imputed_MSBB$NuMeta$SAMPLE)


# SAMPLE - Sample
TMT_emory_BA36_met$SAMPLE <- TMT_emory_BA36_met$SampleID
TMT_emory_BA36_met$Sample <- row.names(TMT_emory_BA36_met)

# REGION
TMT_emory_BA36_met$REGION <- 'PHG'

# CERAD
TMT_emory_BA36_met$CERAD <- TMT_emory_BA36_met$CERAD.std 
# Braak
TMT_emory_BA36_met$Age <- TMT_emory_BA36_met$Age.Death.Censored

# Gender
TMT_emory_BA36_met$Gender <- as.numeric(TMT_emory_BA36_met$Sex) 
TMT_emory_BA36_met[TMT_emory_BA36_met$Gender == 0, ]$Gender <- 2

# PMI - Minutes*
TMT_emory_BA36_met$PMI <- as.numeric(TMT_emory_BA36_met$PMI.minutes)

# Apo1 Apo2
TMT_emory_BA36_met$Apo1 <- do.call( rbind,strsplit(gsub('e','',TMT_emory_BA36_met$ApoE.Genotype), '[.]'))[,1]
TMT_emory_BA36_met$Apo2 <- do.call( rbind,strsplit(gsub('e','',TMT_emory_BA36_met$ApoE.Genotype), '[.]'))[,2]

# Batch 
TMT_emory_BA36_met$Batch <- TMT_emory_BA36_met$batch

#NP
TMT_emory_BA36_met$NP <- TMT_emory_BA36_met$'NP1 provided as CERAD'

#CDR
#RACE
TMT_emory_BA36_met$RACE <- TMT_emory_BA36_met$Race

# SampleForCombat
TMT_emory_BA36_met$SampleForCombat <- TMT_emory_BA36_met$JohnsonDx

TMT_emory_BA36_met$CDR <- as.numeric(TMT_emory_BA36_met$CDR)
TMT_emory_BA36_met$Braak <- as.numeric(TMT_emory_BA36_met$Braak)
TMT_emory_BA36_met$CERAD <- as.numeric(TMT_emory_BA36_met$CERAD)

TMT_emory_BA36_met$SampleForCombat <- 0
TMT_emory_BA36_met[ TMT_emory_BA36_met$CDR >= 1 &  TMT_emory_BA36_met$Braak >= 4 &  TMT_emory_BA36_met$CERAD >= 2, ]$SampleForCombat <- 'AD'
TMT_emory_BA36_met[ TMT_emory_BA36_met$CDR < 1 &  TMT_emory_BA36_met$Braak < 4 &  TMT_emory_BA36_met$CERAD < 2, ]$SampleForCombat <- 'CT'

#Control
TMT_emory_BA36_met$Control <- 0
TMT_emory_BA36_met[ as.character(TMT_emory_BA36_met$SampleForCombat) == 'CT', ]$Control <- 1
#AD
TMT_emory_BA36_met$AD <- 0
TMT_emory_BA36_met[ as.character(TMT_emory_BA36_met$SampleForCombat) == 'AD', ]$AD <- 1

# bbscore
TMT_emory_BA36_met$bbscore <- TMT_emory_BA36_met$Braak

# Uncensore ages
ages <- as.data.frame(data.table::fread(synapser::synGet('syn10156693')$path))
row.names(ages) <- ages$individualIdentifier
TMT_emory_BA36_met$Age <- ages[TMT_emory_BA36_met$SAMPLE,]$AOD

# Fix the CERAD
MSBB_meta <- read.csv( synapser::synGet( 'syn21893059' )$path, stringsAsFactors = F) %>%
  left_join( read.csv( synapser::synGet( 'syn6101474' )$path, stringsAsFactors = F), by = 'individualID' ) %>%
  left_join( read.csv( synapser::synGet( 'syn22344998' )$path, stringsAsFactors = F ), by = 'specimenID' )

MSBB_meta <- MSBB_meta[ MSBB_meta$assay.x == 'TMT quantitation', ]
MSBB_meta <- MSBB_meta[ !(MSBB_meta$individualID == 'Unknown'), ]
MSBB_meta <- MSBB_meta[ !duplicated(MSBB_meta$individualID),]
row.names(MSBB_meta) <- MSBB_meta$individualID
TMT_emory_BA36_met$CERAD <- MSBB_meta[TMT_emory_BA36_met$SAMPLE,]$CERAD

TMT_emory_BA36_met$Apo1 <- do.call( rbind,strsplit(gsub('','',MSBB_meta[TMT_emory_BA36_met$SAMPLE,]$apoeGenotype), ''))[,1]
TMT_emory_BA36_met$Apo2 <- do.call( rbind,strsplit(gsub('','',MSBB_meta[TMT_emory_BA36_met$SAMPLE,]$apoeGenotype), ''))[,2]


TMT_emory_BA36_met <- TMT_emory_BA36_met[ ,
  c("SAMPLE", "Sample", "REGION", "Control", "AD", "CERAD", "Braak", "Age",
    "Gender", "PMI", "Apo1", "Apo2", "Batch", "RACE", "CDR", "NP", "bbscore", 
    "SampleForCombat")]

Imputed_MSBB$NuMeta[ Imputed_MSBB$NuMeta$SampleForCombat == 'CT',]$AD <- 0
Imputed_MSBB$NuMeta[ Imputed_MSBB$NuMeta$SampleForCombat == 'AD',]$Control <- 0

MSBB_meta <- read.csv( synapser::synGet( 'syn21893059' )$path, stringsAsFactors = F) %>%
  left_join( read.csv( synapser::synGet( 'syn6101474' )$path, stringsAsFactors = F), by = 'individualID' ) %>%
  left_join( read.csv( synapser::synGet( 'syn22344998' )$path, stringsAsFactors = F ), by = 'specimenID' )


MSBB_meta <- MSBB_meta[ as.character(MSBB_meta$assay.x) == 'label free mass spectrometry', ]
MSBB_meta <- MSBB_meta[ !duplicated(MSBB_meta$individualID),]
row.names(MSBB_meta) <- MSBB_meta$individualID

Imputed_MSBB$NuMeta <- Imputed_MSBB$NuMeta[ !is.na(Imputed_MSBB$NuMeta$SAMPLE),]
Imputed_MSBB$NuMeta$CERAD <- MSBB_meta[ Imputed_MSBB$NuMeta$SAMPLE, ]$CERAD
Imputed_MSBB$NuMeta$Braak <- MSBB_meta[ Imputed_MSBB$NuMeta$SAMPLE, ]$Braak
Imputed_MSBB$NuMeta$bbscore <- MSBB_meta[ Imputed_MSBB$NuMeta$SAMPLE, ]$Braak
Imputed_MSBB$NuMeta$CDR <- MSBB_meta[ Imputed_MSBB$NuMeta$SAMPLE, ]$CDR

Imputed_MSBB$NuMeta <- Imputed_MSBB$NuMeta[!grepl('GIS', Imputed_MSBB$NuMeta$SampleForCombat),]
Imputed_MSBB$NuMeta$SampleForCombat <- 0
Imputed_MSBB$NuMeta[ Imputed_MSBB$NuMeta$CDR >= 1 &  Imputed_MSBB$NuMeta$Braak >= 4 &  Imputed_MSBB$NuMeta$CERAD >= 2, ]$SampleForCombat <- 'AD'
Imputed_MSBB$NuMeta[ Imputed_MSBB$NuMeta$CDR < 1 &  Imputed_MSBB$NuMeta$Braak < 4 &  Imputed_MSBB$NuMeta$CERAD < 2, ]$SampleForCombat <- 'CT'

Imputed_MSBB$NuMeta$AD <- 0
Imputed_MSBB$NuMeta[ Imputed_MSBB$NuMeta$SampleForCombat == 'AD',]$AD <- 1
Imputed_MSBB$NuMeta$Control <- 0
Imputed_MSBB$NuMeta[ Imputed_MSBB$NuMeta$SampleForCombat == 'CT',]$Control <- 1

Imputed_MSBB$NuMeta$Apo1 <- do.call(rbind,strsplit(as.character(MSBB_meta[Imputed_MSBB$NuMeta$SAMPLE,]$apoeGenotype), ''))[,1]
Imputed_MSBB$NuMeta$Apo2 <- do.call(rbind,strsplit(as.character(MSBB_meta[Imputed_MSBB$NuMeta$SAMPLE,]$apoeGenotype), ''))[,2]
Imputed_MSBB$NuMeta$Apo1 <- as.numeric(Imputed_MSBB$NuMeta$Apo1)
Imputed_MSBB$NuMeta$Apo2 <- as.numeric(Imputed_MSBB$NuMeta$Apo2)

TMT_emory_BA36_met$Apo1 <- as.numeric(TMT_emory_BA36_met$Apo1)
TMT_emory_BA36_met$Apo2 <- as.numeric(TMT_emory_BA36_met$Apo2)

TMT_emory_BA36_met$SampleForCombat <- 0
TMT_emory_BA36_met[ TMT_emory_BA36_met$CDR >= 1 &  TMT_emory_BA36_met$Braak >= 4 &  TMT_emory_BA36_met$CERAD >= 2, ]$SampleForCombat <- 'AD'
TMT_emory_BA36_met[ TMT_emory_BA36_met$CDR < 1 &  TMT_emory_BA36_met$Braak < 4 &  TMT_emory_BA36_met$CERAD < 2, ]$SampleForCombat <- 'CT'

#Control
TMT_emory_BA36_met$Control <- 0
TMT_emory_BA36_met[ as.character(TMT_emory_BA36_met$SampleForCombat) == 'CT', ]$Control <- 1
#AD
TMT_emory_BA36_met$AD <- 0
TMT_emory_BA36_met[ as.character(TMT_emory_BA36_met$SampleForCombat) == 'AD', ]$AD <- 1


foo <- rbind( TMT_emory_BA36_met[, c("SAMPLE", "Control", "AD", "CERAD", "Braak", "Age",
                             "Gender", "PMI", "Apo1", "Apo2", "RACE", "CDR", "NP", "bbscore", 
                             "SampleForCombat")],
       Imputed_MSBB$NuMeta[,c("SAMPLE", "Control", "AD", "CERAD", "Braak", "Age",
                   "Gender", "PMI", "Apo1", "Apo2", "RACE", "CDR", "NP", "bbscore", 
                   "SampleForCombat")])

foo <- foo[!grepl('GIS', foo$SampleForCombat),]

table(table(foo$SAMPLE))

foo <- foo[ !duplicated(foo),]
table(table(foo$SAMPLE))

head(table(foo$SAMPLE)[table(foo$SAMPLE)>1])
foo[foo$SAMPLE =='AMPAD_MSSM_0000032714',]
Imputed_MSBB$NuMeta[Imputed_MSBB$NuMeta$SAMPLE %in%'AMPAD_MSSM_0000032714',]
TMT_emory_BA36_met[TMT_emory_BA36_met$SAMPLE %in%'AMPAD_MSSM_0000032714',]

tmt_meta_data$TMT_msbb_BA36 <- TMT_emory_BA36_met
#### rosmap BA6+BA37
TMT_rosmap_BA6_37_met <- readxl::read_excel( synapser::synGet('syn25006835')$path, sheet = 1, col_names = T) %>%
  tibble::column_to_rownames(var = "batch.channel")

# Toss GIS
TMT_rosmap_BA6_37_met <- TMT_rosmap_BA6_37_met[ !grepl('GIS', TMT_rosmap_BA6_37_met$Diagnosis), ]
TMT_rosmap_BA6_37_met <- TMT_rosmap_BA6_37_met[ !grepl('Exclude', TMT_rosmap_BA6_37_met$Diagnosis), ]
TMT_rosmap_BA6_37_met <- TMT_rosmap_BA6_37_met[!(TMT_rosmap_BA6_37_met$Outlier.within.Region==TRUE),]

# Remove the 3 missing indv (6 samples) with no ROSMAP indv data
TMT_rosmap_BA6_37_met <- TMT_rosmap_BA6_37_met[!(TMT_rosmap_BA6_37_met$msex)=='NA',]

#Fix the leading zero issue in the projectIDs
TMT_rosmap_BA6_37_met$projid <- stringr::str_pad(TMT_rosmap_BA6_37_met$projid, 8, pad = "0")

#Ammend Ages over 90....
AgesUnCen <- as.data.frame(readxl::read_excel( synapser::synGet('syn18693152')$path, sheet = 'syn_all', col_names = T))
row.names(AgesUnCen) <- AgesUnCen$projid
# samples missing age of death
missing <- TMT_rosmap_BA6_37_met[ !(TMT_rosmap_BA6_37_met$projid %in% AgesUnCen$projid),]$projid
missing <- missing[!duplicated(missing)]

TMT_rosmap_BA6_37_met$age_death <- as.numeric(TMT_rosmap_BA6_37_met$age_death)
TMT_rosmap_BA6_37_met[ TMT_rosmap_BA6_37_met$projid %in% row.names(AgesUnCen), ]$age_death <-
  as.numeric(AgesUnCen[ TMT_rosmap_BA6_37_met[ TMT_rosmap_BA6_37_met$projid %in% row.names(AgesUnCen), ]$projid,]$age_death)

# Dataframe for comparison
ba9_sync <- Meta_D[,c('Sample', 'Gender', 'Age', 'Batch', 
                                          'PMI', 'ApoE', 'CERAD', 'CogDx',
                                          'Braak', 'AD', 'Control', 'Other', 
                                          'SampleForCombat','projid' )]
ba9_sync$projid <- stringr::str_pad(ba9_sync$projid, 8, pad = "0")
ba9_sync <- ba9_sync[,colnames(ba9_sync)[!(colnames(ba9_sync) %in% c('Sample','Batch') )]]


table(table(c(TMT_rosmap_BA6_37_met$projid,ba9_sync$projid)))

# Sample
TMT_rosmap_BA6_37_met$Sample <- row.names(TMT_rosmap_BA6_37_met)

# Gender
TMT_rosmap_BA6_37_met$Gender <- TMT_rosmap_BA6_37_met$msex
TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$Gender == 0, ]$Gender <- 2

# Age
TMT_rosmap_BA6_37_met$Age <- TMT_rosmap_BA6_37_met$age_death

# Batch - no batch
# PMI
TMT_rosmap_BA6_37_met$PMI <- as.numeric(TMT_rosmap_BA6_37_met$pmi)

# ApoE
TMT_rosmap_BA6_37_met$ApoE <- TMT_rosmap_BA6_37_met$apoe_genotype
TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$apoe_genotype %in% c(23,33),]$ApoE <- 0

TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$apoe_genotype %in% c(34),]$ApoE <- 1
TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$apoe_genotype %in% c(44),]$ApoE <- 2

table(TMT_rosmap_BA6_37_met$apoe_genotype)
# CERAD
#table(TMT_rosmap_BA6_37_met$CERAD_std)
#table(ba9_sync$CERAD)
TMT_rosmap_BA6_37_met$CERAD <- TMT_rosmap_BA6_37_met$`ceradsc (ROSMAP nonstandard)`

# CogDx
TMT_rosmap_BA6_37_met$CogDx <- TMT_rosmap_BA6_37_met$cogdx

# Braak
TMT_rosmap_BA6_37_met$Braak <- TMT_rosmap_BA6_37_met$braaksc

#Remove excluded samples:
TMT_rosmap_BA6_37_met <- TMT_rosmap_BA6_37_met[ !(TMT_rosmap_BA6_37_met$Diagnosis %in% 'Exclude'),]

# SampleForCombat
TMT_rosmap_BA6_37_met$SampleForCombat <- 'Other' #TMT_rosmap_BA6_37_met$Diagnosis

# Harmonize case-control status

TMT_rosmap_BA6_37_met$SampleForCombat[
  TMT_rosmap_BA6_37_met$cogdx == 1 & TMT_rosmap_BA6_37_met$braaksc <= 3 & TMT_rosmap_BA6_37_met$ceradsc >= 3
  ] <- "CT"
TMT_rosmap_BA6_37_met$SampleForCombat[
  TMT_rosmap_BA6_37_met$cogdx == 4 & TMT_rosmap_BA6_37_met$braaksc >= 4 & TMT_rosmap_BA6_37_met$ceradsc <= 2
  ] <- "AD"

# AD 
TMT_rosmap_BA6_37_met$AD <- 0
TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$SampleForCombat=='AD',]$AD <- 1

# Control 
TMT_rosmap_BA6_37_met$Control <- 0
TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$SampleForCombat=='CT',]$Control <- 1

# Other
TMT_rosmap_BA6_37_met$Other <- 0
TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$SampleForCombat=='Other',]$Other <- 1

# projid
head(TMT_rosmap_BA6_37_met[,colnames(ba9_sync)])
head(ba9_sync)

TMT_rosmap_BA6_37_met$Gender <- as.numeric(TMT_rosmap_BA6_37_met$Gender)
TMT_rosmap_BA6_37_met$ApoE <- as.numeric(TMT_rosmap_BA6_37_met$ApoE)
TMT_rosmap_BA6_37_met$CERAD <- as.numeric(TMT_rosmap_BA6_37_met$CERAD)
TMT_rosmap_BA6_37_met$CogDx <- as.numeric(TMT_rosmap_BA6_37_met$CogDx)
TMT_rosmap_BA6_37_met$Braak <- as.numeric(TMT_rosmap_BA6_37_met$Braak)

TMT_rosmap_BA6_37_filt <- TMT_rosmap_BA6_37_met[ ,c("Sample", "Gender", "Age", "PMI", "ApoE",
                               "CERAD", "CogDx", "Braak", "AD", "Control", 
                               "Other", "SampleForCombat",'projid'), ]

tmt_meta_data$TMT_rosmap_BA6  <- TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$region == 'BA6', c("Sample", "Gender", "Age", "PMI", "ApoE",
                                                                                             "CERAD", "CogDx", "Braak", "AD", "Control", 
                                                                                             "Other", "SampleForCombat"), ]
tmt_meta_data$TMT_rosmap_BA37 <- TMT_rosmap_BA6_37_met[TMT_rosmap_BA6_37_met$region == 'BA37', c("Sample", "Gender", "Age", "PMI", "ApoE",
                                                                                                "CERAD", "CogDx", "Braak", "AD", "Control", 
                                                                                                "Other", "SampleForCombat"), ]
ba9_sync$ApoE <- as.numeric(as.character(ba9_sync$ApoE))

#Fix Missing Age of deat in BA9
ba9_sync[ba9_sync$projid == '94852511',]$Age <- TMT_rosmap_BA6_37_filt[TMT_rosmap_BA6_37_filt$projid == '94852511',]$Age[1]
tmt_meta_data$TMT_rosmap_BA9[ tmt_meta_data$TMT_rosmap_BA9$Sample=='b09_127C',]$Age <- 
  ba9_sync[ba9_sync$projid == '94852511',]$Age 

# Fix Age Death
ba9_sync[is.na(ba9_sync$Age),]$Age <- as.numeric(AgesUnCen[ba9_sync[is.na(ba9_sync$Age),]$projid,]$age_death)
tmt_meta_data$TMT_rosmap_BA9$Age <- ba9_sync[row.names(tmt_meta_data$TMT_rosmap_BA9),]$Age


foo <- rbind(TMT_rosmap_BA6_37_filt[,colnames(ba9_sync)], ba9_sync)
row.names(foo) <- c(1:dim(foo)[1])
#AOD/PMI need to be rounded:
foo$Age <- signif(foo$Age,5)
foo$PMI <- signif(foo$PMI,4)

foo_2  <- foo[!duplicated(foo),]
sink<-tmt_meta_data$TMT_banner_BA9

######### Fix Scaled Wins...
Proteomics <- c('syn25006659','syn25006659', 'syn25006639', 'syn25006631', 'syn25006652', 'syn25006866', 'syn25006867')
names(Proteomics) <- c(
  'TMT_banner_BA9', 'TMT_rosmap_BA9', 'TMT_emory_BA24', 'TMT_emory_BA9',
  'TMT_msbb_BA36', 'TMT_rosmap_BA6','TMT_rosmap_BA37')


exp <- list(
  'TMT_banner_BA9' = NULL,
  'TMT_rosmap_BA9' = NULL, 
  'TMT_emory_BA24' = NULL,
  'TMT_emory_BA9'  = NULL,
  'TMT_msbb_BA36' = NULL,
  'TMT_rosmap_BA6' = NULL,
  'TMT_rosmap_BA37'  = NULL
)

row.names(tmt_meta_data$TMT_banner_BA9) <- gsub('banner.','', row.names(tmt_meta_data$TMT_banner_BA9))

exp[[names(Proteomics)[1]]] <- counts_load_v1( as.character(Proteomics[1]) )
exp[[names(Proteomics)[1]]] <- exp[[names(Proteomics)[1]]][,colnames(exp[[names(Proteomics)[1]]])[grepl('banner.',colnames(exp[[names(Proteomics)[1]]]))]]
colnames(exp[[names(Proteomics)[1]]]) <- gsub('banner.','',colnames(exp[[names(Proteomics)[1]]]))
tmt_meta_data[[names(Proteomics)[1]]] <- tmt_meta_data[[names(Proteomics)[1]]][colnames(exp[[names(Proteomics)[1]]]), ]

exp[[names(Proteomics)[2]]] <- counts_load_v1( as.character(Proteomics[2]) )
exp[[names(Proteomics)[2]]] <- exp[[names(Proteomics)[2]]][,colnames(exp[[names(Proteomics)[2]]])[grepl('rosmap.',colnames(exp[[names(Proteomics)[2]]]))]]
colnames(exp[[names(Proteomics)[2]]]) <- gsub('rosmap.','',colnames(exp[[names(Proteomics)[2]]]))
tmt_meta_data[[names(Proteomics)[2]]] <- tmt_meta_data[[names(Proteomics)[2]]][colnames(exp[[names(Proteomics)[2]]]), ]

exp[[names(Proteomics)[3]]] <- counts_load_v2( as.character(Proteomics[3]) )
exp[[names(Proteomics)[3]]] <- exp[[names(Proteomics)[3]]][,row.names(tmt_meta_data[[names(Proteomics)[3]]])]

exp[[names(Proteomics)[4]]] <- counts_load_v2( as.character(Proteomics[4]) )
exp[[names(Proteomics)[4]]] <- exp[[names(Proteomics)[4]]][,row.names(tmt_meta_data[[names(Proteomics)[4]]])]

exp[[names(Proteomics)[5]]] <- counts_load_v2( as.character(Proteomics[5]) )
exp[[names(Proteomics)[5]]] <- exp[[names(Proteomics)[5]]][,row.names(tmt_meta_data[[names(Proteomics)[5]]])]

exp[[names(Proteomics)[6]]] <- counts_load_v2( as.character(Proteomics[6]) )
exp[[names(Proteomics)[6]]] <- exp[[names(Proteomics)[6]]][,row.names(tmt_meta_data[[names(Proteomics)[6]]])]

exp[[names(Proteomics)[7]]] <- counts_load_v2( as.character(Proteomics[7]) )
exp[[names(Proteomics)[7]]] <- exp[[names(Proteomics)[7]]][,row.names(tmt_meta_data[[names(Proteomics)[7]]])]

# Scale and Wisorize TMT
#DescTools::Winsorize
parallelThreads=7
clusterLocal <- parallel::makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
doParallel::registerDoParallel(clusterLocal)

exp_scaled_wins <- list(
  'TMT_rosmap_BA9' = NULL, 
  'TMT_banner_BA9' = NULL,
  'TMT_emory_BA24' = NULL,
  'TMT_emory_BA9'  = NULL,
  'TMT_msbb_BA36' = NULL,
  'TMT_rosmap_BA6' = NULL,
  'TMT_rosmap_BA37'  = NULL
)

for( nam in names(exp)){
  exp_scaled_wins[[nam]] <- foreach(i=1:nrow(exp[[nam]]), 
                                    .combine=rbind) %dopar% { scale(DescTools::Winsorize(as.numeric(exp[[nam]][i,]),na.rm = TRUE))[,1] }
}
for( nam in names(exp)){
  row.names(exp_scaled_wins[[nam]]) <- row.names(exp[[nam]])
  colnames(exp_scaled_wins[[nam]]) <- colnames(exp[[nam]])
}

Imputed_Banner$AdjEXP <- Imputed_Banner$AdjEXP[ !(grepl( 'CON__', row.names(Imputed_Banner$AdjEXP))),]
Imputed_Banner$ScalWins <- Imputed_Banner$ScalWins[ !(grepl( 'CON__', row.names(Imputed_Banner$ScalWins))),]

# Push All to Synapse
LFQ <- list(Banner = Imputed_Banner,
            Blsa = Imputed_BLSA,
            Mayo = Imputed_Mayo,
            Msbb = Imputed_MSBB
          )
TMT<- list(
  TMT_banner_BA9 = list(
    AdjEXP = exp$TMT_banner_BA9,
    ScalWins = exp_scaled_wins$TMT_banner_BA9,
    NuMeta = tmt_meta_data$TMT_banner_BA9
  ),
  TMT_rosmap_BA9 = list(
    AdjEXP = exp$TMT_rosmap_BA9,
    ScalWins = exp_scaled_wins$TMT_rosmap_BA9,
    NuMeta = tmt_meta_data$TMT_rosmap_BA9
  ),
  TMT_emory_BA24 = list(
    AdjEXP = exp$TMT_emory_BA24,
    ScalWins = exp_scaled_wins$TMT_emory_BA24,
    NuMeta = tmt_meta_data$TMT_emory_BA24
  ),
  TMT_emory_BA9 = list(
    AdjEXP = exp$TMT_emory_BA9,
    ScalWins = exp_scaled_wins$TMT_emory_BA9,
    NuMeta = tmt_meta_data$TMT_emory_BA9
  ),
  TMT_msbb_BA36 = list(
    AdjEXP = exp$TMT_msbb_BA36,
    ScalWins = exp_scaled_wins$TMT_msbb_BA36,
    NuMeta = tmt_meta_data$TMT_msbb_BA36
  ),
  TMT_rosmap_BA6 = list(
    AdjEXP = exp$TMT_rosmap_BA6,
    ScalWins = exp_scaled_wins$TMT_rosmap_BA6,
    NuMeta = tmt_meta_data$TMT_rosmap_BA6
  ),
  TMT_rosmap_BA37 = list(
    AdjEXP = exp$TMT_rosmap_BA37,
    ScalWins = exp_scaled_wins$TMT_rosmap_BA37,
    NuMeta = tmt_meta_data$TMT_rosmap_BA37
  )
)

save( LFQ,TMT, file = 'Proteomics_Data.RData')

syns_used <- c('syn25006659', 'syn25006659', 'syn25006639', 'syn25006631', 'syn25006652', 
'syn25006866', 'syn25006867','syn18693152', 'syn25006835', 'syn21893059', 'syn6101474', 
'syn22344998', 'syn10156693', 'syn25006648', 'syn25006627','syn25006635',  'syn21266454', 
'syn23573928', 'syn3191087', 'syn23583548', 'syn21323366', 'syn21323404', 'syn25006658', 
'syn26403186', 'syn24828685', 'syn24828684', 'syn24828683', 'syn24828686')


parentId <- 'syn26427434'
activityName = 'New Proteomics Data For Analysis';
activityDescription = 'TMT and LFQ  Data for Omics Scoreing';
thisFileName <- 'ProteomicsProcessing_V2.R'

# Github link
thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))


ENRICH_OBJ <-  synapser::synStore( 
  synapser::File(
    path = 'Proteomics_Data.RData',
    name = 'TMT and LFQ Proteomics Data',
    parentId=parentId 
    ),
  activity = synapser::Activity( 
    used = syns_used, 
    name = activityName, 
    executed = thisFile, 
    description = activityDescription
  )
)

file.remove('Proteomics_Data.RData')
################################################################################
### ENSG to Uniprot|HNGC name


#### Redo ENSG Tag to the Proteomic Data:
ensembl=biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
x <- biomaRt::getBM(
  attributes=c("hgnc_symbol","uniprotsptrembl", "ensembl_gene_id", 'external_gene_name', 'uniprot_gn_symbol'),
  values="*",mart=ensembl
)
x$data_base <- 'tr'
y <- biomaRt::getBM(
  attributes=c("hgnc_symbol","uniprotswissprot",  "ensembl_gene_id", 'external_gene_name','uniprot_gn_symbol'),
  values="*",mart=ensembl
)
y$data_base <- 'sp'

#write.table(x,file="c:/uniprot_sptremblTOhgnc_symbol",sep="\t")
#write.table(y,file="c:/uniprot_swissprotTOhgnc_symbol",sep="\t")

#######################################################################################################################
####  UniProt ID Conversion

#Translation <- read.csv(
#  file = '~/Desktop/Projects/TREAT_AD/Proteomics_Id_Translation/Adapt_Erics_code/cleaned_UniprotToGENESYMBOLApr2015_v03.csv',
#  header = T)
# Load the FASTA that Eric used to process the data
Translation <- read.csv( file = synapser::synGet('syn26428681')$path, header = T)
Translation$ProteomicsID <- paste0(Translation$SYMBOL,'|',Translation$UniprotID)

#Remove  the Translation hgnc column this is incorrectly annotated to UniProt
Translation <- Translation[ , colnames(Translation)[!(colnames(Translation) %in% 'hgnc_symbol')]]
########################################################################################################################
## There are translation entries missing Gene Symbol
Translation$fixed <- 'No'
Translation[Translation$SYMBOL==0,]$fixed <- 'Yes'

## Load a Biomart Object from our RNA-Seq pipeline to harmonize ENSGs to current ENSGs
bm <- read.table( file = synapser::synGet('syn26340639')$path, header = T, sep='\t')

####### Split into those with gene symbol in bioMart (N = 85,394) and those not in Biomart (N = 5,017)
# N = 85,265
uniq_bm_gene_names <- Translation[ Translation$SYMBOL %in% names(table(bm$hgnc_symbol)[table(bm$hgnc_symbol) ==1 ]), ]
for( rows in 1:dim(uniq_bm_gene_names)[1]){
  uniq_bm_gene_names[rows,]$EnsemblGene <- bm[ bm$hgnc_symbol == uniq_bm_gene_names[rows,]$SYMBOL,]$ensembl_gene_id
}

# N = 129
rep_bm_gene_names <- Translation[ Translation$SYMBOL %in% names(table(bm$hgnc_symbol)[table(bm$hgnc_symbol) > 1 ]), ]
rep_bm_gene_names$gene_ensg <- paste0(rep_bm_gene_names$SYMBOL, '|', rep_bm_gene_names$EnsemblGene)

rep_bm_gene_names[ rep_bm_gene_names$EnsemblGene == 'ENSG00000206549',]$EnsemblGene <- 'ENSG00000283706'
rep_bm_gene_names[ rep_bm_gene_names$EnsemblGene == 'ENSG00000116957',]$EnsemblGene <- 'ENSG00000116957' 
  
rep_bm_gene_names[ rep_bm_gene_names$EnsemblGene == '',]$EnsemblGene <- 'ENSG00000228623'
  

new_fasta <- rbind(uniq_bm_gene_names,
                   rep_bm_gene_names[ ,colnames(rep_bm_gene_names)[ !(colnames(rep_bm_gene_names) %in% 'gene_ensg')],]
)


#Not in BioMart Object: (N = 5,017)
notin_bm_gene_names <- Translation[ !(Translation$SYMBOL %in% names(table(bm$hgnc_symbol))), ]

notin_bm_gene_names_fix <- notin_bm_gene_names[
  notin_bm_gene_names$EnsemblGene %in% bm$ensembl_gene_id,
]

for( rows in 1:dim(notin_bm_gene_names_fix)[1]){
  notin_bm_gene_names_fix[rows,]$SYMBOL <- bm[ bm$ensembl_gene_id == notin_bm_gene_names_fix[rows,]$EnsemblGene,]$hgnc_symbol
}

new_fasta <- rbind(new_fasta,notin_bm_gene_names_fix)
new_fasta[ new_fasta$SYMBOL == 'TBCE',]$EnsemblGene <- 'ENSG00000285053'

All_names <- c(row.names(LFQ$Banner$AdjEXP), 
  row.names(LFQ$Blsa$ScalWins),
  row.names(LFQ$Mayo$ScalWins),
  row.names(LFQ$Msbb$ScalWins),
  row.names(TMT$TMT_banner_BA9$ScalWins),
  row.names(TMT$TMT_rosmap_BA9$ScalWins),
  row.names(TMT$TMT_emory_BA24$ScalWins),
  row.names(TMT$TMT_emory_BA9$ScalWins),
  row.names(TMT$TMT_msbb_BA36$ScalWins),
  row.names(TMT$TMT_rosmap_BA6$ScalWins),
  row.names(TMT$TMT_rosmap_BA37$ScalWins)
)

total_names <- All_names[!duplicated(All_names)]

table(total_names %in% new_fasta$ProteomicsID)
table(total_names %in% remainder$ProteomicsID)

Omics_Translation <- new_fasta[,c('ProteomicsID', 'UniprotID', 'SYMBOL', 'EnsemblGene')]

total_names[!(total_names%in%Omics_Translation$ProteomicsID)]

Omics_addons <- as.data.frame(matrix(
  NA,
  length(total_names[!(total_names%in%Omics_Translation$ProteomicsID)]),
  4
))
colnames(Omics_addons) <- colnames(Omics_Translation)
Omics_addons$ProteomicsID <- total_names[!(total_names%in%Omics_Translation$ProteomicsID)]

######### Fix One at a time manually....
Omics_addons[1,2:4] <- c('O00139', 'KIF2A', 'ENSG00000068796')
Omics_addons[2,2:4] <- c('Q9UPP5', 'KIAA1107', NA)
Omics_addons[3,2:4] <- c('Q8NCW5', 'NAXE', 'ENSG00000163382')
Omics_addons[4,2:4] <- c('B4DLN1', 'AC139530.3', 'ENSG00000262660')
Omics_addons[5,2:4] <- c('E9PAV3', 'NACA', 'ENSG00000196531')
Omics_addons[6,2:4] <- c('Q9C0E8', 'LNPK', 'ENSG00000144320')
Omics_addons[7,2:4] <- c('F5H5P2', 'AC011462.1', 'ENSG00000255730')
Omics_addons[8,2:4] <- c('F8W031', 'AC073896.1', 'ENSG00000144785')
Omics_addons[9,2:4] <- c('H3BN98', 'AC138811.2', 'ENSG00000260342')
Omics_addons[10,2:4] <- c('H0YIV9', 'AL021546.1', 'ENSG00000111780')
Omics_addons[11,2:4] <- c('A2RTX5', 'TARSL2', 'ENSG00000185418')
Omics_addons[12,2:4] <- c('O00142', 'TK2', 'ENSG00000166548')
Omics_addons[13,2:4] <- c('A5YM72', 'CARNS1', 'ENSG00000172508')
Omics_addons[14,2:4] <- c('A6NHL2', 'TUBAL3', 'ENSG00000178462')
Omics_addons[15,2:4] <- c('P30711', 'GSTT1', 'ENSG00000277656')
Omics_addons[16,2:4] <- c('A3KMH1', 'VWA8', 'ENSG00000102763')
Omics_addons[17,2:4] <- c('Q5VZK9', 'CARMIL1', 'ENSG00000079691')
Omics_addons[18,2:4] <- c('P01619', 'IGKV3-20', 'ENSG00000239951')
Omics_addons[19,2:4] <- c('H0Y858', 'AC097637.1', 'ENSG00000173366')
Omics_addons[20,2:4] <- c('A5YKK6', 'CNOT1', 'ENSG00000125107')
Omics_addons[21,2:4] <- c('A2RRP1', 'NBAS', 'ENSG00000151779')
Omics_addons[22,2:4] <- c('O60260', 'PRKN', 'ENSG00000185345')
Omics_addons[23,2:4] <- c('P01766', 'IGHV3-13', 'ENSG00000211942')
Omics_addons[24,2:4] <- c('P05067', 'APP', 'ENSG00000142192')
Omics_addons[25,2:4] <- c('A6NHL2', 'TUBAL3', 'ENSG00000178462')
Omics_addons[26,2:4] <- c('O00139', 'KIF2A', 'ENSG00000068796')
Omics_addons[27,2:4] <- c('Q5TGZ0', 'MICOS10', 'ENSG00000173436')
Omics_addons[28,2:4] <- c('O00142', 'TK2', 'ENSG00000166548')
Omics_addons[29,2:4] <- c('A6NHX0', 'CASTOR2', 'ENSG00000274070')
Omics_addons[30,2:4] <- c('A1Z1Q3', 'MACROD2', 'ENSG00000172264')
Omics_addons[31,2:4] <- c('A5PKW4', 'PSD', 'ENSG00000059915')
Omics_addons[32,2:4] <- c('A7KAX9', 'ARHGAP32', 'ENSG00000134909')
Omics_addons[33,2:4] <- c('A5YKK6', 'CNOT1', 'ENSG00000125107')
Omics_addons[34,2:4] <- c('A5YM72', 'CARNS1', 'ENSG00000172508')
Omics_addons[35,2:4] <- c('O60241', 'ADGRB2', 'ENSG00000121753')
Omics_addons[36,2:4] <- c('P55196', 'AFDN', 'ENSG00000130396')
Omics_addons[37,2:4] <- c('Q96D05', 'FAM241B', 'ENSG00000171224')
Omics_addons[38,2:4] <- c('Q96PC5', 'MIA2', 'ENSG00000150527')
Omics_addons[39,2:4] <- c('A6NHR9', 'SMCHD1', 'ENSG00000101596')
Omics_addons[40,2:4] <- c('P55196', 'AFDN', 'ENSG00000130396')
Omics_addons[41,2:4] <- c('Q92890', 'UFD1', 'ENSG00000070010')
Omics_addons[42,2:4] <- c('A6NIZ1', 'AC113404.3', 'ENSG00000254893')
Omics_addons[43,2:4] <- c('G3V599', 'MIA2', 'ENSG00000150527')
Omics_addons[44,2:4] <- c('Q86YV5', 'PRAG1', 'ENSG00000275342')
Omics_addons[45,2:4] <- c('Q9BR61', 'ACBD6', 'ENSG00000230124')
Omics_addons[46,2:4] <- c('A4D161', 'FAM221A', 'ENSG00000188732')
Omics_addons[47,2:4] <- c('A6XGL0', 'YJEFN3', 'ENSG00000250067')
Omics_addons[48,2:4] <- c('A9UHW6', 'MIF4GD', 'ENSG00000125457')
Omics_addons[49,2:4] <- c('P01602', 'IGKV1-5', 'ENSG00000243466')
Omics_addons[50,2:4] <- c('P06312', 'IGKV4-1', 'ENSG00000211598')
Omics_addons[51,2:4] <- c('P01619', 'IGKV3-20', 'ENSG00000239951')
Omics_addons[52,2:4] <- c('P01764', 'IGHV3-23', 'ENSG00000211949')
Omics_addons[53,2:4] <- c('P01619', 'IGKV3-20', 'ENSG00000239951')
Omics_addons[54,2:4] <- c('P80748', 'IGLV3-21', 'ENSG00000211662')
Omics_addons[55,2:4] <- c('P01717', 'IGLV3-25', 'ENSG00000211659')
Omics_addons[56,2:4] <- c('P01714', 'IGLV3-19', 'ENSG00000211663')
Omics_addons[57,2:4] <- c('Q69YL0', 'NCBP2AS2', 'ENSG00000270170')
Omics_addons[58,2:4] <- c('P01762', 'IGHV3-11', 'ENSG00000211941')
Omics_addons[59,2:4] <- c('P01742', 'IGHV1-69', 'ENSG00000211973')
Omics_addons[60,2:4] <- c('P01768', 'IGHV3-30', 'ENSG00000270550')
Omics_addons[61,2:4] <- c('Q674X7', 'KAZN', 'ENSG00000189337')
Omics_addons[62,2:4] <- c('P01703', 'IGLV1-40', 'ENSG00000211653')
Omics_addons[63,2:4] <- c('A5YM72', 'CARNS1', 'ENSG00000172508')
Omics_addons[64,2:4] <- c('P01593', 'IGKV1D-33', 'ENSG00000239975')
Omics_addons[65,2:4] <- c('P06312', 'IGKV4-1', 'ENSG00000211598')
Omics_addons[66,2:4] <- c('P06331', 'IGHV4-34', 'ENSG00000211956')
Omics_addons[67,2:4] <- c('P01742', 'IGHV1-69', 'ENSG00000211973')
Omics_addons[68,2:4] <- c('P23083', 'IGHV1-2', 'ENSG00000211934')
Omics_addons[69,2:4] <- c('P01766', 'IGHV3-13', 'ENSG00000211942')
Omics_addons[70,2:4] <- c('Q641Q2', 'WASHC2A', 'ENSG00000099290')
Omics_addons[71,2:4] <- c('P01780', 'IGHV3-7', 'ENSG00000211938')
Omics_addons[72,2:4] <- c('P01593', 'IGKV1D-33', 'ENSG00000239975')
Omics_addons[73,2:4] <- c('P01619', 'IGKV3-20', 'ENSG00000239951')
Omics_addons[74,2:4] <- c('P01764', 'IGHV3-23', 'ENSG00000211949')
Omics_addons[75,2:4] <- c('A0A0B4J1W5', 'GPR75-ASB3', 'ENSG00000270898')
Omics_addons[76,2:4] <- c('P01611', 'IGKV1D-12', 'ENSG00000278857')
Omics_addons[77,2:4] <- c('P04430', 'IGKV1-16', 'ENSG00000240864')
Omics_addons[78,2:4] <- c('P01619', 'IGKV3-20', 'ENSG00000239951')
Omics_addons[79,2:4] <- c('P01700', 'IGLV1-47', 'ENSG00000211648')
Omics_addons[80,2:4] <- c('P04207', 'IGKV3-15', 'ENSG00000244437')
Omics_addons[81,2:4] <- c('P04434', NA, NA)
Omics_addons[82,2:4] <- c('A5PLN9', 'TRAPPC13', 'ENSG00000113597')
Omics_addons[83,2:4] <- c('P01700', 'IGLV1-47', 'ENSG00000211648')
############################ Finalize: 

proteomics <- as.data.frame(rbind(Omics_Translation,Omics_addons))
unannotated <- remainder[ !(remainder$ProteomicsID %in% proteomics$ProteomicsID), ]

size <- dim(bm[ !(bm$ensembl_gene_id %in% proteomics$EnsemblGene), ])[1]
transcriptomics_only <- data.frame(matrix(NA,size,4))
colnames(transcriptomics_only) <- colnames(proteomics)
transcriptomics_only$EnsemblGene <- bm[ !(bm$ensembl_gene_id %in% proteomics$EnsemblGene), ]$ensembl_gene_id
transcriptomics_only$SYMBOL <- bm[ !(bm$ensembl_gene_id %in% proteomics$EnsemblGene), ]$hgnc_symbol

total_features <- data.frame(rbind(proteomics,transcriptomics_only))

feature_translation <- list(
  Proteomics = proteomics,
  Total_features = total_features,
  Transcriptomics = bm,
  unannotated_proteomics = unannotated
)

#Push to synapse
save(feature_translation, file = 'Feature_Translation.RData')

syns_used <- c('syn26428681', 'syn26340639')

parentId <- 'syn26428679'
activityName = 'Proteomics to Transcriptomics Translation';
activityDescription = 'Contains all entities that translate between proteomics and transcriptomics';
thisFileName <- 'ProteomicsProcessing_V2.R'

# Github link
thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))


ENRICH_OBJ <-  synapser::synStore( 
  synapser::File(
    path = 'Feature_Translation.RData',
    name = 'Proteomics to Transcriptomics ID Translation',
    parentId=parentId 
  ),
  activity = synapser::Activity( 
    used = syns_used, 
    name = activityName, 
    executed = thisFile, 
    description = activityDescription
  )
)

file.remove('Feature_Translation.RData')

