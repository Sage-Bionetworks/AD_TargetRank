#This script takes a user specified Yaml file and creates a markdown file from it
#To run: Rscript code/01-Transcript_Profiler_Initializer.R <Configuration YAML FILE>
#eg. Rscript code/01-Initializer.r configs/Initial_49_Targets.yaml

library(yaml)
library(tidyr)
library(dplyr)
library(plyr)
args <- commandArgs(trailingOnly=TRUE)

#READ IN Config
#configuration<-'configs/Initial_49_Targets.yaml'
configuration <- args[1]
config <- read_yaml(configuration)

setwd( config$filedir )

#Load Synapse and our DE Data:
reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

#DE Results
MSBB_DE <- data.table::fread( syn_temp$get('syn22092611')$path )
colnames(Mayo_DE)[ colnames(Mayo_DE) == 'Tissue.ref' ] <- 'Tissue'
Mayo_DE <- data.table::fread( syn_temp$get('syn22093747')$path )
Rosmap_DE <- data.table::fread( syn_temp$get('syn22097286')$path )
colnames(Rosmap_DE)[ colnames(Rosmap_DE) == 'Region' ] <- 'Tissue'

#Normalized Expression
MSBB_EXP <- as.data.frame(data.table::fread( syn_temp$get('syn22092607')$path ))
row.names(MSBB_EXP) <- MSBB_EXP$ensembl_gene_id
MSBB_EXP <- MSBB_EXP[,colnames(MSBB_EXP) != 'ensembl_gene_id']

Mayo_EXP <- as.data.frame( data.table::fread( syn_temp$get('syn22093746')$path ) )
row.names(Mayo_EXP) <- Mayo_EXP$ensembl_gene_id
Mayo_EXP <- Mayo_EXP[,colnames(Mayo_EXP) != 'ensembl_gene_id']

Rosmap_EXP <- as.data.frame( data.table::fread( syn_temp$get('syn22097285')$path ) )
row.names(Rosmap_EXP) <- Rosmap_EXP$ensembl_gene_id
Rosmap_EXP <- Rosmap_EXP[,colnames(Rosmap_EXP) != 'ensembl_gene_id']

#Covariates
MSBB_Cov <- data.table::fread( syn_temp$get('syn22092576')$path )
row.names(MSBB_Cov) <- MSBB_Cov$SampleID
Mayo_Cov <- data.table::fread( syn_temp$get('syn22093739')$path )
row.names(Mayo_Cov) <- Mayo_Cov$SampleID
Rosmap_Cov <- data.table::fread( syn_temp$get('syn22097275')$path )
row.names(Rosmap_Cov) <- Rosmap_Cov$SampleID

#Get Transcripts for Gene List
Genes <- as.vector( read.table(config$genelistfile)$V1 )

library("biomaRt")
gene_id <- 'ensembl_transcript_id'

if( config$genelisttype == "GNAME" ){
  filters <- 'hgnc_symbol'
  attrs <- c(filters, "ensembl_transcript_id", "ensembl_gene_id")
}else{
  if( config$genelisttype == "ENSG" ){
    filters <- 'ensembl_gene_id'
    attrs <- c(filters, "ensembl_transcript_id", "hgnc_symbol")
  }else{
    warning( "GENE NAME TYPE NOT SPECIFIED CORRECTLY Please Use ENSG or GNAME")
  }
}

host <- 'ensembl.org'
organism <- 'hsa'
attrs <- c(filters, "ensembl_transcript_id", "ensembl_gene_id")
ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host)

ds <- biomaRt::listDatasets(ensembl)[, "dataset"]
ds <- grep(paste0("^", organism), ds, value = TRUE)
if (length(ds) == 0) {
  stop(paste("Mart not found for:", organism))
} else if (length(ds) > 1) {
  message("Found several marts")
  sapply(ds, function(d) message(paste(which(ds == d), d, sep = ": ")))
  n <- readline(paste0("Choose mart (1-", length(ds), ") : "))
  ds <- ds[as.integer(n)]
}
ensembl <- biomaRt::useDataset(ds, mart = ensembl)

Trans <- getBM(filters = 'hgnc_symbol',
               attributes = attrs,
               values = Genes,
               mart = ensembl )


#Filter the Larget tables and write out the tiny tidy data
##Make into temp files and store in temp file directory
features <- Trans$ensembl_transcript_id

#DE Results - store
MSBB_DE <- MSBB_DE[ MSBB_DE$ensembl_gene_id %in% features, ] 
write.table( as.data.frame(MSBB_DE), file = 'Transcript_Profiles/Temp_Files/MSBB_DE.tsv', col.names=T, sep ='\t', row.names = F, quote = F )
Mayo_DE <- Mayo_DE[ Mayo_DE$ensembl_gene_id %in% features, ] 
write.table( as.data.frame(Mayo_DE), file = 'Transcript_Profiles/Temp_Files/Mayo_DE.tsv', col.names=T, sep ='\t', row.names = F, quote = F )
Rosmap_DE <- Rosmap_DE[ Rosmap_DE$ensembl_gene_id %in% features, ] 
write.table( as.data.frame(Rosmap_DE), file = 'Transcript_Profiles/Temp_Files/RosMap_DE.tsv', col.names=T, sep ='\t', row.names = F, quote = F )

#Normalized Expression - store
MSBB_EXP <- MSBB_EXP[ row.names(MSBB_EXP)[ row.names(MSBB_EXP) %in% features ], ]
write.table( as.data.frame(MSBB_EXP), file = 'Transcript_Profiles/Temp_Files/MSBB_Exp.tsv', col.names=T, sep ='\t', row.names = T, quote = F )
Mayo_EXP <- Mayo_EXP[ row.names(Mayo_EXP)[ row.names(Mayo_EXP) %in% features ], ]
write.table( as.data.frame(Mayo_EXP), file = 'Transcript_Profiles/Temp_Files/Mayo_Exp.tsv', col.names=T, sep ='\t', row.names = T, quote = F )
Rosmap_EXP <- Rosmap_EXP[ row.names(Rosmap_EXP)[ row.names(Rosmap_EXP) %in% features ], ]
write.table( as.data.frame(Rosmap_EXP), file = 'Transcript_Profiles/Temp_Files/RosMap_Exp.tsv', col.names=T, sep ='\t', row.names = T, quote = F )

#Write Trans to temp file
write.table( as.data.frame(Trans), file = 'Transcript_Profiles/Temp_Files/Translate.tsv', col.names=T, sep ='\t', row.names = F, quote = F )

## Create a Gene Specific Markdown file for each gene and execute it. Then push to synapse
#CREATE the RMD to RUN
genes <- names( table(Trans$hgnc_symbol) )
#name <- 'ACE'

if( dir.exists( paste0( '~/AD_TargetRank/Transcript_Profiles/', config$runname )) ){
}else{
  system(paste0( 'mkdir ~/AD_TargetRank/Transcript_Profiles/', config$runname ))
}

#Wrangle MetaData for Neuropath data

METADATA_TECH_ID <- 'syn6100548'
METADATA_TECH <- syn_temp$get(METADATA_TECH_ID, version = 7)$path %>%
  read.table(sep=',',header=T, row.names=1, stringsAsFactors = F) %>%
  group_by(sampleIdentifier) %>%
  top_n(1, -rRNA.rate) %>%
  dplyr::select(sampleIdentifier, BrodmannArea, barcode, individualIdentifier, batch, RIN) %>%
  unique %>%
  ddply(.(barcode), .fun = function(x){
    x = dplyr::mutate(x, batch = paste(batch, collapse = ''), RIN = sum(RIN, na.rm = T))  
  }) 

# Get 90+ age
METADATA_AGE_ID = 'syn10156693'
METADATA_AGE <- syn_temp$get(METADATA_AGE_ID, version = 3)$path %>%
  read.table(sep='\t', header=T, stringsAsFactors = F)
# Get clinical metadata
METADATA_CLINICAL_ID <- 'syn6101474'
METADATA_CLINICAL <- read.csv( syn_temp$get(METADATA_CLINICAL_ID)$path,  header = T, stringsAsFactors = F)

row.names(METADATA_CLINICAL) <- METADATA_CLINICAL$individualIdentifier

###Left Here:
colnames(METADATA_CLINICAL)[ colnames(METADATA_CLINICAL) == 'individualID'] <- 'individualIdentifier'
colnames(METADATA_CLINICAL)[ colnames(METADATA_CLINICAL) == 'ageDeath'] <- 'AOD'
row.names(METADATA_TECH) <- METADATA_TECH$sampleIdentifier
METADATA <- dplyr::full_join(METADATA_TECH, METADATA_CLINICAL)

colnames(METADATA)[ colnames(METADATA) %in% 'pmi' ] <- 'PMI'
colnames(METADATA)[ colnames(METADATA) %in% 'Braak' ] <- 'bbscore'
colnames(METADATA)[ colnames(METADATA) %in% 'race' ] <- 'RACE'
colnames(METADATA)[ colnames(METADATA) %in% 'CERAD' ] <- 'NP.1'

METADATA = as_tibble(METADATA, .name_repair = "unique") %>%
  dplyr::filter(!is.na(BrodmannArea), 
                !is.na(CDR), 
                !is.na(PMI), 
                !is.na(RIN),
                !is.na(bbscore),
                !is.na(NP.1),
                !is.na(sampleIdentifier),
                RIN >= 5,
                BrodmannArea %in% c('BM10', 'BM22', 'BM36', 'BM44'))
# Fix metadata
ST <- c('M', 'F')
names(ST) <- c('male', 'female')

METADATA$sex <- as.character( ST[as.character(METADATA$sex)] )
METADATA = METADATA %>%
  dplyr::rename(Sex = sex) %>%
  dplyr::mutate(AOD = gsub('\\+','',AOD),
                Tissue = factor(BrodmannArea,
                                levels = c('BM10', 'BM22', 'BM36', 'BM44'),
                                labels = c('BM10' = 'FP', 'BM22' = 'STG', 
                                           'BM36' = 'PHG', 'BM44' = 'IFG')), 
                RIN2 = RIN^2, 
                Sex = factor(Sex, levels = c('M','F'), labels = c('MALE', 'FEMALE')),
                Tissue1 = Tissue)
rownames(METADATA) = METADATA$sampleIdentifier
# Get harmonised case control definition
METADATA$Diagnosis = 'OTHER'
METADATA$Diagnosis[METADATA$CDR <= 0.5 & METADATA$bbscore <= 3 & METADATA$NP.1 <= 1] = 'CONTROL'
METADATA$Diagnosis[METADATA$CDR >= 1 & METADATA$bbscore >= 4 & METADATA$NP.1 >= 2] = 'AD'
METADATA$Tissue.Diagnosis = paste(METADATA$Tissue, METADATA$Diagnosis, sep = '.')
# Fix missing batch
levels(METADATA$batch)[levels(METADATA$batch) == ''] = 'NoBatch'
# Introduce ApoE4
METADATA$APOE4 = 0

METADATA$Apo1 <- do.call( rbind, strsplit( as.character(METADATA$apoeGenotype), "" ))[,1] 
METADATA$Apo2 <- do.call( rbind, strsplit( as.character(METADATA$apoeGenotype), "" ))[,2] 

METADATA$APOE4[METADATA$Apo1 == 4 | METADATA$Apo2 == 4] = 1
METADATA$APOE4[METADATA$Apo1 == 4 & METADATA$Apo2 == 4] = 2
METADATA$Tissue.APOE4 = paste(METADATA$Tissue, METADATA$APOE4, sep = '.')
# Match covariates to expression data

indToRetain = intersect(METADATA$sampleIdentifier, colnames(MSBB_EXP))
indRemoved = setdiff(colnames(MSBB_EXP), METADATA$sampleIdentifier)
METADATA = as.data.frame(METADATA)
rownames(METADATA) = METADATA$sampleIdentifier
METADATA_MSBB = METADATA[indToRetain,]
write.table( as.data.frame(METADATA_MSBB), file = 'Transcript_Profiles/Temp_Files/MSBB_NeuroPath.tsv', col.names=T, sep ='\t', row.names = F, quote = F )


###MAYO NeuroPath Wrangle
# Get clinical metadata
METADATA_TC_ID <- 'syn3817650'
METADATA_TC <- syn_temp$get(METADATA_TC_ID)$path %>%
  read.table(sep=',',header=T, stringsAsFactors=F) %>%
  tidyr::separate(SampleID, c('Donor_ID', 'Tissue'), sep = '_') %>%
  dplyr::mutate(SampleID = paste(Donor_ID, Tissue, sep = '_'))
# Get clinical metadata
METADATA_CE_ID <- 'syn5223705'
METADATA_CE <- syn_temp$get(METADATA_CE_ID)$path %>%
  read.table(sep=',',header=T, stringsAsFactors=F) %>%
  tidyr::separate(SampleID, c('Donor_ID', 'Tissue')) %>%
  dplyr::mutate(SampleID = paste(Donor_ID, Tissue, sep = '_')) %>%
  dplyr::rename(Gender = Sex, FLOWCELL = Flowcell)
# Get picard metrics from synapse
METADATA_PICARD_CE_ID <- 'syn8698214';
METADATA_PICARD_CE <- syn_temp$get(METADATA_PICARD_CE_ID)$path %>%
  read.table(sep='\t', header=T, stringsAsFactors=F) %>% 
  dplyr::rename(ID = sample)
colnames(METADATA_PICARD_CE) = gsub('AlignmentSummaryMetrics__','',colnames(METADATA_PICARD_CE))
colnames(METADATA_PICARD_CE) = gsub('RnaSeqMetrics__','',colnames(METADATA_PICARD_CE))
METADATA_PICARD_TC_ID <- 'syn8698211';
METADATA_PICARD_TC <- syn_temp$get(METADATA_PICARD_TC_ID)$path %>%
  read.table(sep='\t', header=T, stringsAsFactors=F) %>% 
  dplyr::rename(ID = sample)
colnames(METADATA_PICARD_TC) = gsub('AlignmentSummaryMetrics__','',colnames(METADATA_PICARD_TC))
colnames(METADATA_PICARD_TC) = gsub('RnaSeqMetrics__','',colnames(METADATA_PICARD_TC))

METADATA_PICARD_TC <- METADATA_PICARD_TC[, !duplicated(colnames(METADATA_PICARD_TC))]
METADATA_PICARD_CE <- METADATA_PICARD_CE[, !duplicated(colnames(METADATA_PICARD_CE))]
METADATA_PICARD <- dplyr::full_join( METADATA_PICARD_TC, METADATA_PICARD_CE )

# Merge all metadata
colnames(METADATA_PICARD_TC)[ colnames(METADATA_PICARD_TC) == 'ID'] <- 'SampleID'
colnames(METADATA_PICARD_CE)[ colnames(METADATA_PICARD_CE) == 'ID'] <- 'SampleID'

METADATA = list(METADATA_TC, METADATA_CE) %>%
  data.table::rbindlist(use.names = T, fill = T) %>%
  dplyr::inner_join(list(METADATA_PICARD_TC, METADATA_PICARD_CE) %>% 
  data.table::rbindlist(use.names = T, fill = T)) %>%
  as.data.frame

METADATA = as_tibble(METADATA, .name_repair = "unique") %>%
  dplyr::filter(SampleID %in% colnames(Mayo_EXP)) %>%
  dplyr::filter(!is.na(RIN)) %>%
  dplyr::filter(!is.na(PCT_INTRONIC_BASES)) %>%
  dplyr::filter(!is.na(AgeAtDeath)) %>%
  dplyr::filter(!is.na(Source)) %>%
  dplyr::filter(!(SampleID %in% SAMPLES.EXCLUDE)) %>%
  as.data.frame()
# Fix Tissue
METADATA$Tissue = gsub('CER', 'CBE', METADATA$Tissue)
# Fix AgeAtDeath
METADATA = METADATA %>%
  dplyr::mutate(AgeAtDeath = gsub("_or_above", "", AgeAtDeath))
# Fix Diagnosis
METADATA = METADATA %>%
  dplyr::mutate(Source.Diagnosis = Diagnosis) %>%
  dplyr::mutate(Source.Diagnosis = factor(Source.Diagnosis,
                                          levels = c('AD', 'Control', 'Pathologic Aging', 'PSP'),
                                          labels = c('AD', 'CONTROL', 'PATH_AGE', 'PSP'))) %>%
  dplyr::mutate(Diagnosis = gsub("Pathologic Aging", "OTHER", Diagnosis),
                Diagnosis = gsub("PSP", "OTHER", Diagnosis),
                Diagnosis = factor(Diagnosis, 
                                   levels = c("AD", "Control", "OTHER"),
                                   labels = c('AD', 'CONTROL', 'OTHER'))) %>%
  dplyr::mutate(Tissue.Diagnosis = paste(Tissue, Diagnosis, sep = '.'),
                Tissue.SourceDiagnosis = paste(Tissue, Source.Diagnosis, sep = '.'))
# Fix Gender
METADATA = METADATA %>%
  dplyr::mutate(Sex = Gender,
                Sex = factor(Sex, levels = c('F','M'), labels = c('FEMALE','MALE')))
# Get RIN square
METADATA = METADATA %>%
  dplyr::mutate(RIN2 = RIN^2)
# Add apoe4 genotype (0, 1, 2)
METADATA$APOE4 = 0
METADATA$APOE4[METADATA$ApoE %in% c(24, 34)] = 1
METADATA$APOE4[METADATA$ApoE %in% c(44)] = 2
METADATA$Tissue.APOE4 = paste(METADATA$Tissue, METADATA$APOE4, sep = '.')
# Match covariates to expression data
indToRetain = intersect(METADATA$SampleID, colnames(Mayo_EXP))
indRemoved = setdiff(colnames(Mayo_EXP), METADATA$SampleID)
METADATA = as.data.frame(METADATA)
rownames(METADATA) = METADATA$SampleID
METADATA_Mayo = METADATA[indToRetain,]
write.table( as.data.frame(METADATA_Mayo), file = 'Transcript_Profiles/Temp_Files/Mayo_NeuroPath.tsv', col.names=T, sep ='\t', row.names = F, quote = F )


###RosMap NeuroPath Wrangle
# Get clinical metadata
METADATA.CLINICAL_ID <- 'syn3191087'
METADATA.CLINICAL_OBJ <- syn_temp$get(METADATA.CLINICAL_ID, version = 3)
METADATA.CLINICAL <- read.table(METADATA.CLINICAL_OBJ$path,sep=',',header=T)
# Get clinical metadata with uncensored ages
METADATA.CLINICAL_ID1 <- 'syn7116000'
METADATA.CLINICAL_OBJ1 <- syn_temp$get(METADATA.CLINICAL_ID1, version = 1)
METADATA.CLINICAL1 <- read.table(METADATA.CLINICAL_OBJ1$path,sep=',',header=T)
# Get technical covariates
METADATA.TECH_ID <- 'syn4300313'
METADATA.TECH_OBJ <- syn_temp$get(METADATA.TECH_ID, version = 1)
METADATA.TECH <- read.table(METADATA.TECH_OBJ$path,sep='\t',header=T)
# Get picard metrics from synapse
METADATA.PICARD_ID <- 'syn8698240';
METADATA.PICARD <- syn_temp$get(METADATA.PICARD_ID)$path %>%
  data.table::fread() %>%
  dplyr::rename(Sampleid = sample)
# Fix error in technical covariates data
KEY_ID <- 'syn3382527'
KEY <- syn_temp$get(KEY_ID, version = 7)$path %>%
  read.csv %>% 
  dplyr::filter(!is.na(rnaseq_id)) %>%
  dplyr::select(projid, rnaseq_id) %>%
  tidyr::separate(rnaseq_id, c('a','b','batch'), sep = '_') %>% 
  unite(Sampleid, a, b) %>%
  dplyr::select(-batch) %>%
  unique
# Match technical and clinical covariates
METADATA <- METADATA.TECH %>%
  dplyr::left_join(METADATA.PICARD) %>%
  dplyr::select(-projid) %>%
  dplyr::left_join(KEY) %>%
  dplyr::left_join(METADATA.CLINICAL) %>%
  dplyr::select(-age_first_ad_dx, -age_death, -age_at_visit_max) %>%
  dplyr::left_join(METADATA.CLINICAL1)
# Pick higher quality RIN batch for sample 492_120515
METADATA <- METADATA %>%
  dplyr::group_by(Sampleid) %>%
  dplyr::top_n(1, RINcontinuous)
colnames(METADATA) = gsub('AlignmentSummaryMetrics__','',colnames(METADATA))
colnames(METADATA) = gsub('RnaSeqMetrics__','',colnames(METADATA))

names(METADATA)[13] <- "PF_ALIGNED_BASES_alt"
METADATA <- METADATA %>%
  ungroup %>%
  dplyr::filter(Sampleid %in% colnames(Rosmap_EXP)) %>%
  dplyr::filter(!is.na(cogdx), !is.na(braaksc), !is.na(ceradsc)) %>%
  dplyr::filter(!is.na(RINcontinuous)) %>%
  dplyr::filter(!is.na(PCT_INTRONIC_BASES)) %>%
  dplyr::filter(!is.na(pmi)) %>%
  dplyr::filter(!is.na(age_death)) %>%
  as.data.frame()
# Add harmonised case-control status
METADATA$Diagnosis = 'OTHER'
METADATA$Diagnosis[METADATA$cogdx == 1 & METADATA$braaksc <= 3 & METADATA$ceradsc >= 3] = 'CONTROL'
METADATA$Diagnosis[METADATA$cogdx == 4 & METADATA$braaksc >= 4 & METADATA$ceradsc <= 2] = 'AD'
# Add sex variable 
METADATA$Sex = 'FEMALE'
METADATA$Sex[METADATA$msex == 1] = 'MALE'
# Add apoe4 genotype (0, 1, 2)
METADATA$APOE4 = 0
METADATA$APOE4[METADATA$apoe_genotype %in% c(24, 34)] = 1
METADATA$APOE4[METADATA$apoe_genotype %in% c(44)] = 2
# METADATA$APOE4[is.na(METADATA$apoe_genotype)] = NA
# Get square of RIN
METADATA$RINcontinuous2 = METADATA$RINcontinuous^2
# Match covariates to expression data
indToRetain = intersect(METADATA$Sampleid, colnames(Rosmap_EXP))
removedIDs = setdiff(colnames(Rosmap_EXP), METADATA$Sampleid)
Rosmap_EXP = Rosmap_EXP[,indToRetain]
rownames(METADATA) = METADATA$Sampleid
METADATA_Rosmap = METADATA[indToRetain,]
write.table( as.data.frame(METADATA_Rosmap), file = 'Transcript_Profiles/Temp_Files/Rosmap_NeuroPath.tsv', col.names=T, sep ='\t', row.names = F, quote = F )

parentId = config$parentID;
activityName = 'Transcript Analysis';
activityDescription = 'A transcript analysis of all transcripts across an AMP-AD gene';
thisFileName <- '02-Transcript_Master.Rmd'
thisFileName2 <- '01-Transcript_Profiler_Initializer.R'
# Github link
thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink( repository = thisRepo, repositoryPath=paste0( 'code/',thisFileName ))
thisFile2 <- githubr::getPermlink( repository = thisRepo, repositoryPath=paste0( 'code/',thisFileName2 ))

Syns_Used <- c("syn22092611", "syn22093747", "syn22097286", "syn22092607", "syn22093746", "syn22097285",
               "syn22092576", "syn22093739", "syn22097275", "syn6100548", "syn10156693", "syn6101474",
               "syn3817650", "syn5223705", "syn8698214","syn8698211","syn3191087","syn7116000",
               "syn4300313","syn8698240"
              )
CODE <- syn_temp$store(synapseclient$Folder(name = config$runname, parentId = parentId))

for( name in genes ){
  system( paste0( "sed 's/GENEX/", sub("/", "\\\\/", name), "/g' ", config$filedir, "/code/02-Transcript_Master.Rmd > ", config$filedir, "/code/03-", name,".Rmd"))
  eval(parse(text= paste0( 'rmarkdown::render(\'~/AD_TargetRank/code/03-', 
                           name,
                           '.Rmd\', \'pdf_document\', output_dir = \'~/AD_TargetRank/Transcript_Profiles/',
                           config$runname, '\' )' ) ))
  
  #Load to Synapse and scrub from files
  system( paste0( 'rm ~/AD_TargetRank/code/03-', name, '.Rmd' ) )
  system( paste0( 'rm -r ~/AD_TargetRank/Transcript_Profiles/Initial_49_Targets/03-', name, '_files/' ) )
  
  reticulate::use_python("/usr/bin/python", required = TRUE)
  synapseclient <- reticulate::import("synapseclient")
  syn_temp <- synapseclient$Synapse()
  syn_temp$login()
  
  ENRICH_OBJ <- syn_temp$store( synapseclient$File( path= paste0( '~/AD_TargetRank/Transcript_Profiles/', config$runname,'/03-', name,'.pdf'), name = paste0( 'Transcript analysis for ', name ), 
                                                     parentId=CODE$properties$id ), 
                                                     used = Syns_Used,
                                                     activityName = activityName, 
                                                     executed = c(thisFile,thisFile2),
                                                     activityDescription = activityDescription)
}


