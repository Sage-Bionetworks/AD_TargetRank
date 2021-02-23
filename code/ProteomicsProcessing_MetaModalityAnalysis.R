source('utilityFunctions/loadsyndata.R')
system( 'curl http://www.compgen.pitt.edu/GemTools/GemTools.R > GemTools.R')
source('GemTools.R')

library(biomaRt)
library(doParallel)
library(preprocessCore)
library(doParallel)
library(meta)
library(ggfortify)
library(factoextra)
library( dplyr )

cl <- makeCluster( detectCores()-2 )
registerDoParallel(cl)

reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

Trans <- read.csv( syn_temp$get('syn24216770')$path )
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)
Trans$OldPeptideID <- as.character( Trans$OldPeptideID )
Trans$NewPeptideID <- as.character( Trans$NewPeptideID )
Trans$ENSG <- as.character( Trans$ENSG )

load( syn_temp$get('syn24828732')$path ) 
#head(LFQ_Name_Trans)

LFQ_Total_Name <- as.data.frame( rbind( LFQ_Name_Trans$Banner, LFQ_Name_Trans$Mayo, LFQ_Name_Trans$BLSA, LFQ_Name_Trans$MSBB ))
LFQ_TUniq_Name <- LFQ_Total_Name[ !duplicated(LFQ_Total_Name), ]

LFQ_TUniq_Name$UniqID <- as.character(LFQ_TUniq_Name$UniqID)
LFQ_TUniq_Name$GeneName <- as.character(LFQ_TUniq_Name$GeneName)
LFQ_TUniq_Name$UniProtID <- as.character(LFQ_TUniq_Name$UniProtID)
LFQ_TUniq_Name$UniqID[ LFQ_TUniq_Name$UniqID == 'PALM2|Q8IXS6-2'] <- 'PALM2AKAP2|Q8IXS6-2'
LFQ_TUniq_Name$UniqID[ LFQ_TUniq_Name$UniqID == 'PALM2|Q8IXS6'] <- 'PALM2AKAP2|Q8IXS6'
LFQ_TUniq_Name$GeneName[ LFQ_TUniq_Name$GeneName %in% 'PALM2'] <- 'PALM2AKAP2'

LFQ_TUniq_Name$Old <- 0
LFQ_TUniq_Name$New <- 0

LFQ_TUniq_Name[ as.character(LFQ_TUniq_Name$UniqID) %in% as.character(Trans$OldPeptideID),]$Old <- 1
LFQ_TUniq_Name[ as.character(LFQ_TUniq_Name$UniqID) %in% as.character(Trans$NewPeptideID),]$New <- 1

### Input ENSG into the LFQ DATA
LFQ_TUniq_Name$ENSG <- NA
row.names(Trans) <- as.character(Trans$NewPeptideID)

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniqID %in% row.names(Trans),]$ENSG <- Trans[ LFQ_TUniq_Name[ LFQ_TUniq_Name$UniqID %in% row.names(Trans),]$UniqID, ]$ENSG

##### - #####
LFQ_TUniq_Name$NewUniqID <- LFQ_TUniq_Name$UniqID
LFQ_TUniq_Name$NewGeneName<- LFQ_TUniq_Name$GeneName
LFQ_TUniq_Name$NewUniProtID <- LFQ_TUniq_Name$UniProtID

### Get the ENSGs 
##############################################################################################################################
uniProt <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
Values <- LFQ_TUniq_Name[ !grepl( 'ENSG', LFQ_TUniq_Name$ENSG ), ]$UniProtID

attempt <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name","uniprotswissprot"),filter="uniprotswissprot",values=Values
                          ,mart=uniProt, useCache = FALSE)

attempt[ attempt$uniprotswissprot == "C9JLW8", ]$ensembl_gene_id <- 'ENSG00000225663'
attempt[ attempt$uniprotswissprot == "O43251", ]$ensembl_gene_id <- 'ENSG00000100320'
attempt[ attempt$uniprotswissprot == "O94910", ]$ensembl_gene_id <- 'ENSG00000072071'
attempt[ attempt$uniprotswissprot == "P00338", ]$ensembl_gene_id <- 'ENSG00000134333'
attempt[ attempt$uniprotswissprot == "P00751", ]$ensembl_gene_id <- 'ENSG00000243649'
attempt[ attempt$uniprotswissprot == "P01766", ]$ensembl_gene_id <- 'ENSG00000211942'
attempt[ attempt$uniprotswissprot == "P01857", ]$ensembl_gene_id <- 'ENSG00000211896'
attempt[ attempt$uniprotswissprot == "P01871", ]$ensembl_gene_id <- 'ENSG00000211899'
attempt[ attempt$uniprotswissprot == "P04433", ]$ensembl_gene_id <- 'ENSG00000241351'
attempt[ attempt$uniprotswissprot == "P0DMR1", ]$ensembl_gene_id <- 'ENSG00000179412'
attempt[ attempt$uniprotswissprot == "P23468", ]$ensembl_gene_id <- 'ENSG00000153707'
attempt[ attempt$uniprotswissprot == "P26640", ]$ensembl_gene_id <- 'ENSG00000204394'
attempt[ attempt$uniprotswissprot == "P29692", ]$ensembl_gene_id <- 'ENSG00000104529'
attempt[ attempt$uniprotswissprot == "P30046", ]$ensembl_gene_id <- 'ENSG00000099977'
attempt[ attempt$uniprotswissprot == "P30519", ]$ensembl_gene_id <- 'ENSG00000103415'
attempt[ attempt$uniprotswissprot == "P30613", ]$ensembl_gene_id <- 'ENSG00000143627'
attempt[ attempt$uniprotswissprot == "P49589", ]$ensembl_gene_id <- 'ENSG00000110619'
attempt[ attempt$uniprotswissprot == "P61769", ]$ensembl_gene_id <- 'ENSG00000166710'

attempt[ attempt$uniprotswissprot == "P62805", ]$ensembl_gene_id <- 'ENSG00000197837'
attempt[ attempt$uniprotswissprot == "P62805", ]$external_gene_name <- 'H4-16'


attempt[ attempt$uniprotswissprot == "P62807", ]$ensembl_gene_id <- 'ENSG00000180596'
attempt[ attempt$uniprotswissprot == "P62807", ]$external_gene_name <- 'H2BC4'

attempt[ attempt$uniprotswissprot == "P68431", ]$ensembl_gene_id <- 'ENSG00000275714'
attempt[ attempt$uniprotswissprot == "P68431", ]$external_gene_name <- 'H3C1'

attempt[ attempt$uniprotswissprot == "P84243", ]$ensembl_gene_id <- 'ENSG00000163041'
attempt[ attempt$uniprotswissprot == "P84243", ]$external_gene_name <- 'H3-3A'

attempt[ attempt$uniprotswissprot == "Q01538", ]$ensembl_gene_id <- 'ENSG00000196132'
attempt[ attempt$uniprotswissprot == "Q07157", ]$ensembl_gene_id <- 'ENSG00000104067'
attempt[ attempt$uniprotswissprot == "Q15843", ]$ensembl_gene_id <- 'ENSG00000129559'
attempt[ attempt$uniprotswissprot == "Q16695", ]$ensembl_gene_id <- 'ENSG00000168148'
attempt[ attempt$uniprotswissprot == "Q4G0P3", ]$ensembl_gene_id <- 'ENSG00000157423'

attempt[ attempt$uniprotswissprot == "Q6FI13", ]$ensembl_gene_id <- 'ENSG00000203812'
attempt[ attempt$uniprotswissprot == "Q6FI13", ]$external_gene_name <- 'H2AC18'

attempt[ attempt$uniprotswissprot == "Q6P1N0", ]$ensembl_gene_id <- 'ENSG00000132024'

attempt[ attempt$uniprotswissprot == "Q71DI3", ]$ensembl_gene_id <- 'ENSG00000203811'
attempt[ attempt$uniprotswissprot == "Q71DI3", ]$external_gene_name <- 'H3C14'

attempt[ attempt$uniprotswissprot == "Q7L7L0", ]$ensembl_gene_id <- 'ENSG00000181218'
attempt[ attempt$uniprotswissprot == "Q96MC5", ]$ensembl_gene_id <- 'ENSG00000166780'
attempt[ attempt$uniprotswissprot == "Q99733", ]$ensembl_gene_id <- 'ENSG00000205531'
attempt[ attempt$uniprotswissprot == "Q9H313", ]$ensembl_gene_id <- 'ENSG00000167614'
attempt[ attempt$uniprotswissprot == "Q9NP81", ]$ensembl_gene_id <- 'ENSG00000104835'
attempt[ attempt$uniprotswissprot == "Q9P2T1", ]$ensembl_gene_id <- 'ENSG00000100938'
attempt[ attempt$uniprotswissprot == "Q9UIG5", ]$ensembl_gene_id <- 'ENSG00000204540'
attempt[ attempt$uniprotswissprot == "Q9UII2", ]$ensembl_gene_id <- 'ENSG00000130770'
attempt[ attempt$uniprotswissprot == "Q9UL46", ]$ensembl_gene_id <- 'ENSG00000100911'

attempt <- attempt[ !duplicated(attempt), ]
row.names(attempt) <- attempt$uniprotswissprot

for( i in 1:dim(LFQ_TUniq_Name)[1] ){
  if( LFQ_TUniq_Name$UniProtID[i] %in% row.names(attempt)  ){
    LFQ_TUniq_Name$NewGeneName[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$external_gene_name
    LFQ_TUniq_Name$NewUniProtID[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$uniprotswissprot
    
    LFQ_TUniq_Name$ENSG[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$ensembl_gene_id
  }
}

##############################################################################################################################
Values <- LFQ_TUniq_Name[ !grepl( 'ENSG', LFQ_TUniq_Name$ENSG ), ]$UniProtID

attempt <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", 'uniprot_isoform'),filter="uniprot_isoform",
                          values=Values,
                          mart=uniProt, useCache = FALSE
)
attempt <- attempt[ !duplicated(attempt),]

attempt[ attempt$uniprot_isoform == "O15260-2", ]$ensembl_gene_id <- 'ENSG00000148248'
attempt[ attempt$uniprot_isoform == "O60888-3", ]$ensembl_gene_id <- 'ENSG00000112514'
attempt[ attempt$uniprot_isoform == "O94910-2", ]$ensembl_gene_id <- 'ENSG00000072071'
attempt[ attempt$uniprot_isoform == "P10636-2", ]$ensembl_gene_id <- 'ENSG00000186868'
attempt[ attempt$uniprot_isoform == "P10636-5", ]$ensembl_gene_id <- 'ENSG00000186868'
attempt[ attempt$uniprot_isoform == "P10636-7", ]$ensembl_gene_id <- 'ENSG00000186868'
attempt[ attempt$uniprot_isoform == "P23468-7", ]$ensembl_gene_id <- 'ENSG00000153707'
attempt[ attempt$uniprot_isoform == "P29692-2", ]$ensembl_gene_id <- 'ENSG00000104529'
attempt[ attempt$uniprot_isoform == "P35749-4", ]$ensembl_gene_id <- 'ENSG00000133392'
attempt[ attempt$uniprot_isoform == "P43686-2", ]$ensembl_gene_id <- 'ENSG00000013275'
attempt[ attempt$uniprot_isoform == "P46379-2", ]$ensembl_gene_id <- 'ENSG00000204463'
attempt[ attempt$uniprot_isoform == "P46379-4", ]$ensembl_gene_id <- 'ENSG00000204463'
attempt[ attempt$uniprot_isoform == "P49589-2", ]$ensembl_gene_id <- 'ENSG00000110619'
attempt[ attempt$uniprot_isoform == "P49589-3", ]$ensembl_gene_id <- 'ENSG00000110619'
attempt[ attempt$uniprot_isoform == "Q13045-2", ]$ensembl_gene_id <- 'ENSG00000177731'
attempt[ attempt$uniprot_isoform == "Q13045-3", ]$ensembl_gene_id <- 'ENSG00000177731'
attempt[ attempt$uniprot_isoform == "Q15262-4", ]$ensembl_gene_id <- 'ENSG00000152894'
attempt[ attempt$uniprot_isoform == "Q53H96-2", ]$ensembl_gene_id <- 'ENSG00000104524'
attempt[ attempt$uniprot_isoform == "Q53HC0-2", ]$ensembl_gene_id <- 'ENSG00000119242'
attempt[ attempt$uniprot_isoform == "Q5SQI0-5", ]$ensembl_gene_id <- 'ENSG00000137343'
attempt[ attempt$uniprot_isoform == "Q5SQI0-7", ]$ensembl_gene_id <- 'ENSG00000137343'
attempt[ attempt$uniprot_isoform == "Q5VIR6-2", ]$ensembl_gene_id <- 'ENSG00000141252'
attempt[ attempt$uniprot_isoform == "Q6NUK1-2", ]$ensembl_gene_id <- 'ENSG00000085491'
attempt[ attempt$uniprot_isoform == "Q6UWP2-2", ]$ensembl_gene_id <- 'ENSG00000278535'
attempt[ attempt$uniprot_isoform == "Q6XQN6-3", ]$ensembl_gene_id <- 'ENSG00000147813'
attempt[ attempt$uniprot_isoform == "Q86V88-3", ]$ensembl_gene_id <- 'ENSG00000213920'
attempt[ attempt$uniprot_isoform == "Q8IXJ6-2", ]$ensembl_gene_id <- 'ENSG00000068903'
attempt[ attempt$uniprot_isoform == "Q8NE71-2", ]$ensembl_gene_id <- 'ENSG00000204574'
attempt[ attempt$uniprot_isoform == "Q8WUK0-2", ]$ensembl_gene_id <- 'ENSG00000110536'
attempt[ attempt$uniprot_isoform == "Q92804-2", ]$ensembl_gene_id <- 'ENSG00000270647'
attempt[ attempt$uniprot_isoform == "Q96MC5-2", ]$ensembl_gene_id <- 'ENSG00000166780'
attempt[ attempt$uniprot_isoform == "Q96PV0-4", ]$ensembl_gene_id <- 'ENSG00000197283'
attempt[ attempt$uniprot_isoform == "Q96QZ7-3", ]$ensembl_gene_id <- 'ENSG00000151276'
attempt[ attempt$uniprot_isoform == "Q96T51-2", ]$ensembl_gene_id <- 'ENSG00000176783'
attempt[ attempt$uniprot_isoform == "Q99767-2", ]$ensembl_gene_id <- 'ENSG00000034053'
attempt[ attempt$uniprot_isoform == "Q9H305-3", ]$ensembl_gene_id <- 'ENSG00000089486'
attempt[ attempt$uniprot_isoform == "Q9UBQ5-2", ]$ensembl_gene_id <- 'ENSG00000178982'
attempt[ attempt$uniprot_isoform == "Q9UBS5-2", ]$ensembl_gene_id <- 'ENSG00000204681'
attempt[ attempt$uniprot_isoform == "Q9UHD1-2", ]$ensembl_gene_id <- 'ENSG00000110172'
attempt[ attempt$uniprot_isoform == "Q9UHX1-4", ]$ensembl_gene_id <- 'ENSG00000179950'
attempt[ attempt$uniprot_isoform == "Q9UJ68-2", ]$ensembl_gene_id <- 'ENSG00000175806'
attempt[ attempt$uniprot_isoform == "Q9UMZ2-6", ]$ensembl_gene_id <- 'ENSG00000275066'
attempt[ attempt$uniprot_isoform == "Q9UMZ2-9", ]$ensembl_gene_id <- 'ENSG00000275066'
attempt[ attempt$uniprot_isoform == "Q9Y243-2", ]$ensembl_gene_id <- 'ENSG00000117020'
attempt[ attempt$uniprot_isoform == "Q9Y296-2", ]$ensembl_gene_id <- 'ENSG00000196655'

attempt <- attempt[ !duplicated(attempt),]
row.names(attempt) <- attempt$uniprot_isoform
for( i in 1:dim(LFQ_TUniq_Name)[1] ){
  if( LFQ_TUniq_Name$UniProtID[i] %in% row.names(attempt)  ){
    LFQ_TUniq_Name$NewGeneName[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$external_gene_name
    LFQ_TUniq_Name$NewUniProtID[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$uniprot_isoform
    
    LFQ_TUniq_Name$ENSG[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$ensembl_gene_id
  }
}

##############################################################################################################################
Values <- LFQ_TUniq_Name[ !grepl( 'ENSG', LFQ_TUniq_Name$ENSG ), ]$UniProtID

attempt <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", 'uniprotsptrembl'),filter="uniprotsptrembl",
                          values=Values,
                          mart=uniProt, useCache = FALSE
)
attempt <- attempt[ !duplicated(attempt),]
attempt$Old <- attempt$uniprotsptrembl

attempt[ attempt$uniprotsptrembl == "A0A0A6YYL4", ]$ensembl_gene_id <- 'ENSG00000103426'
attempt[ attempt$uniprotsptrembl == "C9J650", ]$ensembl_gene_id <- 'ENSG00000166922'
attempt[ attempt$uniprotsptrembl == "C9J650", ]$uniprotsptrembl <- 'P05408'

attempt[ attempt$uniprotsptrembl == "G3V1N2", ]$ensembl_gene_id <- 'ENSG00000188536'
attempt[ attempt$uniprotsptrembl == "G3V1N2", ]$external_gene_name <- 'HBA2'
attempt[ attempt$uniprotsptrembl == "G3V1N2", ]$uniprotsptrembl <- 'P69905'

attempt[ attempt$uniprotsptrembl == "G3V5A3", ]$ensembl_gene_id <- 'ENSG00000100605'
attempt[ attempt$uniprotsptrembl == "G3V5A3", ]$external_gene_name <- 'ITPK1'
attempt[ attempt$uniprotsptrembl == "G3V5A3", ]$uniprotsptrembl <- 'Q13572'

attempt[ attempt$uniprotsptrembl == "G3XAH0", ]$ensembl_gene_id <- 'ENSG00000184702'
attempt[ attempt$uniprotsptrembl == "G3XAH0", ]$external_gene_name <- 'SEPT5'
attempt[ attempt$uniprotsptrembl == "G3XAH0", ]$uniprotsptrembl <- 'Q99719'

attempt[ attempt$uniprotsptrembl == "H0Y8D1", ]$ensembl_gene_id <- 'ENSG00000204983'
attempt[ attempt$uniprotsptrembl == "H0Y8D1", ]$uniprotsptrembl <- 'P07477'

attempt[ attempt$uniprotsptrembl == "H0YA31", ]$ensembl_gene_id <- 'ENSG00000157881'
attempt[ attempt$uniprotsptrembl == "H0YA31", ]$uniprotsptrembl <- 'Q9NVE7'

attempt[ attempt$uniprotsptrembl == "Q5JNW7", ]$ensembl_gene_id <- 'ENSG00000204264'
attempt[ attempt$uniprotsptrembl == "Q5JNW7", ]$uniprotsptrembl <-'P28062'

attempt[ attempt$uniprotsptrembl == "Q5JP53", ]$ensembl_gene_id <- 'ENSG00000196230'
attempt[ attempt$uniprotsptrembl == "Q5JP53", ]$uniprotsptrembl <- 'P07437'

attempt <- attempt[ !duplicated(attempt),]
row.names(attempt) <- attempt$Old

for( i in 1:dim(LFQ_TUniq_Name)[1] ){
  if( LFQ_TUniq_Name$UniProtID[i] %in% row.names(attempt)  ){
    LFQ_TUniq_Name$NewGeneName[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$external_gene_name
    LFQ_TUniq_Name$NewUniProtID[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$uniprotsptrembl
    
    LFQ_TUniq_Name$ENSG[i] <-  attempt[ LFQ_TUniq_Name$UniProtID[i], ]$ensembl_gene_id
  }
}


LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID == 'Q5JP53',]
##############################################################################################################################
Values <- LFQ_TUniq_Name[ !grepl( 'ENSG', LFQ_TUniq_Name$ENSG ), ]$UniProtID
foo<-data.frame(Values, stringsAsFactors = F)
foo$Adapted <- do.call( rbind, strsplit( Values,'-' ))[,1]

attempt <- biomaRt::getBM(attributes =c("ensembl_gene_id","external_gene_name", 'uniprotswissprot'),filter="uniprotswissprot",
                          values=foo$Adapted,
                          mart=uniProt, useCache = FALSE
)

# attempt[ attempt$uniprotsptrembl == "Q5JNW7", ]$ensembl_gene_id <- ''

attempt[ attempt$uniprotswissprot == 'P34896', ]$ensembl_gene_id <- 'ENSG00000176974'
attempt[ attempt$uniprotswissprot == 'P46459', ]$ensembl_gene_id <- 'ENSG00000073969'
attempt[ attempt$uniprotswissprot == 'P57772', ]$ensembl_gene_id <- 'ENSG00000132394'
attempt[ attempt$uniprotswissprot == 'Q5SQH8', ]$ensembl_gene_id <- 'ENSG00000204564'
attempt[ attempt$uniprotswissprot == 'Q676U5', ]$ensembl_gene_id <- 'ENSG00000085978'
attempt[ attempt$uniprotswissprot == 'Q8N1B4', ]$ensembl_gene_id <- 'ENSG00000223501'
attempt[ attempt$uniprotswissprot == 'Q8N1K5', ]$ensembl_gene_id <- 'ENSG00000172673'
attempt[ attempt$uniprotswissprot == 'Q8WWY3', ]$ensembl_gene_id <- 'ENSG00000105618'
attempt[ attempt$uniprotswissprot == 'Q96CX3', ]$ensembl_gene_id <- 'ENSG00000186446'
attempt[ attempt$uniprotswissprot == 'Q9UBS5', ]$ensembl_gene_id <- 'ENSG00000204681'
attempt[ attempt$uniprotswissprot == 'Q9UPR0', ]$ensembl_gene_id <- 'ENSG00000154822'

attempt <- attempt[ !duplicated(attempt), ]
row.names(attempt) <- attempt$uniprotswissprot

row.names(foo) <- foo$Values

for( i in 1:dim(LFQ_TUniq_Name)[1] ){
  if( LFQ_TUniq_Name$UniProtID[i] %in% row.names(foo)  ){
    
    if( foo[ LFQ_TUniq_Name$UniProtID[i], ]$Adapted %in% row.names(attempt) ){
      ENSG_Fetch <- attempt[ foo[ LFQ_TUniq_Name$UniProtID[i], ]$Adapted, ]
      NameT <- foo[ LFQ_TUniq_Name$UniProtID[i], ]
      
      
      if( ENSG_Fetch$external_gene_name == LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% NameT$Values, ]$GeneName ){
        LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% NameT$Values, ]$ENSG <- ENSG_Fetch$ensembl_gene_id
        
      }else{
        message(paste0( paste0( LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% NameT$Values, ], collapse = " "),
                        paste0( ENSG_Fetch, collapse = " "),
                        paste0( NameT, collapse = " "),
                        collapse = " " ))
      }
    }else{
    }
    
  }
}
##############################################################################################################################
LFQ_TUniq_Name <- LFQ_TUniq_Name[ !grepl('CON_', LFQ_TUniq_Name$UniqID),]
#P15636
#P15636
#P35908
#P01966

Values <- LFQ_TUniq_Name[ !grepl( 'ENSG', LFQ_TUniq_Name$ENSG ), ]$UniProtID

#LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% Values[1], ]
#LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% Values[1], ]

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P62158", ]$ENSG <- 'ENSG00000198668'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P62158", ]$NewUniProtID <- 'P0DP23'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0C4DGS0", ]$ENSG <- 'ENSG00000184983'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0C4DGS0", ]$NewUniProtID <- 'P56556'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A0MTR1", ]$ENSG <- 'ENSG00000140945'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A0MTR1", ]$NewUniProtID <- 'P55290'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X105", ]$ENSG <- 'ENSG00000182985'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X105", ]$NewUniProtID <- 'Q9BY67'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5HY54", ]$ENSG <- 'ENSG00000196924'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5HY54", ]$NewUniProtID <- 'P21333'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X079", ]$ENSG <- 'ENSG00000211896'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X079", ]$NewUniProtID <- 'P01857'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WXS7", ]$NewUniProtID <- 'O43681'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WXS7", ]$NewGeneName <- 'GET3'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WXS7", ]$ENSG <- 'ENSG00000198356'


LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8NCW5-2", ]$ENSG <- 'ENSG00000163382'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8NCW5-2", ]$NewGeneName <- 'NAXE'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E7ENA2", ]$ENSG <- 'ENSG00000184470'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E7ENA2", ]$NewUniProtID <- 'Q9NNW7'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A075B6K9", ]$ENSG <- 'ENSG00000211677'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A075B6K9", ]$NewUniProtID <- 'P0DOY2'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q58FF6", ]$ENSG <- 'ENSG00000282100'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q58FF6", ]$NewUniProtID <- 'Q58FF6'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P0CW22", ]$ENSG <- 'ENSG00000182774'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P0CW22", ]$NewGeneName <- 'RPS17'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "K7EQW8", ]$ENSG <- 'ENSG00000167460'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "K7EQW8", ]$NewUniProtID <- 'P67936'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "H7BXE5", ]$ENSG <- 'ENSG00000130477'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "H7BXE5", ]$NewUniProtID <- 'Q9UPW8'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9C0E8-2", ]$ENSG <- 'ENSG00000144320'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9C0E8-2", ]$NewGeneName <- 'LNPK'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "X6R5Z6", ]$ENSG <- 'ENSG00000168275'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "X6R5Z6", ]$NewUniProtID <- 'Q5JTJ3'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5T063", ]$ENSG <- 'ENSG00000132423'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5T063", ]$NewUniProtID <- 'Q9NZJ6'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X0C8", ]$ENSG <- 'ENSG00000149260'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X0C8", ]$NewUniProtID <- 'O15484'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "D6W648", ]$ENSG <- 'ENSG00000076826'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "D6W648", ]$NewUniProtID <- 'Q9P1Y5'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q86V81", ]$ENSG <- 'ENSG00000183684'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WZZ5", ]$ENSG <- 'ENSG00000087365'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WZZ5", ]$NewUniProtID <- 'Q13435'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "K7EMQ9", ]$ENSG <- 'ENSG00000178982'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "K7EMQ9", ]$NewUniProtID <- 'Q9UBQ5'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A2RTX5-2", ]$ENSG <- 'ENSG00000185418'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A2RTX5-2", ]$NewGeneName <- 'TARS3'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A2RTX5-2", ]$NewUniProtID <- 'A2RTX5'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "U3KQS7", ]$ENSG <- 'ENSG00000130520'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "U3KQS7", ]$NewUniProtID <- 'Q9Y4Z0'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30479", ]$ENSG <- 'ENSG00000234745'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30479", ]$NewUniProtID <- 'P01889'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A6YYC3", ]$ENSG <- 'ENSG00000269955'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A6YYC3", ]$NewGeneName <- 'FMC1-LUC7L2'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A6YYC3", ]$NewUniProtID <- 'A0A0A6YYJ8'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "B1PS43", ]$ENSG <- 'ENSG00000133392'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "B1PS43", ]$NewUniProtID <- 'P35749'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5VZK9-2", ]$ENSG <- 'ENSG00000079691'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5VZK9-2", ]$NewGeneName <- 'CARMIL1'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5VZK9-2", ]$NewUniProtID <- 'Q5VZK9-2'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WXL8", ]$ENSG <- 'ENSG00000211897'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WXL8", ]$NewUniProtID <- 'P01860'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P04206", ]$ENSG <- 'ENSG00000239951'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P04206", ]$NewGeneName <- 'IGKV3-20'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P04206", ]$NewUniProtID <- 'P01619'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "", ]$ENSG <- ''
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "", ]$NewUniProtID <- ''

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X1V9", ]$ENSG <- 'ENSG00000211592'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X1V9", ]$NewUniProtID <- 'P01834'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q96NW7-2", ]$ENSG <- 'ENSG00000033122'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E7EVH9", ]$ENSG <- 'ENSG00000130021'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E7EVH9", ]$NewGeneName <- 'PUDP'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E7EVH9", ]$NewUniProtID <- 'Q08623'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q13748", ]$ENSG <- 'ENSG00000075886'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q13748", ]$NewUniProtID <- 'P0DPH8'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X1A5", ]$ENSG <- 'ENSG00000124214'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X1A5", ]$NewUniProtID <- 'O95793'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q96RT1-3", ]$ENSG <- 'ENSG00000112851'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q96RT1-3", ]$NewGeneName <- 'ERBIN'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O60260-5", ]$ENSG <- 'ENSG00000185345'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O60260-5", ]$NewGeneName <- 'PRKN'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8TCE6-2", ]$ENSG <- 'ENSG00000119979'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8TCE6-2", ]$NewGeneName <- 'DENND10'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A0MS21", ]$ENSG <- 'ENSG00000177614'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A0A0MS21", ]$NewUniProtID <- 'Q8N414'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E9PC74", ]$ENSG <- 'ENSG00000145191'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "E9PC74", ]$NewUniProtID <- 'Q13144'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O60241-3", ]$ENSG <- 'ENSG00000121753'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O60241-3", ]$NewGeneName <- 'BAI2'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P07305-2", ]$ENSG <- 'ENSG00000189060'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P07305-2", ]$NewGeneNameD <- 'H1-0'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P0CG05", ]$ENSG <- 'ENSG00000211677'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P0CG05", ]$NewUniProtID <- 'P0DOY2'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P10316", ]$ENSG <- 'ENSG00000206503'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P10316", ]$NewUniProtID <- 'P04439'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30042", ]$ENSG <- 'ENSG00000160221'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30042", ]$NewGeneName <- 'GATD3A'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30042", ]$NewUniProtID <- 'P0DPI2'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30493", ]$ENSG <- 'ENSG00000234745'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30493", ]$NewUniProtID <- 'P01889'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P55196-2", ]$ENSG <- 'ENSG00000130396'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P55196-2", ]$NewGeneName <- 'AFDN'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5VTU8", ]$ENSG <- 'ENSG00000180389'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5VTU8", ]$NewGeneName <- 'ATP5F1EP2'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q59GN2", ]$ENSG <- 'ENSG00000214289'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8N350-2", ]$ENSG <- 'ENSG00000099625'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8N350-2", ]$NewGeneName <- 'CBARP'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q96D05-2", ]$ENSG <- 'ENSG00000171224'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q96D05-2", ]$NewGeneName <- 'FAM241B'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9BZG9-2", ]$ENSG <- 'ENSG00000180155'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9BZG9-2", ]$NewUniProtID <- 'P0DP58'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9P1Z9-2", ]$ENSG <- 'ENSG00000197816'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WUF0", ]$ENSG <- 'ENSG00000186529'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WUF0", ]$NewUniProtID <- 'Q08477'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P01860", ]$ENSG <- 'ENSG00000211897'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9Y4D8-4", ]$ENSG <- 'ENSG00000173064'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O15320-9", ]$ENSG <- 'ENSG00000150527'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O15320-9", ]$NewGeneName <- 'MIA2'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O15320-9", ]$NewUniProtID <- 'Q96PC5'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "G8JLE9", ]$ENSG <- 'ENSG00000102081'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "G8JLE9", ]$NewUniProtID <- 'Q06787'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "O15197-2", ]$ENSG <- 'ENSG00000106123'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q59G71", ]$ENSG <- 'ENSG00000079308'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q59G71", ]$NewUniProtID <- 'Q9HBL0'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q16181-2", ]$ENSG <- 'ENSG00000122545'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q16181-2", ]$NewGeneName <- 'SEPTIN7'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q6DRA6", ]$ENSG <- 'ENSG00000220323'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q6DRA6", ]$NewGeneName <- 'H2BC19P'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "F8VXZ7", ]$ENSG <- 'ENSG00000167552'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "F8VXZ7", ]$NewUniProtID <- 'Q71U36'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "K7ENF6", ]$ENSG <- 'ENSG00000165868'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "K7ENF6", ]$NewUniProtID <- 'O43301'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P01834", ]$ENSG <- 'ENSG00000211592'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "M0R1E0", ]$ENSG <- 'ENSG00000167578'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "M0R1E0", ]$NewUniProtID <- 'P61018'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q9BYB0-3", ]$ENSG <- 'ENSG00000251322'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q6S5H4-2", ]$ENSG <- 'ENSG00000278522'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q6S5H4-2", ]$NewGeneName <- 'POTEB3'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X082", ]$ENSG <- 'ENSG00000135387'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X082", ]$NewUniProtID <- 'Q14444'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A2A3N6", ]$ENSG <- 'ENSG00000180764'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WUX8", ]$ENSG <- 'ENSG00000155008'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087WUX8", ]$NewUniProtID <- 'Q6UXV4'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q8NHP1", ]$ENSG <- 'ENSG00000211454'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A6NNH0", ]$ENSG <- 'ENSG00000274070'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A6NNH0", ]$NewGeneName <- 'CASTOR2'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A6NNH0", ]$NewUniProtID <- 'A6NHX0'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P16190", ]$ENSG <- 'ENSG00000206503'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P16190", ]$NewUniProtID <- 'P04439'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30466", ]$ENSG <- 'ENSG00000234745'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P30466", ]$NewUniProtID <- 'P01889'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P55196-1", ]$ENSG <- 'ENSG00000130396'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "P55196-1", ]$NewGeneName <- 'AFDN'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X2C0", ]$ENSG <- 'ENSG00000211899'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "A0A087X2C0", ]$NewUniProtID <- 'P01871'

LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5TGZ0-2", ]$ENSG <- 'ENSG00000173436'
LFQ_TUniq_Name[ LFQ_TUniq_Name$UniProtID %in% "Q5TGZ0-2", ]$NewGeneName <- 'MICOS10'


##############################################################################################################################
#Replace SEPT with SEPTIN
LFQ_TUniq_Name[ grepl('SEPT5', LFQ_TUniq_Name$NewGeneName), ]$NewGeneName <- 'SEPTIN5'
LFQ_TUniq_Name$NewUniqID <- paste0( LFQ_TUniq_Name$NewGeneName, '|', LFQ_TUniq_Name$NewUniProtID)

##Look for overlap in TMT - Gained 193 matches Now 60.72247% Match
table(LFQ_TUniq_Name$NewUniqID %in% Trans$OldPeptideID)
table(LFQ_TUniq_Name$NewUniqID %in% Trans$NewPeptideID)

LFQ_TUniq_Name$NewNew <- 0 
LFQ_TUniq_Name$OldNew <- 0

LFQ_TUniq_Name[ as.character(LFQ_TUniq_Name$NewUniqID) %in% as.character(Trans$OldPeptideID),]$OldNew <- 1
LFQ_TUniq_Name[ as.character(LFQ_TUniq_Name$NewUniqID) %in% as.character(Trans$NewPeptideID),]$NewNew <- 1

#ENSG Match..... 96.10573% 
table(LFQ_TUniq_Name$ENSG %in% Trans$ENSG)
# FALSE  TRUE 
#   221  5454 

##############################################################################################################################
# - # Rename the TMT and the LFQ
Trans[ row.names(Log2_Normalized), ]$NewPeptideID
row.names(Log2_Normalized)[ !(row.names(Log2_Normalized) %in% Trans$NewPeptideID) ]
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)

table(Trans$NewPeptideID %in% row.names( Log2_Normalized ))
table(Trans$OldPeptideID %in% row.names( Log2_Normalized ))

grepl(Q6F5E8)

table(table(LFQ_TUniq_Name$ENSG))
table(table(Trans$ENSG))

LFQ_Trans <- LFQ_TUniq_Name[ , c('UniqID','NewUniqID', 'GeneName', 'UniProtID', 'NewGeneName', 'NewUniProtID', 'ENSG') ]
colnames(LFQ_Trans) <- colnames(Trans)

Tran <- data.frame( rbind( as.matrix(Trans), as.matrix(LFQ_Trans)), stringsAsFactors = F)
table(table(Tran$ENSG))
dim(Tran)
#14492     7
Tran <- Tran[ !duplicated(Tran), ]
dim(Tran)
#11092     7

######## Repeat Old IDs....
Tran[ Tran$OldPeptideID == "ASNA1|A0A087WXS7",]$NewPeptideID <- 'GET3|O43681'
Tran[ Tran$OldPeptideID == "ASNA1|A0A087WXS7",]$New_Pep <- 'O43681'

Tran[ Tran$OldPeptideID == "HIST1H4A|P62805",]$ENSG <- 'ENSG00000197837'
Tran[ Tran$OldPeptideID == "HIST1H4A|P62805",]$NewPeptideID <- 'H4-16|P62805'
Tran[ Tran$OldPeptideID == "HIST1H4A|P62805",]$New_Gene <- 'H4-16'

Tran[ Tran$OldPeptideID == "MARC2|Q969Z3",]$Old_Gene <- 'MARC2'

Tran[ Tran$OldPeptideID == "MARCH5|Q9NX47",]$Old_Gene <- 'MARCH5'

Tran[ Tran$OldPeptideID == "SEPT10|Q9P0V9",]$Old_Gene <- 'SEPT10'
Tran[ Tran$OldPeptideID == "SEPT2|Q15019",]$Old_Gene <- 'SEPT2'
Tran[ Tran$OldPeptideID == "SEPT3|Q9UH03-2",]$Old_Gene <- 'SEPT3'
Tran[ Tran$OldPeptideID == "SEPT4|O43236",]$Old_Gene <- 'SEPT4'

Tran[ Tran$OldPeptideID == "SEPT5|G3XAH0",]$Old_Gene <- 'SEPT5'
Tran[ Tran$OldPeptideID == "SEPT5|G3XAH0",]$New_Pep <-'Q99719'
Tran[ Tran$OldPeptideID == "SEPT5|G3XAH0",]$NewPeptideID <-'SEPTIN5|Q99719'

Tran[ Tran$OldPeptideID == "SEPT5|Q99719",]$Old_Gene <- 'SEPT5'
Tran[ Tran$OldPeptideID == "SEPT7|Q16181",]$Old_Gene <- 'SEPT7'
Tran[ Tran$OldPeptideID == "SEPT7|Q16181-2",]$Old_Gene <- 'SEPT7'
Tran[ Tran$OldPeptideID == "SEPT8|Q92599",]$Old_Gene <- 'SEPT8'

Tran <- Tran[ !duplicated(Tran), ]
row.names(Tran) <- Tran$OldPeptideID
#################################

#Winsorize TMT 
pre_Win_TMT <- Log2_Normalized
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)

foo <- matrix( NA, dim(Log2_Normalized)[1], dim(Log2_Normalized)[2] )
colnames(foo) <- colnames(Log2_Normalized)
row.names(foo)<- row.names(Log2_Normalized)

for( i in 1:dim(Log2_Normalized)[1] ){
  foo[i,] <- scale( DescTools::Winsorize( as.numeric(Log2_Normalized[i,]), na.rm = TRUE ) )
}

for( i in 1:dim(foo)[1] ){
  if( !grepl('SEPT5', row.names( foo )[i]) ){
    row.names( foo )[i] <- Tran[ row.names( foo )[i], ]$NewPeptideID
  }
}

row.names( foo )[ grepl('PALM2|Q8IXS6-2', row.names(foo)) ] <- 'PALM2AKAP2|Q8IXS6-2'
row.names( foo )[ grepl('SEPT5|G3XAH0', row.names(foo)) ] <- 'SEPTIN5|G3XAH0'
row.names( foo )[ grepl('SEPT5|Q99719', row.names(foo)) ] <-'SEPTIN5|Q99719'

#################################
# Load the LFQ:
load( syn_temp$get( 'syn24828686' )$path )
load( syn_temp$get( 'syn24828683' )$path )
load( syn_temp$get( 'syn24828684' )$path )
load( syn_temp$get( 'syn24828685' )$path )

LegacyTest<-row.names(Imputed_BLSA$ScalWins)
ReNameR <- function( Mat ){
  
  for( i in 1:length(row.names( Mat )) ){
    if( "PALM2|Q8IXS6" == row.names( Mat )[i] ){
      row.names( Mat )[i] <- 'PALM2AKAP2|Q8IXS6'
      row.names( Mat )[i] <- Tran[ row.names( Mat )[i] , ]$NewPeptideID
    }else{
      row.names( Mat )[i] <- Tran[ row.names( Mat )[i] , ]$NewPeptideID
    }
  }
  return( Mat )
  
}

rn_BLSA <- ReNameR( Imputed_BLSA$ScalWins )

row.names(Imputed_Banner$ScalWins)[ !( row.names(Imputed_Banner$ScalWins) %in% Tran$OldPeptideID)  ]
Imputed_Banner$ScalWins <- Imputed_Banner$ScalWins[ row.names(Imputed_Banner$ScalWins)[ !grepl('CON__',row.names(Imputed_Banner$ScalWins)) ], ]
row.names( Imputed_Banner$ScalWins )[ grepl('ABETA', row.names( Imputed_Banner$ScalWins )) ] <- 'ABETA|P05067'
row.names( Imputed_Banner$ScalWins )[ grepl('PALM2', row.names( Imputed_Banner$ScalWins )) ] <- 'PALM2AKAP2|Q8IXS6'
rn_Banner <- ReNameR( Imputed_Banner$ScalWins )

row.names( Imputed_Mayo$ScalWins )[ grepl('PALM2', row.names( Imputed_Mayo$ScalWins )) ] <- 'PALM2AKAP2|Q8IXS6'
rn_Mayo <- ReNameR( Imputed_Mayo$ScalWins )
row.names( rn_Mayo )[ grepl('BAI2|O60241-3', row.names( rn_Mayo )) ] <- 'ADGRB2|O60241-3'

row.names( Imputed_MSBB$ScalWins )[ grepl('PALM2', row.names( Imputed_MSBB$ScalWins )) ] <- 'PALM2AKAP2|Q8IXS6-2'
rn_MSBB <- ReNameR( Imputed_MSBB$ScalWins )

row.names(rn_Mayo)[ grepl('H1F0|P07305-2', row.names(rn_Mayo)) ] <- 'H1-0|P07305-2'


Peptides <- c( row.names(rn_BLSA),row.names(rn_Banner), row.names(rn_Mayo), row.names(rn_MSBB), row.names(foo) ) 
length(Peptides[ !duplicated(Peptides) ])

Comb <- data.frame( matrix( NA, length(Peptides[ !duplicated(Peptides) ]), sum( dim( foo )[2], dim( rn_BLSA )[2], dim( rn_Banner )[2], dim( rn_Mayo )[2], dim( rn_MSBB  )[2] ) ))
colnames( Comb ) <- c( colnames(rn_BLSA),colnames(rn_Banner), colnames(rn_Mayo), colnames(rn_MSBB), colnames(foo)  )
## "SEPT5|G3XAH0" "SEPT5|Q99719"
row.names( Comb ) <- Peptides[ !duplicated(Peptides) ]


Comb <- as.matrix(Comb)
Comb[ row.names( rn_MSBB ), colnames(rn_MSBB) ] <- rn_MSBB[ row.names( rn_MSBB ), colnames(rn_MSBB) ] 
Comb[ row.names( rn_Mayo ), colnames(rn_Mayo) ] <- rn_Mayo[ row.names( rn_Mayo ), colnames(rn_Mayo) ] 
Comb[ row.names( rn_Banner ), colnames(rn_Banner) ] <- rn_Banner[ row.names( rn_Banner ), colnames(rn_Banner) ] 
Comb[ row.names( rn_BLSA ), colnames(rn_BLSA) ] <- rn_BLSA[ row.names( rn_BLSA ), colnames(rn_BLSA) ] 
Comb[ row.names( foo ), colnames(foo) ] <- foo[ row.names( foo ), colnames(foo) ] 

#Absolutltey Complete
Comp_Comb <- Comb[ complete.cases(Comb), ]

# Find all feature with more than 50% coverage in each cohort
Keeps <- NULL
for( feat in row.names( Comb ) ){
  if( mean(is.na( Comb[ feat, colnames(rn_BLSA) ] )) < .5 & 
      mean(is.na( Comb[ feat, colnames(rn_MSBB) ] )) < .5 & 
      mean(is.na( Comb[ feat, colnames(rn_Banner) ] )) < .5 & 
      mean(is.na( Comb[ feat, colnames(rn_BLSA) ] )) < .5 & 
      mean(is.na( Comb[ feat, colnames(foo) ] )) < .5){
    Keeps<-c(feat,Keeps)
  }
}
#Greater than 50% in all cohorts
GT50_Comb <- Comb[ Keeps, ]

###################################################################################################################################################
# Align Meta Data for the proteomics data:
## TMT 

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

## LFQ 
##########Load Phenotype Meta Data
#BANNER

Combat_Info <- syn_temp$get('syn18914620')
Combat_Info <- read.csv( Combat_Info$path, header=T, stringsAsFactors = F)
rownames(Combat_Info) <- Combat_Info$Sample

Combat_Info <- cbind( Combat_Info, Diagnosis = 0 )
Combat_Info[ Combat_Info$AD == 1 , ]$Diagnosis <- 'AD'
Combat_Info[ Combat_Info$Control == 1 , ]$Diagnosis <- 'CT'
Combat_Info[ grepl( 'bgis', Combat_Info$Sample) == T , ]$Diagnosis <- 'BGIS'
Combat_Info[ grepl( 'egis', Combat_Info$Sample) == T , ]$Diagnosis <- 'EGIS'


##BLSA
BLSA_Met <- syn_temp$get('syn18914694')
BLSA_Met <- read.csv(file = BLSA_Met$path, row.names = 1, stringsAsFactors = F)
colnames(BLSA_Met)[1:2] <- c( "Sample", "SampleName" )
colnames(BLSA_Met)[6] <- "Control"
colnames(BLSA_Met)[10] <- "Braak"
colnames(BLSA_Met)[11] <- "Age"
colnames(BLSA_Met)[12] <- "Gender"

BLSA_Met <- cbind( BLSA_Met, Diagnosis = 0 )
BLSA_Met[ BLSA_Met$AD == 1 , ]$Diagnosis <- 'AD'
BLSA_Met[ BLSA_Met$Control == 1 , ]$Diagnosis <- 'CT'
BLSA_Met[ BLSA_Met$AsymAD == 1 , ]$Diagnosis <- 'AsymAD'
BLSA_Met <- cbind( BLSA_Met, Batch = 1 )
#_#BLSA_Met <- BLSA_Met[ (row.names(BLSA_Met) %in% colnames(BLSA_RDY))==T, ]

BLSA_data <- syn_temp$get('syn18918327')
BLSA_data <- data.frame( readxl::read_excel(BLSA_data$path, sheet = 2), stringsAsFactors = F)

BLSA_Fdata <- as.data.frame( cbind( 'Sample' = BLSA_data$SanpleName, 
                                    'SampleName' = BLSA_data$`Sample#`, 
                                    'BRC' = BLSA_data$`BRC#`, 
                                    'MFG' = BLSA_data$MFG, 
                                    'Precuneus' = BLSA_data$Precuneus, 
                                    'Control' = BLSA_data$CT, 
                                    'AsymAD' = BLSA_data$AsymAD, 
                                    'AD' = BLSA_data$AD, 
                                    'CERAD' = BLSA_data$CERAD, 
                                    'Braak' = BLSA_data$BRAAK, 
                                    'Age' = BLSA_data$AGE, 
                                    'Gender' = BLSA_data$SEX,  
                                    'PMI' = BLSA_data$PMI, 
                                    'ApoE' = BLSA_data$`ApoE Numeric Score` )
                             , stringsAsFactors = F )
row.names(BLSA_Fdata) <- BLSA_Fdata$Sample
BLSA_Fdata[ 1, ]$Gender <- 1
BLSA_Fdata[ 2, ]$Gender <- 0
BLSA_Fdata <- cbind( BLSA_Fdata, Affected = 0 )
BLSA_Fdata[ BLSA_Fdata$AD == 1 , ]$Affected <- 2
BLSA_Fdata[ BLSA_Fdata$Control == 1 , ]$Affected <- 0
BLSA_Fdata[ BLSA_Fdata$AsymAD == 1 , ]$Affected <- 1
BLSA_Fdata <- cbind( BLSA_Fdata, Diagnosis = 0 )
BLSA_Fdata[ BLSA_Fdata$AD == 1 , ]$Diagnosis <- 'AD'
BLSA_Fdata[ BLSA_Fdata$Control == 1 , ]$Diagnosis <- 'CT'
BLSA_Fdata[ BLSA_Fdata$AsymAD == 1 , ]$Diagnosis <- 'AsymAD'
BLSA_Fdata$REGION <- NA
BLSA_Fdata[ BLSA_Fdata$MFG == 1, ]$REGION <- 'MFG'
BLSA_Fdata[ BLSA_Fdata$Precuneus == 1, ]$REGION <- 'Precuneus'


###MAYO
MAYO_meta <- read.csv( syn_temp$get( 'syn23277389' )$path, stringsAsFactors = F) %>%
  left_join( read.csv( syn_temp$get( 'syn20827192' )$path, stringsAsFactors = F), by = 'individualID' ) %>%
  left_join( read.csv( syn_temp$get( 'syn23474101' )$path, stringsAsFactors = F), by = 'specimenID' )

MAYO_meta <- MAYO_meta[ as.character(MAYO_meta$assay.x) == 'label free mass spectrometry', ]

Mayo_Met <- syn_temp$get('syn18914935')
Mayo_Met <- read.csv( file = Mayo_Met$path, row.names = 1, stringsAsFactors = F)
Mayo_Met <- Mayo_Met[ Mayo_Met$SampleNames %in% colnames(Comb), ]
row.names( Mayo_Met ) <- Mayo_Met$Samples_Simple.Proteomics

MAYO_meta$Batch <- MAYO_meta[ MAYO_meta$specimenID %in% row.names( Mayo_Met ), ]$batch

MAYO_meta <- MAYO_meta[ MAYO_meta$specimenID %in% row.names( Mayo_Met ), ]
MAYO_meta$Cor_SName <- Mayo_Met[ MAYO_meta$specimenID, ]$SampleNames
MAYO_meta$ageDeath <- Mayo_Met[ MAYO_meta$specimenID, ]$AgeAtDeath
MAYO_meta$Diagnosis <- Mayo_Met[ MAYO_meta$specimenID, ]$Diagnosis
MAYO_meta$diagnosis <- Mayo_Met[ MAYO_meta$specimenID, ]$Diagnosis
MAYO_meta$pmi <- Mayo_Met[ MAYO_meta$specimenID, ]$PMI

MAYO_meta$Batch <- Mayo_Met[ MAYO_meta$specimenID, ]$batch

###MSBB
MSBB_meta <- read.csv( syn_temp$get( 'syn21893059' )$path, stringsAsFactors = F) %>%
  left_join( read.csv( syn_temp$get( 'syn6101474' )$path, stringsAsFactors = F), by = 'individualID' ) %>%
  left_join( read.csv( syn_temp$get( 'syn22344998' )$path, stringsAsFactors = F ), by = 'specimenID' )


MSBB_meta <- MSBB_meta[ as.character(MSBB_meta$assay.x) == 'label free mass spectrometry', ]
MSBB_meta$Cor_SName <- paste0( do.call( rbind, strsplit( as.character(MSBB_meta$specimenID), '_'))[,1],
                               "_",
                               do.call( rbind, strsplit( as.character(MSBB_meta$specimenID), '_'))[,3],
                               "_",
                               do.call( rbind, strsplit( as.character(MSBB_meta$specimenID), '_'))[,2]
)


MSBB_meta_LFQ <- MSBB_meta[ as.character(MSBB_meta$Cor_SName) %in% colnames(Comb), ]
MSBB_meta_LFQ$ageDeath <- MSBB_meta_LFQ$Cor_SName

tPhen <- syn_temp$get('syn18914939')
tPhen <- read.csv(file=tPhen$path, stringsAsFactors = F)
row.names(tPhen) <- tPhen$RunName

MSBB_meta_LFQ$ageDeath <- tPhen[ MSBB_meta_LFQ$Cor_SName, ]$AOD
MSBB_meta_LFQ <- MSBB_meta_LFQ[  ,c( 'Cor_SName', 'individualID', 'individualIdSource', 
                                     'specimenIdSource', 'tissue', 'BrodmannArea', 'batch', 
                                     'CERAD', 'Braak', 'CDR', 'plaqueMean', 'apoeGenotype', 'pmi', 'ageDeath', 'sex', 'race', 'ethnicity') ]

#############
### Comb LFQ MetaData
####### - PMI
MSBB_meta_LFQ$pmi <- MSBB_meta_LFQ$pmi/60
colnames( Combat_Info)[ colnames( Combat_Info) == 'PMI' ] <- 'pmi'
BLSA_Fdata$pmi <- BLSA_Fdata$PMI


####### - APOE4
Combat_Info$APOE4 <- Combat_Info$ApoE
Combat_Info[ Combat_Info$APOE4 == -1 & (is.na(Combat_Info$APOE4)==F), ]$APOE4 <- 0
Combat_Info[ Combat_Info$APOE4 == -2 & (is.na(Combat_Info$APOE4)==F), ]$APOE4 <- 0
Combat_Info[ Combat_Info$APOE4 == .3 & (is.na(Combat_Info$APOE4)==F), ]$APOE4 <- 0

MSBB_meta_LFQ$APOE4 <- NA
MSBB_meta_LFQ$apoeGenotype <- as.numeric(MSBB_meta_LFQ$apoeGenotype)
MSBB_meta_LFQ[ MSBB_meta_LFQ$apoeGenotype == 22 &  !is.na(MSBB_meta_LFQ$apoeGenotype) , ]$APOE4 <- 0
MSBB_meta_LFQ[ MSBB_meta_LFQ$apoeGenotype == 23 &  !is.na(MSBB_meta_LFQ$apoeGenotype) , ]$APOE4 <- 0
MSBB_meta_LFQ[ MSBB_meta_LFQ$apoeGenotype == 33 &  !is.na(MSBB_meta_LFQ$apoeGenotype) , ]$APOE4 <- 0
MSBB_meta_LFQ[ MSBB_meta_LFQ$apoeGenotype == 24 &  !is.na(MSBB_meta_LFQ$apoeGenotype) , ]$APOE4 <- 1
MSBB_meta_LFQ[ MSBB_meta_LFQ$apoeGenotype == 34 &  !is.na(MSBB_meta_LFQ$apoeGenotype) , ]$APOE4 <- 1
MSBB_meta_LFQ[ MSBB_meta_LFQ$apoeGenotype == 44 &  !is.na(MSBB_meta_LFQ$apoeGenotype) , ]$APOE4 <- 2

MAYO_meta$APOE4 <- NA
MAYO_meta$apoeGenotype <- as.numeric(MAYO_meta$apoeGenotype)
MAYO_meta[ MAYO_meta$apoeGenotype == 23 &  !is.na(MAYO_meta$apoeGenotype) , ]$APOE4 <- 0
MAYO_meta[ MAYO_meta$apoeGenotype == 33 &  !is.na(MAYO_meta$apoeGenotype) , ]$APOE4 <- 0
MAYO_meta[ MAYO_meta$apoeGenotype == 24 &  !is.na(MAYO_meta$apoeGenotype) , ]$APOE4 <- 1
MAYO_meta[ MAYO_meta$apoeGenotype == 34 &  !is.na(MAYO_meta$apoeGenotype) , ]$APOE4 <- 1
MAYO_meta[ MAYO_meta$apoeGenotype == 44 &  !is.na(MAYO_meta$apoeGenotype) , ]$APOE4 <- 2

BLSA_Fdata$APOE4 <- NA
BLSA_Met <- syn_temp$get('syn18914694')
BLSA_Met <- read.csv(file = BLSA_Met$path, row.names = 1, stringsAsFactors = F)
BLSA_Fdata <- BLSA_Fdata[ BLSA_Fdata$Sample %in% colnames(Comb),]
row.names(BLSA_Met) <- BLSA_Met$SampleName
BLSA_Fdata$APOE4 <- BLSA_Met[ BLSA_Fdata$Sample, ]$ApoE
BLSA_Fdata$APOE4[ BLSA_Fdata$APOE4 < 0 & !is.na( BLSA_Fdata$APOE4)] <- 0

Meta_D$APOE4 <- Meta_D$APOE



###### - Tissue
MSBB_meta_LFQ$tissue <- as.character( MSBB_meta_LFQ$tissue )
MSBB_meta_LFQ$tissue <- 'PFC'
MAYO_meta$tissue <- 'TCX'
BLSA_Fdata$tissue <- BLSA_Fdata$REGION
Meta_D$tissue <- 'DLPFC'
Combat_Info$tissue <- 'DLPFC'

####### - Study
MSBB_meta_LFQ$study <- 'MSBB'
MAYO_meta$study <- 'Mayo'
BLSA_Fdata$study <- 'BLSA'
Meta_D$study <- Meta_D$Study
Combat_Info$study <- 'Banner'

####### - Batch
MSBB_meta_LFQ$batch <- paste0( MSBB_meta_LFQ$batch, '_', MSBB_meta_LFQ$study )
MAYO_meta$batch <- paste0( MAYO_meta$Batch, '_', MAYO_meta$study )
BLSA_Fdata$batch <- paste0( 'NA_', BLSA_Fdata$study )
Combat_Info$batch <- paste0( Combat_Info$Batch, '_', Combat_Info$study )
Meta_D$batch <- paste0( Meta_D$batch, '_', Combat_Info$study )

####### - SampleID
BLSA_Fdata$SampleID <- BLSA_Fdata$Sample
Combat_Info$SampleID <- Combat_Info$Sample
MSBB_meta_LFQ$SampleID <- MSBB_meta_LFQ$Cor_SName
MAYO_meta$SampleID <- MAYO_meta$Cor_SName
Meta_D$SampleID <- row.names( Meta_D )

####### - Sex
BLSA_Fdata$Sex <- 'FEMALE'
BLSA_Fdata[ BLSA_Fdata$Gender == 1,]$Sex <- 'MALE'

Combat_Info$Sex <- 'FEMALE'
Combat_Info[ Combat_Info$Gender == 1,]$Sex <- 'MALE'

MSBB_meta_LFQ$Sex <- 'FEMALE'
MSBB_meta_LFQ[ MSBB_meta_LFQ$sex %in% 'male',]$Sex <- 'MALE'

MAYO_meta$Sex <- 'FEMALE'
MAYO_meta[ MAYO_meta$sex %in% 'male',]$Sex <- 'MALE'

Meta_D$Sex <- 'FEMALE'
Meta_D[ Meta_D$msex == 1,]$Sex <- 'MALE'

####### - Diagnosis 
BLSA_Fdata$Diagnosis <- BLSA_Fdata$Diagnosis
MAYO_meta[ MAYO_meta$Diagnosis %in% 'Control', ]$Diagnosis <- 'CT'
Combat_Info$Diagnosis <- Combat_Info$Diagnosis
MSBB_meta_LFQ$Diagnosis <- NA
MSBB_meta_LFQ[ MSBB_meta_LFQ$CDR >= 1 &  MSBB_meta_LFQ$Braak >= 4 &  MSBB_meta_LFQ$CERAD >= 2, ]$Diagnosis <- 'AD'
MSBB_meta_LFQ[ MSBB_meta_LFQ$CDR < 1 &  MSBB_meta_LFQ$Braak < 4 &  MSBB_meta_LFQ$CERAD < 2, ]$Diagnosis <- 'CT'

Meta_D$diagnosis <- as.character(Meta_D$diagnosis )
Meta_D[ Meta_D$diagnosis %in% 'control',]$diagnosis <- 'CT'
Meta_D$Diagnosis  <- Meta_D$diagnosis  

####### - BRAKK
BLSA_Fdata$BRAAK <- NA
Combat_Info$BRAAK <- Combat_Info$Braak
MSBB_meta_LFQ$BRAAK <- MSBB_meta_LFQ$Braak
Meta_D$BRAAK <- Meta_D$braaksc
MAYO_meta$BRAAK <- MAYO_meta$Braak

####### - CERAD
BLSA_Fdata$CERAD <- NA
MAYO_meta$CERAD <- NA
Combat_Info$CERAD <- Combat_Info$CERAD
# - all set - # MSBB_meta_LFQ$CERAD
Meta_D$CERAD <- Meta_D$ceradsc

####### - dcdfx
BLSA_Fdata$dcfdx_lv <- NA
MAYO_meta$dcfdx_lv <- NA
Combat_Info$dcfdx_lv <- NA
MSBB_meta_LFQ$dcfdx_lv <- NA

####### - Age of Death
BLSA_Fdata$AOD <- BLSA_Fdata$Age
#Mayo Capped at 90
MAYO_meta$AOD <- MAYO_meta$ageDeath
Combat_Info$AOD <- Combat_Info$Age
MSBB_meta_LFQ$AOD <- MSBB_meta_LFQ$ageDeath
Meta_D$AOD <- Meta_D$age_death


Comb_Vars <- c( 'SampleID', 'batch', 'study', 'tissue', 'Diagnosis', 'Sex', 'APOE4', 'AOD', 'pmi', 'dcfdx_lv', 'CERAD', 'BRAAK')

META_TOTAL <- rbind( BLSA_Fdata[, Comb_Vars ],
                     MAYO_meta[, Comb_Vars ],
                     Combat_Info[, Comb_Vars ],
                     MSBB_meta_LFQ[, Comb_Vars ],
                     Meta_D[, Comb_Vars ]
)
row.names(META_TOTAL) <- META_TOTAL$SampleID
META_TOTAL <- META_TOTAL[ colnames(Comb), ]

##################################################################################################
############## PCA Check....
pca_res <- princomp( Comp_Comb, cor = FALSE, scores = TRUE )

eigs <- pca_res$sdev^2
Proportion = eigs/sum(eigs)
#Proportion[ Proportion > 0.01 ] 

PC1_Var_exp <- Proportion[1]*100
PC2_Var_exp <- Proportion[2]*100

PCA_Plot <- META_TOTAL
PCA_Plot$PC1 <- pca_res$loadings[,1]
PCA_Plot$PC2 <- pca_res$loadings[,2]

P <- ggplot(PCA_Plot, aes( PC1, PC2) ) + geom_point( aes(colour = study, shape = tissue), cex=.8, alpha = 10/10) 
P <- P + xlab(paste0( "PC1 ", signif( PC1_Var_exp, 3), "%" )) 
P <- P + ylab(paste0( "PC2 ", signif( PC2_Var_exp, 3), "%" )) 
P <- P + ggtitle( paste0( "Proteomics PCA feature N= ", dim(Comp_Comb)[1])) + theme(plot.title = element_text(hjust = 0.5))
P

EC_foo.cluster = clusterGem(gnt = t(Comp_Comb), id = colnames(Comp_Comb), min.dim = 3, max.dim = 11)
Meta_Cust <- META_TOTAL
Meta_Cust$Clusters <- EC_foo.cluster$clusters[ row.names(Meta_Cust) ]

Clusts <- as.matrix( table( Meta_Cust$study, Meta_Cust$Clusters )) 
HM <- heatmap( Clusts )
#HM

#######################################################################################################################################
# - Meta Analysis - Model
#Comb
covar <- META_TOTAL
covar$Tissue <- paste0( covar$study, '_', covar$tissue)
#covar$Tissue <- gsub( 'ROS_DLPFC', 'ROSMAP_DLPFC', gsub( 'MAP_DLPFC', 'ROSMAP_DLPFC', covar$Tissue))

tissue.dx.summary = plyr::ddply(covar, .variables=c('Tissue', 'Diagnosis'), .fun = function(x, y){
  data.frame(peptide_id = row.names(y),
             n = dim(x)[1],
             mn = rowMeans(y[,x$SampleID], na.rm = T),
             sd = apply(y[,x$SampleID], 1, sd, na.rm = T))
}, Comb)

# Perform meta-analysis for AD-CONTROL comparison
meta.anlz.ad_cntrl = plyr::ddply(tissue.dx.summary[complete.cases(tissue.dx.summary),], .variables=c('peptide_id'), .fun = function(x){
  exp.effect = dplyr::filter(x, Diagnosis == 'AD')
  rownames(exp.effect) = exp.effect$Tissue
  cntrl.effect = dplyr::filter(x, Diagnosis == 'CT')
  rownames(cntrl.effect) = cntrl.effect$Tissue
  cntrl.effect = cntrl.effect[rownames(exp.effect), ]
  
  tmp = metacont(exp.effect$n, exp.effect$mn, exp.effect$sd, 
                 cntrl.effect$n, cntrl.effect$mn, cntrl.effect$sd,
                 studlab = exp.effect$Tissue,
                 sm = 'SMD', method.smd = 'Hedges',
                 method.tau = 'REML', control=list(maxiter=1000))
  
  return(data.frame(tmp[c('TE.fixed', 'seTE.fixed', 'lower.fixed', 'upper.fixed', 'zval.fixed', 'pval.fixed',
                          'TE.random', 'seTE.random', 'lower.random', 'upper.random', 'zval.random', 'pval.random',
                          'Q', 'tau', 'H', 'I2')]))
}, .parallel = TRUE, .paropts = list(.packages = c('meta', 'dplyr'))) %>%
  dplyr::mutate(fdr.fixed = p.adjust(pval.fixed, method = 'fdr'),
                fdr.random = p.adjust(pval.random, method = 'fdr'))

p = list()
p[[1]] = ggplot(meta.anlz.ad_cntrl, aes(x = -log10(fdr.fixed), y = -log10(fdr.random)))+geom_point()
p[[2]] = ggplot(meta.anlz.ad_cntrl, aes(y = -log10(fdr.fixed), x = TE.fixed))+geom_point()+geom_hline(yintercept = -log10(0.05), color = 'red')
p[[3]] = ggplot(meta.anlz.ad_cntrl, aes(y = -log10(fdr.random), x = TE.random))+geom_point()+geom_hline(yintercept = -log10(0.05), color = 'red')
ggpubr::ggarrange(plotlist = p, ncol = 3, nrow = 1)

# Code if random effects is significant
meta.anlz.ad_cntrl$Rand_Sig <- 'NO'
meta.anlz.ad_cntrl[ meta.anlz.ad_cntrl$fdr.random < 0.05 ,]$Rand_Sig <- 'YES'

# Code if random effects direction if it is significant
meta.anlz.ad_cntrl$Direction <- 'NONE'
meta.anlz.ad_cntrl[ meta.anlz.ad_cntrl$Rand_Sig == 'YES' & meta.anlz.ad_cntrl$TE.random < 0 ,]$Direction <- 'DOWN'
meta.anlz.ad_cntrl[ meta.anlz.ad_cntrl$Rand_Sig == 'YES' & meta.anlz.ad_cntrl$TE.random > 0 ,]$Direction <- 'UP'

table( meta.anlz.ad_cntrl$Direction, meta.anlz.ad_cntrl$Direction )


LF_DE <- cbind( do.call( rbind, strsplit( as.character(meta.anlz.ad_cntrl$peptide_id), '[|]' )), meta.anlz.ad_cntrl )
colnames(LF_DE)[c(1,2)] <- c('GName','ProtID')
LF_DE <- as.data.frame(LF_DE,  stringsAsFactors =F )
LF_DE$GName <- as.character(LF_DE$GName)

Process_Meta <- LF_DE
table(LF_DE$GName)[ table(LF_DE$GName) > 8]
LF_DE[ LF_DE$GName == 'TPM1', ] 

Process_Meta$ENSG <- NA 

table( LF_DE[ table(LF_DE$GName)[ table(LF_DE$GName) > 1], ]$Rand_Sig )
# NO YES 
# 888 147 

Dups <- names( table(LF_DE$GName)[ table(LF_DE$GName) > 1] )
Dups <- Dups[ !duplicated(Dups) ]

Col2Av <- c('TE.fixed', 'seTE.fixed', 'lower.fixed', 'upper.fixed', 'zval.fixed', 'pval.fixed', 'TE.random',
            'seTE.random', 'lower.random', 'upper.random', 'zval.random', 'pval.random', 'Q', 'tau', 'H', 'I2', 
            'fdr.fixed', 'fdr.random')

IndiciesToToss <- NULL
Process_Meta$Averaged <- 'NO'
Process_Meta$BiDirectional <- 'NO'
AbigousRows <- NULL
AbigousIDs <- NULL

### Final Name Clean
#SEPT5
#Process_Meta[ grepl('SEPT5', Process_Meta$peptide_id),]
#Process_Meta[ grepl('PALM2', Process_Meta$peptide_id),]
#Process_Meta[ grepl('ADGRB2', Process_Meta$peptide_id) | grepl('BAI2', Process_Meta$peptide_id),]
#Process_Meta[ grepl('BAI2', Process_Meta$peptide_id) | grepl('BAI2', Process_Meta$peptide_id),]

#   'ENSG00000184702', 'ENSG00000157654', 'ENSG00000121753', 'ENSG00000189060'
## SEPT5 - ENSG00000184702
#  PALM2 - ENSG00000157654
## PALM2AKAP2 - ENSG00000157654
#  ADGRB2 - ENSG00000121753 
## BAI2 - ENSG00000121753 
#  ADGRB2 - ENSG00000189060 
## BAI2 - ENSG00000189060 




for( i in 1:length(Dups)){
  Examine <- row.names(Process_Meta[ Process_Meta$GName == Dups[i], ])
  #At Least one Peptide is Significant
  if( 'YES' %in% Process_Meta[ Process_Meta$GName == Dups[i], ]$Rand_Sig ){
    TempLook <- Process_Meta[ Process_Meta$GName == Dups[i], ]
    IndiciesToToss <- c(IndiciesToToss, row.names( TempLook[ TempLook$Rand_Sig == 'NO', ] ) )
    #Are there Multiple Significant Peptides?
    if( as.numeric(table(TempLook$Rand_Sig)['YES']) > 1 ){
      yTs <- TempLook[ TempLook$Rand_Sig == 'YES', ]
      
      #Find the Ambigous Direction Protien Changes
      if( length( table(yTs$Direction)) > 1 ){
        AbigousRows <- c( AbigousRows,row.names( yTs ) )
        AbigousIDs <- c( AbigousIDs, yTs[1,]$GName )
        message( paste0('Error ', yTs[1,]$GName,' is ambiguos direction, taking Max Value') )
        yTs$BiDirectional <- 'YES'
        #yTs[ yTs$pval.random < 0.05, ]
        Keep <- row.names( yTs[ yTs$pval.random < 0.05 & 
                                  abs( yTs$TE.random ) == max(abs( yTs$TE.random ) ), ])
        
        Process_Meta[ row.names(yTs), ] <- yTs[row.names(yTs),]
        IndiciesToToss <- c(IndiciesToToss, row.names(TempLook)[ row.names(TempLook) != Keep])
      }else{
        for( COL in Col2Av){
          eval(parse(text=paste0( 'yTs$', COL, ' <- mean(yTs$', COL, ')' )))
          IndiciesToToss <- c(IndiciesToToss, row.names(yTs)[2:length(row.names(yTs))] )
          yTs$Averaged <- 'YES'
          Process_Meta[ row.names(yTs), ] <- yTs[row.names(yTs),]
        }
      }
    }else{
      IndiciesToToss <- c(IndiciesToToss, row.names(TempLook[ TempLook$Rand_Sig == 'NO', ]) )
    }
    
  }else{
    IndiciesToToss <- c(IndiciesToToss, Examine[2:length(Examine)] )
  }
}

#35 ambigous - or 9 tissueXstudy (ROSMAP) - or 7 tissueXstudy (ROS sep from MAP)
Ambigous <- c('AAK1', 'AGAP1', 'ARF3', 'ASPH', 'ATXN2', 'BMERB1', 'CADPS', 'CARNS1', 'COPS6', 'CSE1L', 'DCLK1', 'ENAH',
              'EPHB1', 'ERBIN', 'FAM219A', 'GNAL', 'GSTZ1', 'HIP1', 'LDHB', 'MACF1','MAP4', 'MYADM', 'MYL6', 'NCAM1', 'PEX11B', 
              'PNPLA6', 'PPP3R1', 'PRKACB', 'RAP1B', 'RPL31', 'SLC9A6', 'SYNGAP1', 'TLN2', 'TRAF2', 'VCL')
# 0.3172301% All Peptides
# 0.4127845% All Genes
# 0.9083831% All Sig Peptides
# 1.003728% All Sig Genes
Ambigous_7 <- c( "ASPH", "BMERB1", "ERBIN", "GNAL", "MAP4", "NCAM1", "PRKACB", "RPL31", "TLN2" )

#IndiciesToToss <- c(IndiciesToToss,
#                   row.names(Process_Meta[ Process_Meta$GName=='TUBAL3' & Process_Meta$peptide_id=='TUBAL3|A6NHL2-2' ,  ])[2])
Process_Meta <- Process_Meta[ (row.names(Process_Meta) %in% IndiciesToToss) == F , ]
table(grepl('ENSG', Process_Meta$ENSG))



Process_Meta$peptide_id <- as.character( Process_Meta$peptide_id )

TTran <- Tran[ , c('NewPeptideID', 'ENSG') ]
TTran <- TTran[ !duplicated(TTran), ]
row.names(TTran) <- TTran$NewPeptideID

Process_Meta$ENSG <- TTran[ Process_Meta$peptide_id, ]$ENSG

Process_Meta <- Process_Meta[ order(-abs(Process_Meta$TE.random)), ]


##################################################################################################
#### Rank Model:
Process_Meta <- Process_Meta[ order(abs(Process_Meta$TE.random)), ]
Process_Meta$y <- c( 1:dim(Process_Meta)[1] )/dim(Process_Meta)[1]
Process_Meta <- Process_Meta[ order(-abs(Process_Meta$TE.random)), ]
Process_Meta$Type <- 'Actual Rank'
Process_Meta$TE.random.abs<- abs( Process_Meta$TE.random )

mylogit <- glm(Process_Meta$y  ~ abs( Process_Meta$TE.random.abs), data = Process_Meta, family = "binomial")

#Make Weight Template
Work <- Process_Meta
#Work$y <- rank(abs( Work$logFC)) / max(rank(abs( Work$logFC)))
Work$y <- rank(abs( Work$TE.random.abs)) / max(rank(abs( Work$TE.random.abs)))
Work$Type <- "Actual Rank"
#plot( log(abs( Work$TE.random.abs)), Work$y )

#Fit Logistic Model:
mylogit <- glm(Work$y  ~ abs( Work$TE.random.abs), data = Work, family = "binomial")

Work$abs_TE <- abs( Work$TE.random.abs)

Work2 <- Work  
Work2$Type <- "Predicted Weight"
#Work2$y <- y <- 1/(1+exp( -( mylogit$coefficients[1]+mylogit$coefficients[2]* abs( Work$logFC) ) ))
Work2$y <- 1/(1+exp( -( mylogit$coefficients[1]+mylogit$coefficients[2]* abs( Work$TE.random.abs) ) ))
#plot( log(abs( Work2$TE.random.abs)), Work2$y )

tWork <- as.data.frame(rbind( Work ,Work2 ))
tWork$Sig <- ifelse(tWork$fdr.random < 0.05, "YES", "NO" )

#Push Non-Significant Weights to zero
tWork <- tWork[ order(-abs(tWork$y)), ]
tWork$y_Final <- tWork$y
tWork[ tWork$Sig == 'NO',]$y_Final <- 0


plot( log(abs( tWork[ tWork$Sig == 'YES' & tWork$Type == 'Predicted Weight', ]$TE.random.abs)), 
      tWork[ tWork$Sig == 'YES' & tWork$Type == 'Predicted Weight', ]$y 
)



# Save DataFrame and push to Synapse
parentId = 'syn22351719';
activityName = 'Meta analysis of differential Proteomics Data';
activityDescription = 'Fixed and random effect meta-analysis of AMP-AD Proteomics Data (4 brain regions )';
thisFileName <- 'ProteomicsProcessing_MetaModalityAnalysis.R'

# Github link
thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))

CODE <- syn_temp$store(synapseclient$Folder(name = "MetaAnalysis", parentId = parentId))
Syns_Used <- c( 'syn24216770', 'syn21266454', 'syn24828732', 'syn24828686', 'syn24828683', 'syn24828684', 'syn24828685',
                'syn21323404', 'syn21323366', 'syn23583548', 'syn3191087', 'syn23573928', 'syn18914620', 'syn18914694',
                'syn18918327', 'syn23277389', 'syn20827192', 'syn23474101', 'syn18914935', 'syn22344998', 'syn6101474',
                'syn21893059', 'syn18914939')

# Write results to files
data.table::fwrite(LF_DE, file = 'Raw_Proteomics_meta.anlz.ad_cntrl.tsv', sep = '\t', row.names = F, quote = F)
ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='Raw_Proteomics_meta.anlz.ad_cntrl.tsv', name = 'Proteomics AD-Control meta-analysis across 4 brain regions', parentId=CODE$properties$id ), used = Syns_Used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)

data.table::fwrite(Process_Meta, file = 'Processed_Proteomics_meta.anlz.ad_cntrl.tsv', sep = '\t', row.names = F, quote = F)
ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='Raw_Proteomics_meta.anlz.ad_cntrl.tsv', name = 'Proteomics AD-Control meta-analysis across 4 brain regions', parentId=CODE$properties$id ), used = Syns_Used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)

# Write results to files
data.table::fwrite(tWork, file = 'Proteomics_meta_Weights.ad_cntrl.tsv', sep = '\t', row.names = F, quote = F)
ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='Proteomics_meta_Weights.ad_cntrl.tsv', name = 'Proteomics AD-Control meta-analysis Weights', parentId=CODE$properties$id ), used = Syns_Used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)


