### Proteomics Meta Analysis Version 2
library(dplyr)
source('utilityFunctions/loadsyndata.R')
system( 'curl http://www.compgen.pitt.edu/GemTools/GemTools.R > GemTools.R')
source('GemTools.R')
synapser::synLogin()
load(synapser::synGet('syn26427476')$path)
load(synapser::synGet('syn26433992')$path)

#Make Sample Names Truly unique across the data sets:
append_deets <- function( df, string, col=NULL,cnames=FALSE){
  if(is.null(col)){
    row.names(df) <- paste0(string, '_', row.names(df))
  }else{
    if(isTRUE(cnames)){
      colnames(df) <- paste0(string, '_', colnames(df))
    }else{
      df[,col] <- paste0(string, '_', as.character(df[,col]))
    }
  }
  return(df)
}

LFQ$Banner$NuMeta <- append_deets(df=LFQ$Banner$NuMeta, string = 'lfq_banner')
LFQ$Banner$NuMeta <- append_deets(df=LFQ$Banner$NuMeta, string = 'lfq_banner', col='Sample')
LFQ$Banner$ScalWins <- append_deets(df=LFQ$Banner$ScalWins, string = 'lfq_banner', col=FALSE, cnames=TRUE)

LFQ$Blsa$NuMeta <- append_deets(df=LFQ$Blsa$NuMeta, string = 'lfq_blsa')
LFQ$Blsa$NuMeta <- append_deets(df=LFQ$Blsa$NuMeta, string = 'lfq_blsa', col='Sample')
LFQ$Blsa$ScalWins <- append_deets(df=LFQ$Blsa$ScalWins, string = 'lfq_blsa', col=FALSE, cnames=TRUE)

LFQ$Mayo$NuMeta <- append_deets(df=LFQ$Mayo$NuMeta, string = 'lfq_mayo')
LFQ$Mayo$NuMeta <- append_deets(df=LFQ$Mayo$NuMeta, string = 'lfq_mayo', col='Sample')
LFQ$Mayo$ScalWins <- append_deets(df=LFQ$Mayo$ScalWins, string = 'lfq_mayo', col=FALSE, cnames=TRUE)

LFQ$Msbb$NuMeta <- append_deets(df=LFQ$Msbb$NuMeta, string = 'lfq_msbb')
LFQ$Msbb$NuMeta <- append_deets(df=LFQ$Msbb$NuMeta, string = 'lfq_msbb', col='Sample')
LFQ$Msbb$ScalWins <- append_deets(df=LFQ$Msbb$ScalWins, string = 'lfq_msbb', col=FALSE, cnames=TRUE)

TMT$TMT_banner_BA9$NuMeta$Sample <- row.names(TMT$TMT_banner_BA9$NuMeta)
TMT$TMT_banner_BA9$NuMeta <- append_deets(df=TMT$TMT_banner_BA9$NuMeta, string = 'tmt_banner_ba9')
TMT$TMT_banner_BA9$NuMeta <- append_deets(df=TMT$TMT_banner_BA9$NuMeta, string = 'tmt_banner_ba9', col='Sample')
TMT$TMT_banner_BA9$ScalWins <- append_deets(df=TMT$TMT_banner_BA9$ScalWins, string = 'tmt_banner_ba9', col=FALSE, cnames=TRUE)

TMT$TMT_rosmap_BA9$NuMeta$Sample <- row.names(TMT$TMT_rosmap_BA9$NuMeta)
TMT$TMT_rosmap_BA9$NuMeta <- append_deets(df=TMT$TMT_rosmap_BA9$NuMeta, string = 'tmt_rosmap_ba9')
TMT$TMT_rosmap_BA9$NuMeta <- append_deets(df=TMT$TMT_rosmap_BA9$NuMeta, string = 'tmt_rosmap_ba9', col='Sample')
TMT$TMT_rosmap_BA9$ScalWins <- append_deets(df=TMT$TMT_rosmap_BA9$ScalWins, string = 'tmt_rosmap_ba9', col=FALSE, cnames=TRUE)

TMT$TMT_emory_BA24$NuMeta <- append_deets(df=TMT$TMT_emory_BA24$NuMeta, string = 'tmt_emory_ba24')
TMT$TMT_emory_BA24$NuMeta <- append_deets(df=TMT$TMT_emory_BA24$NuMeta, string = 'tmt_emory_ba24', col='Sample')
TMT$TMT_emory_BA24$ScalWins <- append_deets(df=TMT$TMT_emory_BA24$ScalWins, string = 'tmt_emory_ba24', col=FALSE, cnames=TRUE)

TMT$TMT_emory_BA9$NuMeta <- append_deets(df=TMT$TMT_emory_BA9$NuMeta, string = 'tmt_emory_ba9')
TMT$TMT_emory_BA9$NuMeta <- append_deets(df=TMT$TMT_emory_BA9$NuMeta, string = 'tmt_emory_ba9', col='Sample')
TMT$TMT_emory_BA9$ScalWins <- append_deets(df=TMT$TMT_emory_BA9$ScalWins, string = 'tmt_emory_ba9', col=FALSE, cnames=TRUE)

TMT$TMT_msbb_BA36$NuMeta <- append_deets(df=TMT$TMT_msbb_BA36$NuMeta, string = 'tmt_msbb_ba36')
TMT$TMT_msbb_BA36$NuMeta <- append_deets(df=TMT$TMT_msbb_BA36$NuMeta, string = 'tmt_msbb_ba36', col='Sample')
TMT$TMT_msbb_BA36$ScalWins <- append_deets(df=TMT$TMT_msbb_BA36$ScalWins, string = 'tmt_msbb_ba36', col=FALSE, cnames=TRUE)

TMT$TMT_rosmap_BA6$NuMeta <- append_deets(df=TMT$TMT_rosmap_BA6$NuMeta, string = 'tmt_rosmap_ba6')
TMT$TMT_rosmap_BA6$NuMeta <- append_deets(df=TMT$TMT_rosmap_BA6$NuMeta, string = 'tmt_rosmap_ba6', col='Sample')
TMT$TMT_rosmap_BA6$ScalWins <- append_deets(df=TMT$TMT_rosmap_BA6$ScalWins, string = 'tmt_rosmap_ba6', col=FALSE, cnames=TRUE)

TMT$TMT_rosmap_BA37$NuMeta <- append_deets(df=TMT$TMT_rosmap_BA37$NuMeta, string = 'tmt_rosmap_ba37')
TMT$TMT_rosmap_BA37$NuMeta <- append_deets(df=TMT$TMT_rosmap_BA37$NuMeta, string = 'tmt_rosmap_ba37', col='Sample')
TMT$TMT_rosmap_BA37$ScalWins <- append_deets(df=TMT$TMT_rosmap_BA37$ScalWins, string = 'tmt_rosmap_ba37', col=FALSE, cnames=TRUE)

all_features <-  c(
  row.names(LFQ$Banner$ScalWins),
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
all_features <- all_features[!duplicated(all_features)]

expression <- dplyr::full_join(as.data.frame(LFQ$Banner$ScalWins) %>% tibble::rownames_to_column('rn'),
                               as.data.frame(LFQ$Blsa$ScalWins) %>% tibble::rownames_to_column('rn'),
                               by="rn") %>%
  dplyr::full_join(as.data.frame(LFQ$Mayo$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(LFQ$Msbb$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_banner_BA9$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_rosmap_BA9$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_emory_BA24$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_emory_BA9$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_msbb_BA36$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_rosmap_BA6$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_rosmap_BA37$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  tibble::column_to_rownames('rn')

keep_cols <- c( 'Sample', 'Age',  'Gender', 'PMI', 'AD', 'Control', 'Braak', 'CERAD', 'SampleForCombat')
LFQ$Mayo$NuMeta$CERAD <- NA

TMT$TMT_rosmap_BA9$NuMeta$Gender <- as.factor(TMT$TMT_rosmap_BA9$NuMeta$Gender)
TMT$TMT_emory_BA24$NuMeta$Age <- as.numeric(TMT$TMT_emory_BA24$NuMeta$Age)
TMT$TMT_emory_BA24$NuMeta$Gender <- as.factor(TMT$TMT_emory_BA24$NuMeta$Gender)
TMT$TMT_emory_BA24$NuMeta$PMI <- as.numeric(TMT$TMT_emory_BA24$NuMeta$PMI)
TMT$TMT_emory_BA24$NuMeta$Braak <- as.numeric(TMT$TMT_emory_BA24$NuMeta$Braak)

TMT$TMT_emory_BA9$NuMeta$Age <- as.numeric(TMT$TMT_emory_BA9$NuMeta$Age)
TMT$TMT_emory_BA9$NuMeta$Gender <- as.factor(TMT$TMT_emory_BA9$NuMeta$Gender)
TMT$TMT_emory_BA9$NuMeta$PMI <- as.numeric(TMT$TMT_emory_BA9$NuMeta$PMI)
TMT$TMT_emory_BA9$NuMeta$Braak <- as.numeric(TMT$TMT_emory_BA9$NuMeta$Braak)

TMT$TMT_msbb_BA36$NuMeta$Gender <- as.factor(TMT$TMT_msbb_BA36$NuMeta$Gender)
TMT$TMT_rosmap_BA6$NuMeta$Gender <- as.factor(TMT$TMT_rosmap_BA6$NuMeta$Gender)
TMT$TMT_rosmap_BA37$NuMeta$Gender <- as.factor(TMT$TMT_rosmap_BA37$NuMeta$Gender)

LFQ$Banner$NuMeta$data <- 'LFQ'
LFQ$Blsa$NuMeta$data <- 'LFQ'
LFQ$Mayo$NuMeta$data <- 'LFQ'
LFQ$Msbb$NuMeta$data <- 'LFQ'
TMT$TMT_banner_BA9$NuMeta$data <- 'TMT'
TMT$TMT_rosmap_BA9$NuMeta$data <- 'TMT'
TMT$TMT_emory_BA24$NuMeta$data <- 'TMT'
TMT$TMT_emory_BA9$NuMeta$data <- 'TMT'
TMT$TMT_msbb_BA36$NuMeta$data <- 'TMT'
TMT$TMT_rosmap_BA6$NuMeta$data <- 'TMT'
TMT$TMT_rosmap_BA37$NuMeta$data <- 'TMT'

LFQ$Banner$NuMeta$region <- 'TCX'
LFQ$Blsa$NuMeta$region <- 'MFG'
LFQ$Mayo$NuMeta$region <- 'TCX'
LFQ$Msbb$NuMeta$region <- 'AntPFC'
TMT$TMT_banner_BA9$NuMeta$region <- 'DLPFC'
TMT$TMT_rosmap_BA9$NuMeta$region <- 'DLPFC'
TMT$TMT_emory_BA24$NuMeta$region <- 'BA24'
TMT$TMT_emory_BA9$NuMeta$region <- 'DLPFC'
TMT$TMT_msbb_BA36$NuMeta$region <- 'PHG'
TMT$TMT_rosmap_BA6$NuMeta$region <- 'BA6'
TMT$TMT_rosmap_BA37$NuMeta$region <- 'BA37'

LFQ$Banner$NuMeta$site <- 'banner'
LFQ$Blsa$NuMeta$site <- 'blsa'
LFQ$Mayo$NuMeta$site <- 'mayo'
LFQ$Msbb$NuMeta$site <- 'msbb'
TMT$TMT_banner_BA9$NuMeta$site <- 'banner'
TMT$TMT_rosmap_BA9$NuMeta$site <- 'rosmap'
TMT$TMT_emory_BA24$NuMeta$site <- 'emory'
TMT$TMT_emory_BA9$NuMeta$site <- 'emory'
TMT$TMT_msbb_BA36$NuMeta$site <- 'msbb'
TMT$TMT_rosmap_BA6$NuMeta$site <- 'rosmap'
TMT$TMT_rosmap_BA37$NuMeta$site <- 'rosmap'

keep_cols <- c(keep_cols,'data','region', 'site')
meta <- dplyr::bind_rows(
  LFQ$Banner$NuMeta[,keep_cols],
  LFQ$Blsa$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(LFQ$Mayo$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(LFQ$Msbb$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(TMT$TMT_banner_BA9$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(TMT$TMT_rosmap_BA9$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(TMT$TMT_emory_BA24$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(TMT$TMT_emory_BA9$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(TMT$TMT_msbb_BA36$NuMeta[,keep_cols]) %>%
  dplyr::bind_rows(TMT$TMT_rosmap_BA6$NuMeta[,keep_cols])%>%
  dplyr::bind_rows(TMT$TMT_rosmap_BA37$NuMeta[,keep_cols])

expression <- dplyr::full_join(as.data.frame(LFQ$Banner$ScalWins) %>% tibble::rownames_to_column('rn'),
                               as.data.frame(LFQ$Blsa$ScalWins) %>% tibble::rownames_to_column('rn'),
                               by="rn") %>%
  dplyr::full_join(as.data.frame(LFQ$Mayo$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(LFQ$Msbb$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_banner_BA9$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_rosmap_BA9$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_emory_BA24$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_emory_BA9$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_msbb_BA36$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_rosmap_BA6$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  dplyr::full_join(as.data.frame(TMT$TMT_rosmap_BA37$ScalWins) %>% tibble::rownames_to_column('rn'), by="rn") %>%
  tibble::column_to_rownames('rn')


meta <- meta[colnames(expression),]

# PCA Comps
############## PCA Check....
pca_res <- princomp( expression[complete.cases(expression),], cor = F, scores = TRUE )
eigs <- pca_res$sdev^2
Proportion = eigs/sum(eigs)

PC1_Var_exp <- Proportion[1]*100
PC2_Var_exp <- Proportion[2]*100

PCA_Plot <- meta
PCA_Plot$PC1 <- pca_res$loadings[,1]
PCA_Plot$PC2 <- pca_res$loadings[,2]

P <- ggplot(PCA_Plot, aes( PC1, PC2) ) + geom_point( aes(colour = data, shape = region), cex=.8, alpha = 10/10) 
P <- P + xlab(paste0( "PC1 ", signif( PC1_Var_exp, 3), "%" )) 
P <- P + ylab(paste0( "PC2 ", signif( PC2_Var_exp, 3), "%" )) 
P <- P + ggtitle( paste0( "Proteomics PCA fully covered features N= ", dim(expression[complete.cases(expression),])[1])) + theme(plot.title = element_text(hjust = 0.5))
P

EC_foo.cluster = clusterGem(gnt = t(expression[complete.cases(expression),]), id = colnames(expression), min.dim = 3, max.dim = 11)
Meta_Cust <- meta
Meta_Cust$Clusters <- EC_foo.cluster$clusters[ row.names(Meta_Cust) ]

Clusts <- as.matrix( table( Meta_Cust$region, Meta_Cust$Clusters )) 
HM <- heatmap( Clusts )
Clusts <- as.matrix( table( Meta_Cust$data, Meta_Cust$Clusters )) 
HM <- heatmap( Clusts )

#######################################################################################################################################
# - Meta Analysis - Model
#Comb
covar <- meta
covar$cohort <- paste0( covar$data, '_', covar$site, '_', covar$region)
#covar$Tissue <- gsub( 'ROS_DLPFC', 'ROSMAP_DLPFC', gsub( 'MAP_DLPFC', 'ROSMAP_DLPFC', covar$Tissue))

table(covar$SampleForCombat)
covar <- covar[ covar$SampleForCombat %in% c('AD','CT'), ]
covar$Diagnosis <- covar$SampleForCombat

Comb <- expression[,row.names(covar)]

tissue.dx.summary = plyr::ddply(covar, .variables=c('cohort', 'Diagnosis'), .fun = function(x, y){
  data.frame(peptide_id = row.names(y),
             n = dim(x)[1],
             mn = rowMeans(y[,x$Sample], na.rm = T),
             sd = apply(y[,x$Sample], 1, sd, na.rm = T))
}, as.matrix(Comb))

# Perform meta-analysis for AD-CONTROL comparison
meta.anlz.ad_cntrl = plyr::ddply(tissue.dx.summary[complete.cases(tissue.dx.summary),], .variables=c('peptide_id'), .fun = function(x){
  exp.effect = dplyr::filter(x, Diagnosis == 'AD')
  rownames(exp.effect) = exp.effect$cohort
  cntrl.effect = dplyr::filter(x, Diagnosis == 'CT')
  rownames(cntrl.effect) = cntrl.effect$cohort
  cntrl.effect = cntrl.effect[rownames(exp.effect), ]
  
  tmp = metacont(exp.effect$n, exp.effect$mn, exp.effect$sd, 
                 cntrl.effect$n, cntrl.effect$mn, cntrl.effect$sd,
                 studlab = exp.effect$cohort,
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

# Annotate with ENSG and Process Duplicates:

load(synapser::synGet('syn26433992')$path)
trans <- feature_translation$Proteomics
row.names(trans) <- trans$ProteomicsID

meta.anlz.ad_cntrl$UniprotID <- trans[meta.anlz.ad_cntrl$peptide_id,]$UniprotID
meta.anlz.ad_cntrl$SYMBOL <- trans[meta.anlz.ad_cntrl$peptide_id,]$SYMBOL
meta.anlz.ad_cntrl$EnsemblGene <- trans[meta.anlz.ad_cntrl$peptide_id,]$EnsemblGene
##################################################################################################


#### Rank Model:
Process_Meta <- meta.anlz.ad_cntrl
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

tWork$plotname <- NA
label_name <- c( "MTCH1", "MTOR", "APPL1", "BAX", "STMN2", "SLC9A9", "TOMM40L", "GABRA2",  "MTCH2"  )
tWork[ tWork$SYMBOL %in% label_name, ]$plotname <- tWork[ tWork$SYMBOL %in% label_name, ]$SYMBOL

A <- tWork[ tWork$Type == 'Actual Rank', ]
row.names(A) <- A$ENSG
colnames(A)[colnames(A)=='y_Final'] <- 'y_Final_Actual'
colnames(A)[colnames(A)=='y'] <- 'y_Actual'

B <- tWork[ tWork$Type == 'Predicted Weight', ]
row.names(B) <- B$ENSG
colnames(B)[colnames(B)=='y_Final'] <- 'y_Final_Predicted'
colnames(B)[colnames(B)=='y'] <- 'y_Predicted'

plotframe <- as.data.frame(cbind(A[,c('peptide_id', 'EnsemblGene','SYMBOL', 'TE.random', 'seTE.random', 'fdr.random', 'abs_TE', 'y_Actual','y_Final_Actual', 'plotname')], B[row.names(A), c('y_Predicted','y_Final_Predicted')]))

plotframe <- plotframe[ plotframe$fdr.random < 0.05, ]
plotframe$logTE <- log(abs(plotframe$abs_TE))

library(ggplot2)
library(ggrepel)
ggplot2::ggplot(data=plotframe) + 
  ggplot2::geom_point( ggplot2::aes(x=logTE, y=y_Predicted)) +
  ggrepel::geom_label_repel(data = plotframe, 
                            ggplot2::aes(x=logTE, y=y_Predicted, label = plotname, size = NULL, color = NULL),
                            nudge_x = .45,
                            nudge_y = -.05,
                            segment.size  = 0.2,
                            segment.color = "grey50",
                            direction     = "x"
  ) +
  ggplot2::xlab("Log( absolute(Effect Size))") + 
  ggplot2::ylab("Assigned weight") +
  ggplot2::labs(title="Weights of Differentially Expressed Proteins") + 
  ggplot2::theme(plot.title =  ggplot2::element_text(hjust = 0.5))




plot( log(abs( tWork[ tWork$Sig == 'YES' & tWork$Type == 'Predicted Weight', ]$TE.random.abs)), 
      tWork[ tWork$Sig == 'YES' & tWork$Type == 'Predicted Weight', ]$y, pch =16, cex = .3, las =1,
      xlab = "Log2( absolute(Effect Size))", ylab = "Assigned weight", main = "Weights of Differentially Expressed Proteins"
)

###################
# Save DataFrame and push to Synapse
parentId = 'syn22351719';
activityName = 'Meta analysis of differential Proteomics Data';
activityDescription = 'Fixed and random effect meta-analysis of AMP-AD Proteomics Data (4 brain regions )';
thisFileName <- 'Proteomics_Metaanalysis_V2.R'

# Github link
thisRepo <- githubr::getRepo(repository = "Sage-Bionetworks/AD_TargetRank", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))

CODE <- synapser::synGet('syn22414618')
Syns_Used <- c( 'syn24216770', 'syn21266454', 'syn24828732', 'syn24828686', 'syn24828683', 'syn24828684', 'syn24828685',
                'syn21323404', 'syn21323366', 'syn23583548', 'syn3191087', 'syn23573928', 'syn18914620', 'syn18914694',
                'syn18918327', 'syn23277389', 'syn20827192', 'syn23474101', 'syn18914935', 'syn22344998', 'syn6101474',
                'syn21893059', 'syn18914939', 'syn26427476', 'syn26433992'
)

# Write results to files
data.table::fwrite(
  meta.anlz.ad_cntrl, 
  file = 'Raw_Proteomics_meta.anlz.ad_cntrl.tsv', 
  sep = '\t', row.names = F, quote = F
)
ENRICH_OBJ <-  synapser::synStore(
  synapser::File(
    path='Raw_Proteomics_meta.anlz.ad_cntrl.tsv', 
    name = 'Proteomics AD-Control meta-analysis across 4 brain regions', 
    parentId=CODE$properties$id ), 
  used = Syns_Used, activityName = activityName, 
  executed = thisFile, activityDescription = activityDescription)

data.table::fwrite(
  Process_Meta, 
  file = 'Processed_Proteomics_meta.anlz.ad_cntrl.tsv', 
  sep = '\t', row.names = F, quote = F
)
ENRICH_OBJ <-  synapser::synStore(
  synapser::File(
    path='Raw_Proteomics_meta.anlz.ad_cntrl.tsv', 
    name = 'Proteomics AD-Control meta-analysis across 4 brain regions', 
    parentId=CODE$properties$id ), used = Syns_Used, 
  activityName = activityName, executed = thisFile, 
  activityDescription = activityDescription
)

# Write results to files
data.table::fwrite(
  tWork, 
  file = 'Proteomics_meta_Weights.ad_cntrl.tsv', 
  sep = '\t', row.names = F, quote = F
)
ENRICH_OBJ <-  synapser::synStore(
  synapser::File(
    path='Proteomics_meta_Weights.ad_cntrl.tsv', 
    name = 'Proteomics AD-Control meta-analysis Weights', 
    parentId=CODE$properties$id ), 
  used = Syns_Used, 
  activityName = activityName, 
  executed = thisFile, 
  activityDescription = activityDescription)




###################

keeps <- c( 'peptide_id', 'UniprotID', 'SYMBOL', 'EnsemblGene','TE.random',
            'seTE.random', 'lower.random', 'upper.random', 'zval.random', 
            'pval.random', 'Q','tau', 'H', 'I2', 'fdr.random', 'Rand_Sig',
            'Direction')

meta_analysis_values <- meta.anlz.ad_cntrl[, keeps]


