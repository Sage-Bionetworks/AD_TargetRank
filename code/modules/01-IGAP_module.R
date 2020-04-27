#Load IGAP files
#IGAP <- read.table( synfile = "/Users/jgockley/Desktop/Projects/AMP-AD_Targets/Candidate_Target_Profiler/RefrenceData/IGAP_stage_1_2_combined.txt", sep = "\t", header = T)
#IGAP1 <- read.table( file = "/Users/jgockley/Desktop/Projects/AMP-AD_Targets/Candidate_Target_Profiler/RefrenceData/IGAP_stage_1.txt", sep = "\t", header = T)
syn_temp$login()
IGAP <- read.table( syn_temp$get('syn21534587')$path, sep = "\t", header = T)
IGAP1 <- read.table( syn_temp$get('syn21534586')$path, sep = "\t", header = T)

Tab$LowestPValCombined <- 1
Tab$TotalSNPsCombined <- 0

Tab$LowestPValPhase1 <- 1
Tab$TotalSNPsPhase1 <- 0

#New Dataframe for the IGAP results that you can track the pval distributions
TotalIGAP <- NULL

#Look at the dataframe 
MedianScoreIGAP <- NULL

#Look at the dataframe 
BestScoreIGAP <- NULL

#save.image("TempA.RData")

for( i in 1:dim(Tab)[1]){
  
  if( isTRUE(as.character(Tab$chromosome_name)[i] %in% as.character(c(1:22))) ){
    Temp <- IGAP[ (as.numeric(as.character(IGAP$Position)) > (as.numeric(as.character(Tab$start_position))[i]-1e6)) &
                    (as.numeric(as.character(IGAP$Position)) < (as.numeric(as.character(Tab$start_position))[i]+1e6)) &
                    as.numeric(as.character(IGAP$Chromosome)) == as.numeric(as.character(Tab$chromosome_name)[i])  , ]
    if( dim(Temp)[1] == 0){
      Temp[1,] <- 0
      Temp$GeneName <- Tab$hgnc_symbol[i]
    }else{
      Temp$GeneName <- Tab$hgnc_symbol[i]
    }
    TotalIGAP <- as.data.frame(rbind(TotalIGAP,Temp))
    
    TempV1 <- IGAP1[ (as.numeric(as.character(IGAP1$Position)) > (as.numeric(as.character(Tab$start_position))[i]-1e6)) &
                       (as.numeric(as.character(IGAP1$Position)) < (as.numeric(as.character(Tab$start_position))[i]+1e6)) &
                       as.numeric(as.character(IGAP1$Chromosome)) == as.numeric(as.character(Tab$chromosome_name)[i])  , ] 
    if( dim(TempV1)[1]>0 ){
      TempV1$GeneName <- Tab$hgnc_symbol[i]
      #####
      TotalIGAP <- as.data.frame(rbind(TotalIGAP,TempV1))
    }else{}
    
    if( median(c(1:length(TempV1$Pvalue))) == floor(median(c(1:length(TempV1$Pvalue)))) & median(c(1:length(TempV1$Pvalue))) == ceiling(median(c(1:length(TempV1$Pvalue)))) ) {
      MedianScoreIGAP <- as.data.frame(rbind(MedianScoreIGAP,TempV1[ median(TempV1$Pvalue) == as.numeric(TempV1$Pvalue), ]))
    }else{
      
      val <- sort(TempV1$Pvalue)[ floor(median(c(1:length(TempV1$Pvalue)))) ]
      MedianScoreIGAP <- as.data.frame(rbind(MedianScoreIGAP,TempV1[ val == as.numeric(TempV1$Pvalue), ]))
    }
    
    if( dim(Temp)[1] == 0 ){
      
    }else{
      Tab$LowestPValCombined[i] <- min( Temp$Pvalue )
      Tab$TotalSNPsCombined[i] <- paste0("N=", dim(Temp)[1] )
    }
    if( dim(TempV1)[1] == 0 ){
    }else{
      Tab$LowestPValPhase1[i] <- min( TempV1$Pvalue )
      Tab$TotalSNPsPhase1[i] <- paste0("N=", dim(TempV1)[1] )
    }
  }else{
    #Temp <- as.data.frame(matrix(0, 1, 9))
    #colnames(Temp) <- c( 'Chromosome', 'Position', 'MarkerName', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue', 'GeneName' ) 
    #Temp$GeneName <- Tab[ i, ]$GeneName
    #Temp$Pvalue <- 1
    #BestScoreIGAP <- as.data.frame(rbind(BestScoreIGAP,Temp[ 1, ]))
    
  }
  if( dim(Temp)[1] == 0){
    Temp <- as.data.frame(matrix(0, 1, 10))
    colnames(Temp) <- c( 'Chromosome', 'Position', 'MarkerName', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue', 'Rank', 'GeneName' ) 
    Temp$GeneName <- Tab[ i, ]$hgnc_symbol
    Temp$Pvalue <- 1
    BestScoreIGAP <- as.data.frame(rbind(BestScoreIGAP,Temp[ 1, ]))
  }else{
    BestScoreIGAP <- as.data.frame(rbind(BestScoreIGAP,Temp[ min(Temp$Pvalue) == as.numeric(Temp$Pvalue), ]))
  }
}
colnames(Tab) <- c("ENSG", "GeneName", "Chr", "Start", "End", "Coordhg19", "Intervalhg19", "LowestPValCombined", "TotalSNPsCombined", "LowestPValPhaseI", "TotalSNPsPhaseI" )
row.names(Tab) <- Tab$ENSG

TotalIGAP$LogP <- -log(TotalIGAP$Pvalue)

#Replace Infinite values with max + one
if( Inf %in% TotalIGAP$LogP ){
  TotalIGAP[ TotalIGAP$LogP == Inf, ]$LogP <- max(TotalIGAP$LogP[TotalIGAP$LogP!=Inf]+1)
}else{}
p <- ggplot(TotalIGAP, aes(x=GeneName, y=LogP)) + 
  geom_violin( ) + geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_jitter(shape=16, size=.1, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle = 45) ) + ylab('-Log(P-Value)')  
#p 

setEPS()
postscript(file = paste0(plotdir,'/TotalIGAP.eps') )
ggplot(TotalIGAP, aes(x=GeneName, y=LogP)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  geom_jitter(shape=16, size=.1, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle = 45) ) + ylab('-Log(P-Value)') 
dev.off()


IGAP$Rank<-rank(-IGAP$Pvalue)
FOO<-IGAP

#Replace Infinite values with max + one
if( Inf %in% FOO$Pvalue ){
  FOO[ FOO$Pvalue == Inf,]$Pvalue <- max(FOO$Pvalue[FOO$Pvalue!=Inf]+1)
}else{
  
}

FOO$Pvalue <- -log(FOO$Pvalue)
#Replace Infinite values with max + one and zero vals w/ min/100
FOO[ FOO$Pvalue == Inf,]$Pvalue <- max(FOO$Pvalue[FOO$Pvalue!=Inf]+1)
#FOO[ FOO$Pvalue == 0,]$Pvalue <- min(FOO$Pvalue[FOO$Pvalue>0])/100

FOO$RankModel <- (FOO$Rank/max(FOO$Rank))

Pr <- ggplot( FOO, aes(x=Pvalue, y=RankModel)) +
  geom_point() + xlab("-log(P-Value") + geom_smooth(method = "glm", method.args = list(family = quasibinomial(link = 'logit')), se = FALSE)
#Pr

setEPS()
postscript(file = paste0(plotdir,'/IGAP_WeightModel.eps') )
ggplot( FOO, aes(x=Pvalue, y=RankModel)) +
  geom_point() + xlab("-log(P-Value") + geom_smooth(method = "glm", method.args = list(family = quasibinomial(link = 'logit')), se = FALSE)
dev.off()

#Fit Logistic Model:
FOO$Rank<-rank(FOO$Pvalue)
FOO$RankModel <- (FOO$Rank/dim(FOO)[1])
mylogit <- glm(RankModel ~ Pvalue, data = FOO, family = "binomial")
family = quasibinomial(link = 'logit')
#FOO$RankModel

LogisticScorer <- function( X ){
  y <- 1/(1+exp( -( as.numeric(mylogit$coefficients[1]) + as.numeric(mylogit$coefficients[2]) * X) ))
  return(y)
}

#LogisticScorer <- function( X ){
#  y <- 1/(1+exp( -( -3.305501+0.420491*X) ))
#  return(y)
#}

FOO$Weights <- LogisticScorer(FOO$Pvalue)

Temp <- FOO[ ,c(1:10)]
Temp$Type <- "Actual Rank"
colnames(Temp)[10] <- "y"
temp <- FOO[ ,c(1:9,11)] 
temp$Type <- "Predicted Weight"
colnames(temp)[10] <- "y"

Total <- rbind(Temp,temp)
Pt <- ggplot( Total, aes(x=Pvalue, y=y, col=Type)) + scale_color_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_point() + xlab("-log(P-Value)") + geom_smooth(method = "glm", method.args = list(family = quasibinomial(link = 'logit')), colour="black", size=0.5, se = FALSE)  
#Pt

setEPS()
postscript(file = paste0(plotdir,'/Fitted_IGAP_WeightModel.eps') )
ggplot( Total, aes(x=Pvalue, y=y, col=Type)) + scale_color_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_point() + xlab("-log(P-Value)") + geom_smooth(method = "glm", method.args = list(family = quasibinomial(link = 'logit')), colour="black", size=0.5, se = FALSE)  
dev.off()

#Remove the pval ties
BestScoreIGAP <- BestScoreIGAP[ !duplicated(BestScoreIGAP$GeneName), ]
BestScoreIGAP$Wts <- LogisticScorer(-log(BestScoreIGAP$Pvalue))
BestScoreIGAP$LogP <- -log(BestScoreIGAP$Pvalue)

Py <- ggplot( BestScoreIGAP, aes(x=LogP, y=Wts, col=GeneName)) +
  geom_point() + xlab("-log(P-Value)") + ylab("IGAP Weight") + 
  geom_smooth( data=Total, aes(x=Pvalue, y=y, col=Type), method = "glm", method.args = list(family = quasibinomial(link = 'logit')), colour="black", size=0.5, se = FALSE) + 
  geom_label_repel(aes(label = GeneName), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + theme(legend.position = "none")

#Py
setEPS()
postscript(file = paste0(plotdir,'/Best_IGAP_Weight_Score.eps') )
ggplot( BestScoreIGAP, aes(x=LogP, y=Wts, col=GeneName)) +
  geom_point() + xlab("-log(P-Value)") + ylab("IGAP Weight") + 
  geom_smooth( data=Total, aes(x=Pvalue, y=y, col=Type), method = "glm", method.args = list(family = quasibinomial(link = 'logit')), colour="black", size=0.5, se = FALSE) + 
  geom_label_repel(aes(label = GeneName), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')  
dev.off()

lay <- rbind(c(1,2),c(3,4))

g1 <- ggplotGrob( p)
g2 <- ggplotGrob( Pr)
g3 <- ggplotGrob( Pt)
g4 <- ggplotGrob( Py )

gs = list( g1, g2, g3, g4 )

plot( grid.arrange(grobs = gs, layout_matrix = lay, top = "IGAP Weight Models Applied to Lowest SNP Within 1MB of Combined Phase I + II)") )

##Plot Data Table
Tab %>% 
  mutate(
    LowestPValCombined = cell_spec(LowestPValCombined, bold = ifelse(LowestPValCombined < 0.05, TRUE, FALSE ), color = ifelse(LowestPValCombined < 0.05, "red", ifelse(LowestPValCombined == 1, "white", "black") )),
    LowestPValPhaseI = cell_spec(LowestPValPhaseI, bold = ifelse(LowestPValPhaseI < 0.05, TRUE, FALSE ), color = ifelse(LowestPValPhaseI < 0.05, "red", ifelse(LowestPValPhaseI == 1, "white", "black") ))
  ) %>%
  select( ENSG, GeneName, Start, End, Chr, Coordhg19, Intervalhg19, LowestPValCombined, TotalSNPsCombined, LowestPValPhaseI, TotalSNPsPhaseI)  %>%
  kable( escape = F,  format = "html") %>%
  kable_styling("striped", full_width = F)


save_kable(
  Tab %>% 
    mutate(
      LowestPValCombined = cell_spec(LowestPValCombined, bold = ifelse(LowestPValCombined < 0.05, TRUE, FALSE ), color = ifelse(LowestPValCombined < 0.05, "red", ifelse(LowestPValCombined == 1, "white", "black") )),
      LowestPValPhaseI = cell_spec(LowestPValPhaseI, bold = ifelse(LowestPValPhaseI < 0.05, TRUE, FALSE ), color = ifelse(LowestPValPhaseI < 0.05, "red", ifelse(LowestPValPhaseI == 1, "white", "black") ))
    ) %>%
    select( ENSG, GeneName, Start, End, Chr, Coordhg19, Intervalhg19, LowestPValCombined, TotalSNPsCombined, LowestPValPhaseI, TotalSNPsPhaseI)  %>%
    kable( escape = F,  format = "html") %>%
    kable_styling("striped", full_width = F),
  
  file=paste0(tabledir,'/IGAPStats.pdf')
)

