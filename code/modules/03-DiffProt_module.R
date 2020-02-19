Genes <- Tab$ENSG

Prot <- read.csv( file = syn_temp$get('syn21534585')$path, header = T)
Prot <- Prot[ complete.cases(Prot[,c("Log2_FC", "CI_Upr", "CI_Lwr", "PVal", "Cor_PVal")]), ]
Prot <- Prot[complete.cases(Prot),]

TargetP <- Prot[ Prot$ENSG %in% Genes,]

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

TotP %>% 
  mutate(
    AntPFC = cell_spec(AntPFC, bold = ifelse(TotPval$AntPFC < 0.05, TRUE, FALSE ),color = ifelse(TotPval$AntPFC == 0, "white", ifelse(TotPval$AntPFC < 0.05, ifelse(TotP$AntPFC < 0, "red", "green"), "grey" ))),
    DLPFC = cell_spec(DLPFC, bold = ifelse(TotPval$DLPFC < 0.05, TRUE, FALSE ), color = ifelse(TotPval$DLPFC == 0, "white", ifelse(TotPval$DLPFC < 0.05, ifelse(TotP$DLPFC < 0, "red", "green"), "grey" ))),
    MFG = cell_spec(MFG, bold = ifelse(TotPval$MFG < 0.05, TRUE, FALSE ),color = ifelse(TotPval$MFG == 0, "white", ifelse(TotPval$MFG < 0.05, ifelse(TotP$MFG < 0, "red", "green"), "grey" ))),
    TCX = cell_spec(TCX, bold = ifelse(TotPval$TCX < 0.05, TRUE, FALSE ),color = ifelse(TotPval$TCX == 0, "white", ifelse(TotPval$TCX < 0.05, ifelse(TotP$TCX < 0, "red", "green"), "grey" )))
  ) %>%
  select(GeneName, UniProt, ENSG, AntPFC, DLPFC, MFG, TCX)  %>%
  kable( escape = F) %>%
  kable_styling("striped", full_width = F)

save_kable(
  TotP %>% 
    mutate(
      AntPFC = cell_spec(AntPFC, bold = ifelse(TotPval$AntPFC < 0.05, TRUE, FALSE ),color = ifelse(TotPval$AntPFC == 0, "white", ifelse(TotPval$AntPFC < 0.05, ifelse(TotP$AntPFC < 0, "red", "green"), "grey" ))),
      DLPFC = cell_spec(DLPFC, bold = ifelse(TotPval$DLPFC < 0.05, TRUE, FALSE ), color = ifelse(TotPval$DLPFC == 0, "white", ifelse(TotPval$DLPFC < 0.05, ifelse(TotP$DLPFC < 0, "red", "green"), "grey" ))),
      MFG = cell_spec(MFG, bold = ifelse(TotPval$MFG < 0.05, TRUE, FALSE ),color = ifelse(TotPval$MFG == 0, "white", ifelse(TotPval$MFG < 0.05, ifelse(TotP$MFG < 0, "red", "green"), "grey" ))),
      TCX = cell_spec(TCX, bold = ifelse(TotPval$TCX < 0.05, TRUE, FALSE ),color = ifelse(TotPval$TCX == 0, "white", ifelse(TotPval$TCX < 0.05, ifelse(TotP$TCX < 0, "red", "green"), "grey" )))
    ) %>%
    select(GeneName, UniProt, ENSG, AntPFC, DLPFC, MFG, TCX)  %>%
    kable( escape = F) %>%
    kable_styling("striped", full_width = F),
  
  file=paste0(tabledir,'/Diff_Protiens.pdf')
)



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
    mylogit <- glm(Work$y  ~ abs( Work$logFC), data = FOO, family = "binomial")
    
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
    
    #Plot the model
    if( Tissue == 'CBE' ){
      eval(parse( text=paste0( 
        'Pr <- ggplot( data=Modelz$', Tissue, '[ Modelz$', Tissue, '$Type == \'Predicted Weight\' & (Modelz$', Tissue, '$ensembl_gene_id %in% Genes) == T , ], aes(x=log_abs_PVal, y=y, col=Sig)) +
              geom_point() + 
              geom_smooth( data=Modelz$', Tissue, '[ Modelz$', Tissue, '$Type == \'Predicted Weight\', ], aes(x=log_abs_PVal, y=y), 
                    method = "glm", 
                    method.args = list(family = "binomial"), 
                    colour="black", size=0.5, se = FALSE)  + 
              xlab("log(abs(LogFC))") + ylab("Model Weight") + ggtitle("',Tissue,'") +
              theme(plot.title = element_text(vjust = -40, hjust = .9, size = 18, face="bold") ) +
              geom_label_repel( data=Modelz$', Tissue, '[ Modelz$', Tissue, '$Type == \'Predicted Weight\' & (Modelz$', Tissue, '$ensembl_gene_id %in% Genes) == T & Modelz$', Tissue, '$Sig == \'YES\', ], aes(label = hgnc_symbol, col = Sig),
                    box.padding   = 0.35, 
                    point.padding = 0.5,
                    segment.color = \'grey50\') + 
              scale_colour_manual( values = COLS ) + 
              theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
                  axis.text.x=element_blank(), axis.ticks.x=element_blank() )'
      )))
      eval(parse( text=paste0( 'Plotz[["', Tissue, '"]] <- Pr' ) ))
    }else{
      eval(parse( text=paste0( 
        'Pr <- ggplot( data=Modelz$', Tissue, '[ Modelz$', Tissue, '$Type == \'Predicted Weight\' & (Modelz$', Tissue, '$ensembl_gene_id %in% Genes) == T , ], aes(x=log_abs_PVal, y=y, col=Sig)) +
              geom_point() + 
              geom_smooth( data=Modelz$', Tissue, '[ Modelz$', Tissue, '$Type == \'Predicted Weight\', ], aes(x=log_abs_PVal, y=y), 
                    method = "glm", 
                    method.args = list(family = "binomial"), 
                    colour="black", size=0.5, se = FALSE)  + 
              xlab("log(abs(LogFC))") + ylab("Model Weight") + ggtitle("',Tissue,'") +
              theme(plot.title = element_text(vjust = -40, hjust = .9, size = 18, face="bold") ) +
              geom_label_repel(data=Modelz$', Tissue, '[ Modelz$', Tissue, '$Type == \'Predicted Weight\' & (Modelz$', Tissue, '$ensembl_gene_id %in% Genes) == T & Modelz$', Tissue, '$Sig == \'YES\', ], aes(label = hgnc_symbol, col = Sig),
                    box.padding   = 0.35, 
                    point.padding = 0.5,
                    segment.color = \'grey50\') + 
              scale_colour_manual( values = COLS ) + 
              theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
                  axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())'
      )))
      eval(parse( text=paste0( 'Plotz[["', Tissue, '"]] <- Pr' ) ))
      
      ggsave( filename = paste0('DiffExpress_',Tissue,'.eps'),
              plot = Pr,
              device = "eps",
              path = plotdir
      )
      
    }
  }
  return( list( Plotz=Plotz, Modelz=Modelz ) )
}

#Plot/Run/Build the DE Models
Pro_ALL_Mod <- ModelConstructProt( Prot )
#Pro_ALL_Mod$Plotz$AntPFC
#Pro_ALL_Mod$Plotz$TCX
#Pro_ALL_Mod$Plotz$DLPFC
#Pro_ALL_Mod$Plotz$MFG

lay <- rbind(c(1,2),
             c(3,4))

g1 <- ggplotGrob( Pro_ALL_Mod$Plotz$AntPFC )
g2 <- ggplotGrob( Pro_ALL_Mod$Plotz$TCX )
g3 <- ggplotGrob( Pro_ALL_Mod$Plotz$DLPFC )
g4 <- ggplotGrob( Pro_ALL_Mod$Plotz$MFG )

gs = list( g1, g2, g3, g4 )

grid.arrange(grobs = gs, layout_matrix = lay)
