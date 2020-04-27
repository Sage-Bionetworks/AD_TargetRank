DE <- read.table( syn_temp$get('syn21534583')$path, header=T, sep ="\t")
Genes <- Tab$ENSG
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
    mylogit <- glm(Work$y  ~ abs( Work$logFC), data = FOO, family = "binomial")
    
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
    if( Tissue == 'CBE' ){
      eval(parse( text=paste0( 
        'Pr <- ggplot( data=Models$', Tissue, '[ Models$', Tissue, '$Type == \'Predicted Weight\' & (Models$', Tissue, '$ensembl_gene_id %in% Genes) == T , ], aes(x=log_abs_PVal, y=y, col=Sig)) +
              geom_point() + 
              geom_smooth( data=Models$', Tissue, '[ Models$', Tissue, '$Type == \'Predicted Weight\', ], aes(x=log_abs_PVal, y=y), 
                    method = "glm", 
                    method.args = list(family = "binomial"), 
                    colour="black", size=0.5, se = FALSE)  + 
              xlab("log(abs(LogFC))") + ylab("Model Weight") + ggtitle("',Tissue,'") +
              theme(plot.title = element_text(vjust = -60, hjust = .9, size = 18, face="bold") ) +
              geom_label_repel(aes(label = hgnc_symbol, col = Sig),
                    box.padding   = 0.35, 
                    point.padding = 0.5,
                    segment.color = \'grey50\') + 
              scale_colour_manual( values = COLS ) + 
              theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
                  axis.text.x=element_blank(), axis.ticks.x=element_blank() )'
      )))
      eval(parse( text=paste0( 'Plots[["', Tissue, '"]] <- Pr' ) ))
    }else{
      eval(parse( text=paste0( 
        'Pr <- ggplot( data=Models$', Tissue, '[ Models$', Tissue, '$Type == \'Predicted Weight\' & (Models$', Tissue, '$ensembl_gene_id %in% Genes) == T , ], aes(x=log_abs_PVal, y=y, col=Sig)) +
              geom_point() + 
              geom_smooth( data=Models$', Tissue, '[ Models$', Tissue, '$Type == \'Predicted Weight\', ], aes(x=log_abs_PVal, y=y), 
                    method = "glm", 
                    method.args = list(family = "binomial"), 
                    colour="black", size=0.5, se = FALSE)  + 
              xlab("log(abs(LogFC))") + ylab("Model Weight") + ggtitle("',Tissue,'") +
              theme(plot.title = element_text(vjust = -60, hjust = .9, size = 18, face="bold") ) +
              geom_label_repel(aes(label = hgnc_symbol, col = Sig),
                    box.padding   = 0.35, 
                    point.padding = 0.5,
                    segment.color = \'grey50\') + 
              scale_colour_manual( values = COLS ) + 
              theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
                  axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())'
      )))
      eval(parse( text=paste0( 'Plots[["', Tissue, '"]] <- Pr' ) ))
      
      ggsave( filename = paste0('DiffExpress_', Sex, '_',Tissue,'.eps'),
              plot = Pr,
              device = "eps",
              path = plotdir
      )
      
    }
  }
  return( list( Plots=Plots, Models=Models ) )
}

#Plot/Run/Build the DE Models

ALL_Mod <- ModelConstruct( DE, 'AD-CONTROL', 'Diagnosis', 'ALL' )
Male_Mod <- ModelConstruct( DE, NA, NA, 'MALE' )
Femal_Mod <- ModelConstruct( DE, NA, NA, 'FEMALE' )

REL <- DE[ (DE$ensembl_gene_id %in% Genes) == T , ]
colnames(REL)[ colnames(REL) == 'ensembl_gene_id'] <- 'ensemblgeneid'
colnames(REL)[ colnames(REL) == 'hgnc_symbol'] <- 'hgncsymbol'

ALL <- REL[ REL$Sex == 'ALL' & REL$Comparison == "AD-CONTROL" & REL$Model == 'Diagnosis',  c('hgncsymbol','ensemblgeneid','Tissue','logFC','adj.P.Val')] 
Males <-  REL[ REL$Sex == 'MALE' ,  c('hgncsymbol','ensemblgeneid','Tissue','logFC','adj.P.Val')]
Females <-  REL[ REL$Sex == 'FEMALE',  c('hgncsymbol','ensemblgeneid','Tissue','logFC','adj.P.Val')]

colnames(Tab)[ colnames(Tab) == 'Gene_Name' ] <- 'GeneName'
DE_Plots <- function( INPUT ){
  #Plots the table of the DE of the gene list
  #'@INPUT The input differential expression objects e.g ALL
  Gen <- as.character(INPUT$ensemblgeneid)[ !duplicated(as.character(INPUT$ensemblgeneid))]
  DF <- as.data.frame(matrix(0,length(Gen),9))
  row.names( DF ) <- Gen
  colnames( DF ) <- c( 'ENSG', 'GeneName', Tissues)
  DF$ENSG <- Gen
  DF$GeneName <- Tab[ DF$ENSG, ]$GeneName
  
  for( row in row.names(DF)){
    for( col in colnames(DF)[ colnames(DF) %in% c("CBE", "DLPFC", "FP", "IFG", "PHG", "STG", "TCX") ]){
      #Test if the the gene and tissue exist:
      if( row %in% as.character( INPUT[  INPUT$Tissue == col, ]$ensemblgeneid ) ){
        DF[ row, col ] <- signif( INPUT[ INPUT$ensemblgeneid == row & INPUT$Tissue == col, ]$logFC, 3)
      }else{
        DF[ row, col ] <- NA
      }
    }
  }
  
  #P Vals
  colnames(Tab)[ (colnames(Tab) %in% 'Gene_Name') == T ] <- 'GeneName'
  pDF <- as.data.frame(matrix(0,length(Gen),9))
  row.names( pDF ) <- Gen
  colnames( pDF ) <- c( 'ENSG', 'GeneName', Tissues)
  pDF$ENSG <- Gen
  pDF$GeneName <- Tab[ pDF$ENSG , ]$GeneName
  
  for( row in row.names(pDF)){
    for( col in colnames(pDF)[ colnames(pDF) %in% c("CBE", "DLPFC", "FP", "IFG", "PHG", "STG", "TCX") ]){
      if( row %in% as.character( INPUT[  INPUT$Tissue == col, ]$ensemblgeneid ) ){
        pDF[ row, col ] <- INPUT[ INPUT$ensemblgeneid == row & INPUT$Tissue == col, ]$adj.P.Val
      }else{
        pDF[ row, col ] <- 40
      }
    }
  }

  return( DF %>% 
            mutate(
              CBE = cell_spec( CBE, bold = ifelse(pDF$CBE < 0.05, TRUE, FALSE ), color = ifelse(pDF$CBE < 0.05, ifelse(DF$CB < 0, "red", "green"), rgb(211/255,211/255,211/255)) ) ,
              DLPFC = cell_spec(DLPFC, bold = ifelse(pDF$DLPFC < 0.05, TRUE, FALSE ),color = ifelse(pDF$DLPFC < 0.05, ifelse(DF$DLPFC < 0, "red", "green"), rgb(211/255,211/255,211/255) )),
              FP = cell_spec(FP, bold = ifelse(pDF$FP < 0.05, TRUE, FALSE ),color = ifelse(pDF$FP < 0.05, ifelse(DF$FP < 0, "red", "green"), rgb(211/255,211/255,211/255) )),
              IFG = cell_spec(IFG, bold = ifelse(pDF$IFG < 0.05, TRUE, FALSE ),color = ifelse(pDF$IFG < 0.05, ifelse(DF$IFG < 0, "red", "green"), rgb(211/255,211/255,211/255) )),
              PHG = cell_spec(PHG, bold = ifelse(pDF$PHG < 0.05, TRUE, FALSE ), color = ifelse(pDF$PHG < 0.05, ifelse(DF$PHG < 0, "red", "green"), rgb(211/255,211/255,211/255) )),
              STG = cell_spec(STG, bold = ifelse(pDF$STG < 0.05, TRUE, FALSE ),color = ifelse(pDF$STG < 0.05, ifelse(DF$STG < 0, "red", "green"), rgb(211/255,211/255,211/255) )),
              TCX = cell_spec(TCX, bold = ifelse(pDF$TCX < 0.05, TRUE, FALSE ),color = ifelse(pDF$TCX < 0.05, ifelse(DF$TCX < 0, "red", "green"), rgb(211/255,211/255,211/255) ))
            ) %>%
            kable( escape = F, row.names = FALSE) %>%
            kable_styling("striped", full_width = FALSE, position="left") )
}


DE_Plots(ALL)
save_kable(
  DE_Plots(ALL),
  file=paste0(tabledir,'/DiffExp_ALL.pdf')
)

lay <- rbind(c(1,2,3,4,7),
             c(1,2,5,6,7))

g1 <- ggplotGrob(ALL_Mod$Plots$CBE)
g2 <- ggplotGrob(ALL_Mod$Plots$DLPFC)
g3 <- ggplotGrob( ALL_Mod$Plots$FP)
g4 <- ggplotGrob(ALL_Mod$Plots$IFG)
g5 <- ggplotGrob(ALL_Mod$Plots$PHG)
g6 <- ggplotGrob( ALL_Mod$Plots$STG)
g7 <- ggplotGrob(ALL_Mod$Plots$TCX)

gs = list( g1, g2, g3, g4, g5, g6, g7 )

plot( grid.arrange(grobs = gs, layout_matrix = lay, left = "Calculated Weight", bottom = "Absolute Value of Log(Fold Change)") )

#### Female Cases Versus Controls
DE_Plots(Females)
save_kable(
  DE_Plots(Females),
  file=paste0(tabledir,'/DiffExp_Females.pdf')
)

lay <- rbind(c(1,2,3,4,7),
             c(1,2,5,6,7))

g1 <- ggplotGrob(Femal_Mod$Plots$CBE)
g2 <- ggplotGrob(Femal_Mod$Plots$DLPFC)
g3 <- ggplotGrob( Femal_Mod$Plots$FP)
g4 <- ggplotGrob(Femal_Mod$Plots$IFG)
g5 <- ggplotGrob(Femal_Mod$Plots$PHG)
g6 <- ggplotGrob( Femal_Mod$Plots$STG)
g7 <- ggplotGrob(Femal_Mod$Plots$TCX)

gs = list( g1, g2, g3, g4, g5, g6, g7 )

grid.arrange(grobs = gs, layout_matrix = lay)

#### Male Cases Versus Controls
DE_Plots(Males)
save_kable(
  DE_Plots(Males),
  file=paste0(tabledir,'/DiffExp_Males.pdf')
)

lay <- rbind(c(1,2,3,4,7),
             c(1,2,5,6,7))

g1 <- ggplotGrob(Male_Mod$Plots$CBE)
g2 <- ggplotGrob(Male_Mod$Plots$DLPFC)
g3 <- ggplotGrob( Male_Mod$Plots$FP)
g4 <- ggplotGrob(Male_Mod$Plots$IFG)
g5 <- ggplotGrob(Male_Mod$Plots$PHG)
g6 <- ggplotGrob( Male_Mod$Plots$STG)
g7 <- ggplotGrob(Male_Mod$Plots$TCX)

gs = list( g1, g2, g3, g4, g5, g6, g7 )

grid.arrange(grobs = gs, layout_matrix = lay)

