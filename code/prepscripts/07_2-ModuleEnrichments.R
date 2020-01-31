library(parallel)
library(doParallel)
library(data.table)

reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()
#syn21554905

Its <- c( 'M62', 'EmoryTargets', 'Moesin', 'pQTL', 'Retromer' )

Cors <- c( 'syn21544714', 'syn21546998', 'syn21546838', 'syn21546262', 'syn21546567' )
names(Cors) <- c( 'M62', 'EmoryTargets', 'Moesin', 'pQTL', 'Retromer' )

Lists <- c('genelists/ModuleM62.txt', 
           'genelists/EmoryTargets.txt', 
           'genelists/Moesin.txt',
           'genelists/pQTL.txt',
           'genelists/Retromer.txt')
names(Lists) <- c( 'M62', 'EmoryTargets', 'Moesin', 'pQTL', 'Retromer' )

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

AvgParCor <- list() 
altAvgParCor <- list() 
for( i in Its){
  Genes <- as.character( read.table( Lists[i], header=F, sep='\t' )[,1] )

  results = read.csv(syn_temp$tableQuery(paste0( "select * from syn21554905 where Target_GeneName in ('", 
                                               paste(Genes,collapse="','") ,
                                               "') and Rank < 21"))$filepath)[,c("Target_Gene", 
                                                                                 "Target_GeneName", 
                                                                                 "Hit_Gene", 
                                                                                 "Hit_GeneName", 
                                                                                 "Rank", 
                                                                                 "set")]
  Keeps <- c(paste0( as.character(results$Target_Gene), '_', as.character(results$Hit_Gene)))
  Keeps <- Keeps[!duplicated(Keeps)] 
  writeLines(paste0( "Extended Module list ", i,
                     ' contains ', 
                     length( Keeps ),
                     ' genes.'))
  
  
  Partials <- read.table( syn_temp$get( Cors[i] )$path, header = T, sep ='' )
  #Keeps <- as.character(results$Hit_Gene)[!duplicated(as.character(results$Hit_Gene))]
  #Partials <- Partials[ as.character(Partials$Hit_Gene) %in% as.character(results$Hit_Gene), ]
  
  Partials$Comparison <- paste0( as.character(Partials$Target_Gene), '_', as.character(Partials$Hit_Gene) )
  
  Saves <- Partials[ as.character(Partials$Comparison) %in% Keeps, ]
  #Saves[ as.character(Saves$Comparison) == 'ENSG00000117592_ENSG00000162407', ]
  
  AvgPars <- as.data.frame( matrix( NA, length(Keeps), dim(Partials)[2] ) )
  colnames(AvgPars) <- colnames(Partials)
  AvgPars$Comparison <- Keeps
  row.names(AvgPars) <- AvgPars$Comparison
  
  #for( temp in row.names(AvgPars) ){
  RUNNEr <- function( temp, Partials){
      tempM <- Partials[ Partials$Comparison == temp, ]
      ROW <- c( as.character(tempM$Target_Gene)[1], as.character(tempM$Target_GeneName)[1], as.character(tempM$Seed_Set)[1],
                         paste(as.character(tempM$Hit_Tissue), collapse=','), as.character(tempM$Hit_Data_Set)[1], as.character(tempM$Hit_Gene)[1], 
                         as.character(tempM$Hit_GeneName[1]), mean(tempM$Hit_Probability), mean(tempM$Spearman_Correlation), temp
                       )
      return(ROW)
  }
  
  mark <- Sys.time()
  x <- foreach(i=Keeps[1:length(Keeps)], .combine=rbind) %dopar% RUNNEr( i, Partials=Partials)
  Sys.time() - mark
  
  #Export the DF
  row.names(x) <- x[,10]
  colnames(x) <- colnames(Partials)
  x<-as.data.frame(x)
  x$Hit_Data_Set <- NA
  colnames(x)[8] <- 'Avg_Hit_Probability'
  colnames(x)[9] <- 'Avg_Spearman_Correlation'
  eval( parse( text= paste0( 'AvgParCor$', i, ' <- x' ))) 
  eval( parse( text= paste0( 'write.table(x, file=\'', i, '_AvgCors.tsv\', row.names=F, col.names=F, quote=F, sep=\'\t\' )')))
  
  #Boxplot:
  X <- x[ (as.character(x$Target_GeneName) == as.character(x$Hit_GeneName)) == F , ] 
  summary( as.numeric(as.character(X$Avg_Hit_Probability )))
  eval( parse( text= paste0( 'altAvgParCor$', i, ' <- X' ))) 
  SUM <- summary(as.numeric(as.character(X$Avg_Hit_Probability)))
  writeLines(paste0( "Partrial Correlation Averages across Tissue for ", i,
                     ' Mean: ', 
                     signif(SUM[4], digits = 3) ,
                     ' IQR(', signif(SUM[2], digits = 3),'-', signif(SUM[5], digits = 3),')' ))
  
}


PLT <- as.data.frame( rbind( altAvgParCor$M62, 
                             altAvgParCor$EmoryTargets, 
                             altAvgParCor$Moesin,
                             altAvgParCor$pQTL,
                             altAvgParCor$Retromer
                      ))

library(ggplot2)
PLT$Avg_Hit_Probability <- as.numeric(as.character(PLT$Avg_Hit_Probability))
ggplot(PLT, aes(x=Seed_Set, y=Avg_Hit_Probability, col=Seed_Set)) +
  geom_violin() + geom_boxplot(width=0.1) + theme(legend.position = "none") +
  labs(title = "Distribution of Average Hit Probabilities Across 7 Tissues") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Seed Set") +
  ylab("Average Hit Probability") 

# Enrichr mechanism enrichment
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(enrichR)
library(stringi)
library(gridExtra)

dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), Space)
#dbs <- c("Drug_Perturbations_from_GEO_2014", "Tissue_Protein_Expression_from_ProteomicsDB","Drug_Perturbations_from_GEO_down", "Reactome_2016", "KEGG_2016", "WikiPathways_2016", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "WikiPathways_2019_Human", "KEGG_2019_Human")                                  
dbs <- c("Drug_Perturbations_from_GEO_2014", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "WikiPathways_2019_Human", "KEGG_2019_Human")                                  
#ModuleName<-'Retromer'
EnRICH <- function( ModuleName ){
  
  enriched <- eval(parse(text=paste0('enrichr( as.character(altAvgParCor$', ModuleName ,'$Hit_GeneName), dbs)')))
  
  #enriched$Drug_Perturbations_from_GEO_2014$Adjusted.P.value[1] < 0.05
  #enriched$GO_Biological_Process_2018$Adjusted.P.value[1] < 0.05
  #enriched$GO_Cellular_Component_2018$Adjusted.P.value[1] < 0.05
  #enriched$WikiPathways_2019_Human$Adjusted.P.value[1] < 0.05
  #enriched$KEGG_2019_Human$Adjusted.P.value[1] < 0.05
  
  
  plts <- c("KEGG_2019_Human", "WikiPathways_2019_Human", 
            "GO_Cellular_Component_2018", "GO_Biological_Process_2018"
          )
  
  Lyst <- list()
  for( Path in plts[1:length(plts)]){
    if( eval(parse(text= paste0( 'dim(enriched$',Path,'[ enriched$',
                                 Path,'$Adjusted.P.value < 0.05, ])[1] > 0'))) ){
      foo <- eval(parse(text= paste0( 'enriched$',Path,'[ enriched$',
                                      Path,'$Adjusted.P.value < 0.05, ]' )))
      foo$term <- NA
      for( i in 1:dim(foo)[1]){
        foo[i,]$term <- stri_reverse(stri_split_fixed(stri_reverse(foo[i,]$Term)," ",n = 2)[[1]])[2]
      }
      foo$term <- gsub('Homo sapiens', '', foo$term) 
      
      foo$LogP <- -1*log(foo$Adjusted.P.value)
      foo$rank <- as.numeric(c(dim(foo)[1]:1) )
    
      POL <- do.call( rbind, strsplit(foo$Overlap, '/'))
      foo$PercOverLap <- as.numeric(POL[,1])/as.numeric(POL[,2])
      foo$term <- paste0(foo$term, " N=", as.character(POL[,2]) )
      p <- ggplot( data=foo, aes(y=rank, x=LogP, size=Odds.Ratio, col=PercOverLap) ) + geom_point(  ) +
          scale_y_discrete(name = element_blank(), limits = foo$term ) + 
          xlab('-Log( Adj_PVal )') + labs(title = Path) + 
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_color_gradient(low="blue", high="red")
      eval(parse(text=paste0( 'Lyst$', Path, ' <- p '  )))
    }else{}
  }
  return( Lyst)
}

m62 <- EnRICH( 'M62')
EmTar <- EnRICH( 'EmoryTargets')
Moe <- EnRICH( 'Moesin')
pqtl <- EnRICH( 'pQTL')
Ret <- EnRICH( 'Retromer')

grid.arrange(Ret$KEGG_2019_Human, 
             #Ret$WikiPathways_2019_Human, 
             nrow = 2
            )
grid.arrange(#Ret$GO_Cellular_Component_2018,
             Ret$GO_Biological_Process_2018,
             nrow = 2
            )

#########################################

Lake <- read.table( file=syn_temp$get('syn21563069')$path, sep='\t', header = T )
tGenes <- as.character(Partials$Hit_GeneName)[ !duplicated(as.character(Partials$Hit_GeneName)) ]

CellTypes <- table(Lake$Cluster)

PermNum <- c( length( as.character(altAvgParCor$M62$Hit_GeneName)[!duplicated(as.character(altAvgParCor$M62$Hit_GeneName))] ) ,
              length( as.character(altAvgParCor$EmoryTargets$Hit_GeneName)[!duplicated(as.character(altAvgParCor$EmoryTargets$Hit_GeneName))] ) ,
              length( as.character(altAvgParCor$Moesin$Hit_GeneName)[!duplicated(as.character(altAvgParCor$Moesin$Hit_GeneName))] ) ,
              length( as.character(altAvgParCor$pQTL$Hit_GeneName)[!duplicated(as.character(altAvgParCor$pQTL$Hit_GeneName))] ) ,
              length( as.character(altAvgParCor$Retromer$Hit_GeneName)[!duplicated(as.character(altAvgParCor$Retromer$Hit_GeneName))] ) 
            )
names(PermNum) <- c( 'M62', 'EmoryTargets', 'Moesin', 'pQTL', 'Retromer' )

COLz <- c( paste0( names(PermNum)[1], "_", names(CellTypes) ), 
                      paste0( names(PermNum)[2], "_", names(CellTypes) ),
                      paste0( names(PermNum)[3], "_", names(CellTypes) ),
                      paste0( names(PermNum)[4], "_", names(CellTypes) ),
                      paste0( names(PermNum)[5], "_", names(CellTypes) )
                     )
Perms <- as.data.frame(matrix(0, 10000, length(COLz)))
colnames(Perms) <- COLz

mark <- Sys.time()
for( i in 1:10000){
  perm <- NULL
  for( name in names(PermNum) ){
    Samp <- sample( tGenes, as.numeric(PermNum[name]), replace = FALSE)
    for( Name in names(CellTypes) ){
      perm <- eval(parse( text= paste0( 
        'c(perm, ',name, '_', Name,' = length( Samp[ Samp %in% Lake[ Lake$Cluster == Name, ]$Gene ]) )' 
        )
      ))
      
    }
    
  }
  Perms[i, ] <- perm
}
Sys.time() - mark

#colnames(Perms) <- gsub("\\s+", '', colnames(Perms))
MEANz <- apply(Perms, 2, mean)
SDz <- apply(Perms, 2, sd)

GeneLst <- list( M62 = as.character(altAvgParCor$M62$Hit_GeneName)[!duplicated(as.character(altAvgParCor$M62$Hit_GeneName))] ,
                    EmoryTargets = as.character(altAvgParCor$EmoryTargets$Hit_GeneName)[!duplicated(as.character(altAvgParCor$EmoryTargets$Hit_GeneName))],
                    Moesin = as.character(altAvgParCor$Moesin$Hit_GeneName)[!duplicated(as.character(altAvgParCor$Moesin$Hit_GeneName))],
                    pQTL = as.character(altAvgParCor$pQTL$Hit_GeneName)[!duplicated(as.character(altAvgParCor$pQTL$Hit_GeneName))],
                    Retromer = as.character(altAvgParCor$Retromer$Hit_GeneName)[!duplicated(as.character(altAvgParCor$Retromer$Hit_GeneName))] 
)

ActualsLake <- as.data.frame(matrix( 0, length(GeneLst), length(CellTypes) ))
row.names(ActualsLake) <- names(GeneLst)
colnames(ActualsLake) <- names(CellTypes)

ROW <- 'M62'
COL <- 'Ast'
for( ROW in row.names(ActualsLake)){
  for( COL in colnames(ActualsLake)){
    Test <- eval(parse(text= paste0( 'GeneLst$', ROW )))
    CTyp <- Lake[ Lake$Cluster == COL, ]$Gene
    ActualsLake[ ROW,COL ] <- length( Test[ Test %in% CTyp] )
  }
}

ZVAL <- function( Obs, M, S ){
  return( (Obs - M)/S )
}

ZsLake <- as.data.frame(matrix( 0, length(GeneLst), length(CellTypes) ))
row.names(ZsLake) <- names(GeneLst)
colnames(ZsLake) <- names(CellTypes)

for( ROW in row.names(ZsLake)){
  for( COL in colnames(ZsLake)){
    ZsLake [ ROW,COL ] <- as.numeric( 
      ZVAL( ActualsLake[ ROW,COL ], MEANz[ paste0(ROW, '_', COL) ], SDz[ paste0(ROW, '_', COL) ] )
    )
  }
}

PsLake <- as.data.frame(matrix( 0, length(GeneLst), length(CellTypes) ))
row.names(PsLake) <- names(GeneLst)
colnames(PsLake) <- names(CellTypes)
pvalue2sided <- function(x){ 2*pnorm(-abs(x)) }

for( ROW in row.names(PsLake)){
  for( COL in colnames(PsLake)){
    PsLake[ ROW,COL ] <- pvalue2sided( ZsLake[ ROW,COL ] )
  }
}

CorPsLake <- as.data.frame(matrix( 0, length(GeneLst), length(CellTypes) ))
row.names(CorPsLake) <- names(GeneLst)
colnames(CorPsLake) <- names(CellTypes)
for( ROW in row.names(PsLake)){
    CorPsLake[ ROW, ] <- p.adjust(PsLake[ ROW, ], method='fdr', n=length(PsLake[ ROW, ]))
}

StarPsLake <- as.data.frame(matrix( "", length(GeneLst), length(CellTypes) ))
row.names(StarPsLake) <- names(GeneLst)
colnames(StarPsLake) <- names(CellTypes)
for( ROW in row.names(PsLake)){
  for( COL in colnames(PsLake)){
    if( CorPsLake[ ROW,COL ] < 0.05 ){
      StarPsLake[ ROW,COL ] <- "*"
    }  
  }
}

library(pheatmap)
pheatmap( ZsLake, cluster_rows=F, cluster_cols = F, display_numbers=StarPsLake)
pheatmap( ZsLake, display_numbers=StarPsLake)


Zhang <- read.csv( file=syn_temp$get('syn21534606')$path, header = T )
Zhang <- Zhang[ Zhang$adj.pval < 0.05, ]
Zhang$Module <- as.character(Zhang$Module)


Zrich <- as.data.frame( matrix(0, length(Genes), length(table(Zhang$category)) ))
row.names(Zrich) <- trans[Genes]
colnames(Zrich) <- names(table(Zhang$category))

Lrich <- as.data.frame( matrix(0, length(Genes), length(table(Lake$category)) ))
row.names(Lrich) <- trans[Genes]
colnames(Lrich) <- names(table(Lake$category))


for( gene in Genes ){
  lyst <- as.character( ZoomNet[ ZoomNet$GeneID == gene, ]$Module )
  foo <- Lake[ as.character(Lake$Module) %in% lyst ==T , ]
  for( Cell in as.character(foo$category) ){
    Lrich[ trans[gene], Cell  ] <- Lrich[ trans[gene], Cell  ] + 1
  }
  foo <- Zhang[ as.character(Zhang$Module) %in% lyst ==T , ]
  for( Cell in as.character(foo$category) ){
    Zrich[ trans[gene], Cell  ] <- Zrich[ trans[gene], Cell  ] + 1
  }
}
save_kable(
  Zrich %>% 
    kable( escape = F, row.names = TRUE) %>%
    kable_styling("striped", full_width = F),
  file=paste0(tabledir,'/ZhangSingleCell.pdf')
)
save_kable(
  Lrich %>% 
    kable( escape = F, row.names = TRUE) %>%
    kable_styling("striped", full_width = F),
  file=paste0(tabledir,'/LakeSingleCell.pdf')
)
Zrich %>% 
  kable( escape = F, row.names = TRUE) %>%
  kable_styling("striped", full_width = F)

Lrich %>% 
  kable( escape = F, row.names = TRUE) %>%
  kable_styling("striped", full_width = F)

PAL <- colorRampPalette(c( "white", "yellow", "red"))(n = 50)


setEPS()
postscript(file = paste0(plotdir,'/ZhangCellTypeHeatMap.eps') )
heatmap.2(as.matrix(Zrich), col=PAL, key = FALSE,
          tracecol=NA, cexRow = .9, cexCol=.75)
dev.off()

setEPS()
postscript(file = paste0(plotdir,'/LakeCellTypeHeatMap.eps') )
heatmap.2(as.matrix(Lrich), col=PAL, key = FALSE,
          tracecol=NA, cexRow = .9, cexCol=.75)
dev.off()

heatmap.2(as.matrix(Zrich), col=PAL, key = FALSE,
          tracecol=NA, cexRow = .9, cexCol=.75)

heatmap.2(as.matrix(Lrich), col=PAL, key = FALSE,
          tracecol=NA, cexRow = .9, cexCol=.75)

