library(gplots)
library(superheat)

Genes <- Tab$ENSG
Nets <- read.csv(file = syn_temp$get('syn21534603')$path, header=T)
ZoomNet <- Nets[ (Nets$GeneID %in% Genes) == T, ]

NetPlot <- as.data.frame(matrix(0,length(Genes[!duplicated(Genes)]),5))
row.names(NetPlot) <- Genes[!duplicated(Genes)]
NetPlot$V1 <- row.names(NetPlot)
colnames(NetPlot) <- c("ENSG", "GeneName", "InACluster", "NumberOfClusters", "BrainRegions")

for( g in Genes ){
  if( dim(ZoomNet[ as.character(ZoomNet$GeneID) == NetPlot[g,]$ENSG, ])[1] > 0 ){
    NetPlot[g,]$GeneName <- names( table( as.character( ZoomNet[ as.character(ZoomNet$GeneID) == NetPlot[g,]$ENSG, ]$external_gene_name ) ) )
    NetPlot[g,]$InACluster <- 'YES'
    tmp <- ZoomNet[ as.character(ZoomNet$GeneID) == NetPlot[g,]$ENSG, ]
    NetPlot[g,]$NumberOfClusters <- length(tmp$Module)
    Regions <- c(as.character(tmp$brainRegion))
    NetPlot[g,]$BrainRegions <- paste( Regions[ !duplicated(Regions) ], collapse = ', ')
  }else{
    NetPlot[g,]$GeneName <- as.character( TotP[ TotP$ENSG == g, ]$GeneName )
    NetPlot[g,]$InACluster <- 'NO'
    NetPlot[g, 3:5] <- NA
  }
}

NetPlot %>% 
  select(ENSG, GeneName, InACluster, NumberOfClusters, BrainRegions)  %>%
  kable( escape = F, row.names = FALSE) %>%
  kable_styling("striped", full_width = F)

save_kable(
  NetPlot %>% 
    select(ENSG, GeneName, InACluster, NumberOfClusters, BrainRegions)  %>%
    kable( escape = F, row.names = FALSE) %>%
    kable_styling("striped", full_width = F),
  
  file=paste0(tabledir,'/InNetCluster.pdf')
)

Regions <- as.data.frame( table(ZoomNet$brainRegion) )
colnames(Regions) <- c("Region", "Targets")

Targets <-as.data.frame(Regions)
Targets %>% 
  select(Region, Targets)  %>%
  kable( escape = F, row.names = FALSE) %>%
  kable_styling("striped", full_width = F)

save_kable(
  Targets %>% 
    select(Region, Targets)  %>%
    kable( escape = F, row.names = FALSE) %>%
    kable_styling("striped", full_width = F),
  
  file=paste0(tabledir,'/Targets_by_Region.pdf')
)

Mods <- as.character(ZoomNet$Module)
Mods <- Mods[ !duplicated(Mods) ]

dim( table( ZoomNet[ ZoomNet$brainRegion == 'DLPFC', ]$ModuleName,  ZoomNet[ ZoomNet$brainRegion == 'DLPFC', ]$external_gene_name ) )


HM_Mat <- matrix( 0, length(Genes), length(Mods) )
row.names(HM_Mat) <- as.character(NetPlot$GeneName)
colnames(HM_Mat) <- Mods

trans <- NetPlot$GeneName
names(trans) <- NetPlot$ENSG

for( i in 1:dim(ZoomNet)[1] ){
  HM_Mat[ trans[as.character( ZoomNet[ i, ]$GeneID )], 
          as.character( ZoomNet[ i, ]$Module ) ] <- 1
}

PAL <- colorRampPalette(c( "white", "red"))(n = 2)
heatmap.2(as.matrix(HM_Mat), col=PAL, key = FALSE, tracecol=NA, cexRow = .9, cexCol=.75) 

setEPS()
postscript(file = paste0(plotdir,'/ModuleHeatMap.eps') )
heatmap.2(as.matrix(HM_Mat), col=PAL, key = FALSE,
          tracecol=NA, cexRow = .9, cexCol=.75)
dev.off()

#Figure 2A: https://www.biorxiv.org/content/biorxiv/early/2019/01/03/510420.full.pdf
Consensus_Cluster_A <- c("TCXblue", "IFGyellow", "PHGyellow")
Consensus_Cluster_B <- c("DLPFCblue", "CBEturquoise", "STGblue", "PHGturquoise", "IFGturquoise", "TCXturquoise", "FPturquoise")
Consensus_Cluster_C <- c("IFGbrown", "STGbrown", "DLPFCyellow", "TCXgreen", "FPyellow", "CBEyellow", "PHGbrown" )
Consensus_Cluster_D <- c("DLPFCbrown", "STGyellow", "PHGgreen", "CBEbrown", "TCXyellow", "IFGblue", "FPblue")
Consensus_Cluster_E <- c("FPbrown", "CBEblue", "DLPFCturquoise", "TCXbrown", "STGturquoise", "PHGblue")

ConClust <- as.data.frame( cbind( c(Consensus_Cluster_A, Consensus_Cluster_B, Consensus_Cluster_C, Consensus_Cluster_D, Consensus_Cluster_E),
                                  c( rep("A", length(Consensus_Cluster_A)) ,rep("B", length(Consensus_Cluster_B)), rep("C", length(Consensus_Cluster_C)),
                                     rep("D", length(Consensus_Cluster_D)), rep("E", length(Consensus_Cluster_E))) ))

AllGenes <- Nets$GeneID[ !duplicated(Nets$GeneID) ]

colnames(ConClust) <- c("Module", "ConesusCluster")
row.names(ConClust) <- ConClust$Module

#-#TempClust <- as.data.frame( matrix(0, length(Genes), 5))
TempClust <- as.data.frame( matrix(0, length(AllGenes), 5))
#-#row.names(TempClust) <- Genes
row.names(TempClust) <- AllGenes
colnames(TempClust) <- c("A", "B", "C", "D", "E")

for( i in AllGenes){
  #tmp_cs <- table( as.character( ConClust[ as.character(ZoomNet[ ZoomNet$GeneID == i, ]$Module), ]$ConesusCluster ) )
  tmp_cs <- table( as.character( ConClust[ as.character(Nets[ Nets$GeneID == i, ]$Module), ]$ConesusCluster ) )
  for(y in names(tmp_cs)){
    TempClust[ i, y] <- as.numeric(tmp_cs[y])
  }
}

#-#row.names(TempClust) <- as.character( trans[ row.names(TempClust)] )
tNets <- Nets[ !duplicated(Nets$GeneID), ]
row.names(tNets) <- tNets$GeneID

#Fix the repeats...
PAL <- colorRampPalette(c(  "white", "yellow", "red"))(n = 10)

Clus <- TempClust[names(trans),]
row.names(Clus) <- trans[row.names(Clus)]
HM <- heatmap.2(as.matrix(Clus), col=PAL, key = FALSE, tracecol=NA, cexRow = .9, cexCol=.75)

Tempclust <- TempClust[ rownames(TempClust) %in% Genes, ] 
row.names(Tempclust) <- trans[ row.names(Tempclust) ]

save_kable(
  Tempclust %>% 
    select( A, B, C, D, E)  %>%
    kable( escape = F, row.names = TRUE) %>%
    kable_styling("striped", full_width = F),
  file=paste0(tabledir,'/ConsensousCluster.pdf')
)

Tempclust %>% 
  select( A, B, C, D, E)  %>%
  kable( escape = F, row.names = TRUE) %>%
  kable_styling("striped", full_width = F)
