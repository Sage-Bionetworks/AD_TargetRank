
Lake <- read.csv( file=syn_temp$get('syn21534588')$path, header = T )
Lake <- Lake[ Lake$adj.pval < 0.05, ]
Lake$Module <- as.character(Lake$Module)
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
