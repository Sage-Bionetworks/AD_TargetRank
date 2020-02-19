eqtl <- read.csv( file =syn_temp$get('syn21534584')$path, header=T)

Genes <- Tab$ENSG

uEQs <- eqtl[ (eqtl$ensembl_gene_id %in% Genes) == T,]

uEQs$hasEqtl <- as.character(uEQs$hasEqtl)
Rep <- c("NO", "YES")
names(Rep) <- c('FALSE', 'TRUE')
uEQs$hasEqtl <- as.character(Rep[ as.character(uEQs$hasEqtl) ])

uEQs %>%
  select(ensembl_gene_id, hgnc_symbol, hasEqtl)  %>%
  kable( escape = F, row.names = FALSE ) %>%
  kable_styling("striped", full_width = F)

save_kable(
  uEQs %>%
    select(ensembl_gene_id, hgnc_symbol, hasEqtl)  %>%
    kable( escape = F, row.names = FALSE ) %>%
    kable_styling("striped", full_width = F),
  
  file=paste0(tabledir,'/eQTL_inBrain.pdf')
)
