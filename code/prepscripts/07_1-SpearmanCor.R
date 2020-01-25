library(parallel)
library(doParallel)

ExpressionDS <- c('syn21291908','syn21292041','syn21285564','syn21285564','syn21285564','syn21285564','syn21291908')
names( ExpressionDS ) <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

for( Tissue in c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX') ){
  
  thisRepo <- githubr::getRepo(repository = "jgockley62/AD_TargetRank", ref="branch", refName='master')
  thisFile  <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/prepscripts/07_1-SpearmanCor.R')
  
  LIST <- colnames(x)
  LIST <- combn(LIST,2)
  
  OUT_Cor <- as.data.frame( matrix(0, length(colnames(x)), length(colnames(x))) )
  colnames(OUT_Cor) <- colnames(x)
  row.names(OUT_Cor) <- colnames(x)
  
  OUT_pVal <- as.data.frame( matrix(0, length(colnames(x)), length(colnames(x))) )
  colnames(OUT_pVal) <- colnames(x)
  row.names(OUT_pVal) <- colnames(x)
  
  library(doParallel)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  
  mark <- Sys.time()
  foreach(i=c(1:dim(LIST)[2])) %do% {
    
    mark <- Sys.time()
    foo<-LIST[1,i]
    foo2<-LIST[2,i]
    test <- cor.test( x[,foo], x[,foo2], method = "spearman" )
    
    OUT_Cor[ foo,foo2 ] <- test$estimate
    OUT_Cor[ foo2,foo ] <- test$estimate
    
    #pVal Test
    OUT_pVal[ foo,foo2 ] <- test$p.value
    OUT_pVal[ foo2,foo ] <- test$p.value
    Sys.time()-mark
  }
  Sys.time()-mark
  
  stopCluster(cl)
  
}