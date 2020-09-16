#Final Code to develop DE Proteomic Analysis for Agora
#Note: Output will vary slightly between runs as a result of the bootstraping regression
##Special Thanks to Eric Dammer and the Bioinformatics Team at Emory for an extrodinary amount of help

library("readxl")
library(sva)
library(impute)
library(variancePartition)
library(doParallel)
library(affy)
library(impute)
library(preprocessCore)
library("doParallel")
library(boot)
library(WGCNA)
#install.packages("BiocManager")
#BiocManager::install("WGCNA")
#sudo apt-get install libjpeg-dev
library(Cairo)
#Woof this Sucks
#sudo apt-get install r-base-dev, sudo apt-get install libxt-dev
#sudo apt-get install xfonts-base, sudo apt-get install xauth, sudo apt-get install libjpeg-dev
#sudo apt-get install xvfb, sudo apt-get install libcairo2-dev, sudo apt-get install libgtk2.0-dev
library(NMF)
#install.packages('NMF')
library(biomaRt)
library(variancePartition)
#BiocManager::install("variancePartition")
library(limma)

reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

Process <- function( FILE, Opt ){
  ##Processes the Raw data and collapse to GeneName|UniProtID 
  #'@FILE matrix of raw protien.groups file
  #'@Opt a one or a zero 
  
  options(stringsAsFactors = FALSE)
  proteinGroups <- FILE #FILE
  
  decoyCount = length(which(proteinGroups$Reverse=="+"))
  originalRows = nrow(proteinGroups)
  decoyIndices = which(proteinGroups$Reverse=="+")
  FDR = paste0(round(length(decoyIndices)/nrow(proteinGroups)*100, 2), "% FDR")
  LFQindices = which(grepl("LFQ.intensity.", colnames(proteinGroups)))
  cat(paste0("Imported data has ", originalRows, "x", ncol(proteinGroups), " rows x columns; ", decoyCount, " reverse hits, for a net ", FDR, ".\n", "LFQ intensity is available for ", length(LFQindices), " experiment samples.\n"))
  
  #make a matrix of just LFQ values with rownames as UniqueID (decoys removed, so they are truly unique now)
  exprMat0 <- proteinGroups[ -decoyIndices, c(1,LFQindices) ]
  
  #Split Out Preferred Uniprot ID preference:(sp| > tr| > other)
  IDdf <- as.data.frame( do.call( rbind, strsplit(as.character( exprMat0$Protein.IDs ), ";") ) )
  IDdf2 <- apply(IDdf, 1, function(x) if( sum(as.numeric(grepl("sp\\|", x))) > 0 ) {
    x[which(grepl("sp\\|",x))[1]] } 
    else { 
      if( sum(as.numeric(grepl("tr\\|", x))) > 0 ) { 
        x[which(grepl("tr\\|", x))[1]] 
      } 
      else { x[1] } 
    } 
  )
  
  IDdf2 <- apply( data.frame(col1=IDdf2), 1, function(x) gsub("sp\\|", "", gsub("tr\\|", "", x)) )
  
  uniprotIDs.preferred <- as.data.frame( do.call(rbind, strsplit(IDdf2, "[|]")) )[,1]
  
  #Lookup HGNC Gene Symbols
  lookupUniprotSymbol <- syn_temp$get('syn18914942')
  lookupUniprotSymbol <- read.csv(lookupUniprotSymbol$path, header = T )
  
  symbols.lookup <- lookupUniprotSymbol[ match(uniprotIDs.preferred, lookupUniprotSymbol[,1]), 2 ]
  symbols.lookup[is.na(symbols.lookup)] <- ""
  
  #make concatenated UniqueID for row.names of exprMat0
  uniqueIDs<-paste0( symbols.lookup, "|", uniprotIDs.preferred )
  
  #Set rownames to UniqueIDs and remove original column1
  exprMat0<-exprMat0[ , !colnames(exprMat0)=="Protein.IDs" ]
  
  rownames(exprMat0)<-uniqueIDs
  exprMat0<-as.matrix(exprMat0)
  dim(exprMat0)
  colnames(exprMat0) <- gsub( "LFQ.intensity.", "", colnames(exprMat0) )
  
  #clean up data structures we will not use again
  rm(proteinGroups)
  rm(uniqueIDs)
  rm(IDdf)
  rm(IDdf2)
  
  print(dim(exprMat0))
  
  if( Opt == 0 ){
    exprMat0 <- exprMat0[ -grep( "CON__", row.names(exprMat0) ), ]
    exprMat0 <- exprMat0[ -grep( "0\\|", substr(row.names(exprMat0) , 1, 2) ), ]
    exprMat0 <- exprMat0[ -grep( "\\|", substr(row.names(exprMat0) , 1, 1) ), ]
    
  }else{}
  
  exprMat0 <- exprMat0[ -grep( "0\\|", substr(row.names(exprMat0) , 1, 2) ), ]
  exprMat0 <- exprMat0[ -grep( "\\|", substr(row.names(exprMat0) , 1, 1) ), ]
  
  #Remove rows with over half missing values 
  row_sub = apply(exprMat0, 1, function(row) all( median(as.numeric(row)) != 0 ))
  exprMat0 <- exprMat0[ row_sub, ]
  
  return( as.data.frame( exprMat0 ) )
}

Logger_KNN <- function( Inpt, File ){
  ##Processes the GeneName|UniProt collapsed counts to log2 and imputes
  #'@Inpt the input dataframe eg. Banner_CDat
  #'@File the file name to save eg. BANNER
  foo <- log2(Inpt)
  foo2 <- as.data.frame( apply( foo, 2, function(x){ gsub(-Inf, NA, x) }) )
  row.names(foo2) <- row.names(foo)
  message(paste0( dim(foo2[ complete.cases(foo2),])[1] ))
  exprMatIMP <- impute.knn(as.matrix(foo2))
  exprMat.VAL <- exprMatIMP$data
  message(paste0( dim(exprMat.VAL[ complete.cases(exprMat.VAL),])[1] ))
  return( exprMat.VAL ) 
}

Runner <- function( exprMat.val, CombatInfo, Sample){
  ##Batch Corrects with COMBAT and Bootstrap Regression Covariates
  #'@exprMat.val Expression Matrix
  #'@CombatInfo  Meta Data Matrix
  #'@Sample Sample type
  #exprMat.val <-MSBB_RDY
  #CombatInfo <-MSBB_Pheno
  #Sample <- "MSBB"
  #exprMat.val <- exprMat.VAL
  #CombatInfo <- Combat_Info
  #Sample <- "BANNER"
  
  #exprMat.val <- Mayo_RDY
  #CombatInfo <- Mayo_Met
  #Sample <- "MAYO"
  
  
  if( Sample == "BANNER" ){
    #CombatInfo$Sample <- paste0( do.call( rbind, strsplit( CombatInfo$Sample,'_' ))[,1], '_',  do.call( rbind, strsplit( CombatInfo$Sample,'_' ))[,2])
    CombatInfo <- CombatInfo[ (CombatInfo$Sample %in% colnames(exprMat.val)) == T, ]
    row.names(CombatInfo) <- CombatInfo$Sample
  }
  CombatInfo <- CombatInfo[match(colnames(exprMat.val),rownames(CombatInfo)),]
  
  #<OPTIONAL>#
  #=======================================================================================#
  #Use Variance Partition to estimate variance explained by each covariates Before Combat #
  #=======================================================================================#
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  CombatInfo$Age <- as.numeric(CombatInfo$Age)
  
  if(Sample == "MAYO"){
    form <- ~ Age + (1|Batch)+Braak+(1|Gender)+PMI #+APOE
  }else{
    if(Sample == "BLSA"){
      form <- ~ Age +Braak+CERAD+(1|Gender)+PMI #+APOE
    }else{
      form <- ~ Age + (1|Batch)+Braak+CERAD+(1|Gender)+PMI #+APOE
    }
  }
  
  #_ NA's Dummy Vals _# CombatInfo[ grepl('gis', CombatInfo$Sample) == T, c(2,3,4,7:12)] <- NA
  
  CombatInfo$Sample <- factor(CombatInfo$Sample)
  CombatInfo$Gender <- factor(CombatInfo$Gender)
  if(Sample == "BLSA"){
    CombatInfo <- cbind( CombatInfo, Batch =1)
  }else{
    CombatInfo$Batch<-factor(CombatInfo$Batch)
  }
  
  varPart <- fitExtractVarPartModel( exprMat.val, form, CombatInfo)
  # sort variables (i.e. columns) by median fraction of variance explained
  vp <- sortCols( varPart )
  # Bar plot of variance fractions for the first 10 genes 
  plotPercentBars( vp[1:10,] )
  # violin plot of contribution of each variable to total variance 
  plotVarPart( vp )
  #}  
  colnames(CombatInfo)[ colnames(CombatInfo) == "Diagnosis" ] <- 'SampleForCombat'
  CombatInfo$SampleForCombat <- factor(CombatInfo$SampleForCombat)
  
  #Assemble Model Matrix for Combat
  model = model.matrix(~CombatInfo$SampleForCombat, data=as.data.frame(exprMat.val))
  if( Sample == "BLSA"){
    dat <- exprMat.val
    cleanDatImputed <- dat
  }else{
    if(Sample == "BANNER" || Sample == "MAYO" || Sample == "MSBB" ){
      exprMat.combat <- ComBat(dat=as.matrix(exprMat.val), batch=as.vector(CombatInfo$Batch), mod=model)
      dat <- exprMat.combat
      cleanDatImputed <- dat
    }
    else{
      exprMat.combat <- ComBat(dat=as.matrix(exprMat.val), batch=as.vector(CombatInfo$Batch), mod=model)
      dat <- exprMat.combat
      cleanDatImputed <- dat[ , (c(1:dim(dat)[2]) %in% c(grep('gis',colnames(dat)))) == F]
    }
  }
  cov <- CombatInfo
  
  rownames(cleanDatImputed)<-rownames(dat)
  #*MATCH CASE ORDER
  numericMeta<-cov[ match( colnames(cleanDatImputed)[1:ncol(cleanDatImputed)],cov$Sample ),]
  
  cleanDat <- cleanDatImputed
  
  numericMetaSimple <- numericMeta[ match( colnames(cleanDat),as.character(numericMeta$Sample) ),c(1:dim(numericMeta)[2])]
  
  if( Sample == "BLSA"){
    numericMetaSimple <- numericMetaSimple[ complete.cases(numericMetaSimple), ]
  }else{}
  rownames(numericMetaSimple) <- numericMetaSimple$Sample
  
  if( Sample == 'MAYO' || Sample == 'BANNER' || Sample == 'BLSA' ){
    numericMetaSimple <- numericMetaSimple[,-1]
  }else{
    numericMetaSimple <- numericMetaSimple[,-c(1,2)]
  }
  
  ##########################BOOTSTRAP REGRESSION########################################
  parallelThreads=10
  #Run Local
  clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
  registerDoParallel(clusterLocal)
  
  cleanDat.unreg <- cleanDat
  
  options(stringsAsFactors=FALSE)
  #library(boot)
  boot <- TRUE
  numboot <- 1000
  bs <- function(formula, data, indices) {
    d <- data[indices,] # allows bootstrap function to select samples
    fit <- lm(formula, data=d)
    return(coef(fit))
  }
  
  ## Get the covariate data
  condition <- as.numeric(factor(numericMetaSimple$AD))
  condition.AD <- as.numeric(condition==1)
  
  age = as.numeric(numericMetaSimple$Age)
  sex = as.numeric(factor(numericMetaSimple$Gender))-1
  PMI = as.numeric(numericMetaSimple$PMI)
  
  regvars <- as.data.frame(cbind(condition.AD, age,sex,PMI))
  regvars$condition.AD <- as.factor(regvars$condition.AD)
  #condition.ad,condition.ftdu,condition.tau,condition.pdd,condition.msa,condition.pd,condition.als
  
  ## Run the regression
  normExpr.reg <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(cleanDat))
  rownames(normExpr.reg) <- rownames(cleanDat)
  colnames(normExpr.reg) <- colnames(cleanDat)
  
  ## change it to ncol(regvars)+1 when condition has 2 levels
  coefmat <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(regvars)+2) 
  set.seed(8675309);
  if (parallelThreads > 1) {
    
    if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
      set.seed(8675309)
      cat('[bootstrap-PARALLEL] Working on ORDINARY NONPARAMETRIC BOOTSTRAP regression with ', parallelThreads, ' threads over ', nrow(cleanDat), ' iterations.\n Estimated time to complete:', round(120/parallelThreads*nrow(cleanDat)/2736,1), ' minutes.\n') #intermediate progress printouts would not be visible in parallel mode
      coefmat <- foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% {
        options(stringsAsFactors=FALSE)
        library(boot)
        thisexp <- as.numeric(cleanDat[i,])
        bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                           R=numboot, formula=thisexp~condition.AD +age+sex+PMI) #condition.ad+condition.ftdu+condition.tau+condition.pdd+condition.msa+condition.pd+condition.als
        ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median) 
        #b4_154_52 
        
        bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
        for (n in 1:ncol(bs.results$t)) {
          bs.stats[n] <- median(na.omit(bs.results$t[,n]))
        }
        bs.stats
        #cat('[bootstrap] Done for Protein ',i,'\n') #will not be visible
      }
      ## - residuals.MArrayLM(bs.stats, cleanDat[1,] )
      normExpr.reg <-foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% { (cleanDat[i,]- coefmat[i,3]*regvars[,"age"] - coefmat[i,4]*regvars[,"sex"]- coefmat[i,5]*regvars[,"PMI"]) }
      #Expr.Resid.Condition <- foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% { (cleanDat[i,]- (coefmat[i,1] + coefmat[i,3]*regvars[,"age"] + coefmat[i,4]*regvars[,"sex"] + coefmat[i,5]*regvars[,"PMI"]) )} 
      #row.names(Expr.Resid.Condition) <- row.names( cleanDat )
      #row.names(normExpr.reg) <- row.names( cleanDat )
      #Expr.resid <- cleanDat - normExpr.reg[,colnames(cleanDat)]
      #plot( as.factor(numericMetaSimple[colnames(cleanDat),]$AD), cleanDat['GAPDH|P04406',]  )
      #plot( as.factor(numericMetaSimple[colnames(cleanDat),]$AD), Expr.Resid.Condition['GAPDH|P04406',]  )
      
      #plot( as.factor(numericMetaSimple$AD), Expr.resid['APOE|P02649',]  )
      #plot( as.factor(numericMetaSimple$AD), Expr.resid['RPL36AL|Q969Q0',]  )
      #
      #wilcox.test( normExpr.reg[ 'GAPDH|P04406', ][row.names(numericMetaSimple[ numericMetaSimple$AD == 0, ])], 
      #             normExpr.reg[ 'GAPDH|P04406',][row.names(numericMetaSimple[ numericMetaSimple$AD == 1, ])]
      #)
      #wilcox.test( scale(DescTools::Winsorize(normExpr.reg[ 'GAPDH|P04406', ]))[row.names(numericMetaSimple[ numericMetaSimple$AD == 0, ]),], 
      #             scale(DescTools::Winsorize(normExpr.reg[ 'GAPDH|P04406',]))[row.names(numericMetaSimple[ numericMetaSimple$AD == 1, ]),]
      #)
      
      
    } else { 
      #linear model regression; faster but incomplete regression of Age, Sex, PMI effects, SO NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
      coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
      for (i in 1:nrow(cleanDat)) {
        if (i%%1000 == 0) {print(i)}
        lmmod1 <- lm(as.numeric(cleanDat[i,])~condition.ad+age+sex+PMI,data=regvars) #ALL 3 regression
        #lmmod1 <- lm(as.numeric(cleanDat[i,])~condition +age+sex+PMI,data=regvars)
        ##datpred <- predict(object=lmmod1,newdata=regvars)
        #_# hist(residuals.MArrayLM(lmmod1, Dat_Sink[1,]),breaks=60)
        coef <- coef(lmmod1)
        coefmat[i,] <- coef
        ## The full data - the undesired covariates
        normExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals 
        ## Also equivalent to <- thisexp - coef*var expression above
        cat('Done for Protein ',i,'\n')
      }
      #Fill in later iff needed
      Expr.resid <- NULL
    }
    
  } else {
    #single thread not handled
  }
  
  quantile( cleanDat.unreg[,1], c(0,0.025,0.25,0.5,0.75,0.975,1), na.rm=TRUE)
  quantile( normExpr.reg[,1], c(0,0.025,0.25,0.5,0.75,0.975,1), na.rm=TRUE)
  row.names(normExpr.reg) <- row.names(cleanDatImputed)
  
  cleanDat<-normExpr.reg[, row.names(numericMetaSimple)]
  
  Dat_Sink <- cleanDat
  Meta_Sink <- numericMetaSimple
  
  numericMetaSimple <- numericMetaSimple[ (numericMetaSimple$SampleForCombat %in% c("AD", "CT")) == T , ]
  numericMetaSimple <- numericMetaSimple[ colnames(cleanDat), ]
  numericMetaSimple <- numericMetaSimple[ complete.cases(numericMetaSimple$Gender), ]
  cleanDat <- cleanDat[, row.names(numericMetaSimple)]
  
  # numericMetaSimple <- numericMetaSimple
  varPart.combat <- fitExtractVarPartModel( cleanDat, form, numericMetaSimple)
  # sort variables (i.e. columns) by median fraction of variance explained
  vp.combat <- sortCols( varPart.combat )
  
  # Bar plot of variance fractions for the first 10 genes 
  plotPercentBars( vp.combat[1:10,] )
  # violin plot of contribution of each variable to total variance 
  plotVarPart( vp.combat )
  
  #DescTools::Winsorize
  Expr.resid <-foreach(i=1:nrow(cleanDat), .combine=rbind) %dopar% { scale(DescTools::Winsorize(cleanDat[i,]))[,1] }
  row.names(Expr.resid) <- row.names(Dat_Sink)
  
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
  return( list( AdjEXP = Dat_Sink, ScalWins =Expr.resid, Cov = Meta_Sink, NuMeta = numericMeta))
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
}

Anov <- function( Dat_S, Meta_S, NuMeta ){
  ##Applies ANOVA Testing
  #'@Dat_Sink Exoression data from impute
  #'@numericSimple meta data from impute
  ###########################
  ####                   #### 
  ####      ANOVA        #### 
  ####                   #### 
  ###########################
  #numericMeta <- numericMeta[ complete.cases( numericMeta$Age), ]
  cleanDat <- Dat_S[ ,row.names(Meta_S)]
  numericMetaSimple <- Meta_S
  numericMeta <- NuMeta
  
  metdat <- data.frame(Age=numericMeta[,"Age"],Gender=numericMeta[,"Gender"]) #,LED=numericMeta[,"LED"],Amyloid=numericMeta[,"Amyloid"],
  Grouping <- as.character(numericMeta$SampleForCombat)
  data = as.data.frame(cbind(colnames(cleanDat), Grouping, t(cleanDat)))
  colnames(data)[1:2]<-c("CODE","SampleType")
  #test run gets column headers for output
  i=3
  aov <- aov(data[,i]~SampleType, data=data)
  anovaresult <- anova(aov)
  tuk <- TukeyHSD(aov)
  #ASSUMES NO PAIRWISE COMPARISONS ARE MISSING FOR FIRST cleanDat protein (ROWS OF THIS data frame)--this is the template for comparison columns in ANOVAout
  tukresult1 <- data.frame(tuk$SampleType) 
  j = length(rownames(tukresult1))
  comparisonList <- rownames(tukresult1)
  
  line = c(paste("Protein", "F-Value", "Pr(>F)", sep=","))
  for(a in 1:length(comparisonList)) {
    line = c(paste(line,comparisonList[a],sep=","))
  }
  for (a in 1:length(comparisonList)) {
    line=c(paste(line,paste0("diff ",comparisonList[a]), 
                 sep=",")
    )
  }
  for (a in 1:length(comparisonList)) {
    line=c(paste(line, paste0("lwr ",comparisonList[a]),
                 sep=",")
    )
  }
  for (a in 1:length(comparisonList)) {
    line=c(paste(line, paste0("upr ",comparisonList[a]), sep=",")
    )
  }
  
  ## Fast code with apply
  dataclipped <- data[, 3:ncol(data)] # as.matrix(as.numeric(data[,3:ncol(data)]))
  SampleType <- data$SampleType
  ANOVAout <- apply(dataclipped, 2, function(x) {
    #as.double(x) instead of x corrects:   Error in lm.fit(x, y,... NA/NaN/Inf in 'y'
    x <- data.frame(x = as.double(x), SampleType = SampleType) 
    aov <- aov(x ~ SampleType, data = x)
    anovaresult <- anova(aov)
    tuk <- TukeyHSD(aov)
    tukresult <- data.frame(tuk$SampleType)
    if (length(rownames(tukresult)) == length(rownames(tukresult1))) {
      lwr <- as.vector( tukresult[,"lwr"] )
      upr <- as.vector( tukresult[,"upr"] )
      
      c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]), lwr, upr)
    } else {
      tukresult <- tukresult[match(rownames(tukresult1), rownames(tukresult)), ]
      lwr <- as.vector( tukresult[,"lwr"] )
      upr <- as.vector( tukresult[,"upr"] )
      c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]), lwr, upr)
    }
  })
  ANOVAout <- t(ANOVAout)
  ANOVAcols <- as.vector(data.frame(do.call("rbind", strsplit(as.character(line), "[,]"))))
  colnames(ANOVAout) <- ANOVAcols[2:length(ANOVAcols)]
  ANOVAout <- as.data.frame(ANOVAout)
  
  #T-test p value equivalency for single comparison, add in BH corr p value.
  if(length(unique(Grouping))==2) {
    #get BH FDR for ANOVA/T-test p values
    ANOVAout[,dim(ANOVAout)[2]+1] <- p.adjust(ANOVAout[,2],method="BH",n=dim(ANOVAout)[1]) 
    colnames(ANOVAout)[ dim(ANOVAout)[2] ] <- "FDR (BH)"
    ANOVAout[,dim(ANOVAout)[2]+1] <- p.adjust(ANOVAout[,2],method="fdr",n=dim(ANOVAout)[1]) 
    colnames(ANOVAout)[ dim(ANOVAout)[2] ] <- "FDR (FDR)"
  }
  
  ANOVAout <- cbind( ANOVAout, 
                     'CT-AD_Cor_BF' = p.adjust(ANOVAout$`CT-AD`, method = "bonferroni", n=dim(ANOVAout)[1]),
                     'CT-AD_Cor_BH' = p.adjust(ANOVAout$`CT-AD`, method = "BH", n=dim(ANOVAout)[1]),
                     'CT-AD_Cor_FDR' = p.adjust(ANOVAout$`CT-AD`, method = "fdr", n=dim(ANOVAout)[1])
  )
  
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
  return( list( AdjEXP = cleanDat, Cov = numericMeta, Stats = ANOVAout))
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
}

Test <- function( EXP, Pheno, Template, region ){
  #__Tests Cases V Controls T-Test and Wilcoxon Tests are applied
  #'@EXP The gene expression data frame eg Temp
  #'@Pheno The phenotype file eg MSBB
  #'@Template The Gene Template File eg Master
  #'@region the brain region eg "DLPFC"
  ### eg USAGE: Test( Temp, MSBB, Master )
  
  Template <- Template[ row.names(EXP), ]
  
  Pheno <- Pheno[ (Pheno$Control == 1) | (Pheno$AD == 1), ]
  EXP <- EXP[ , row.names(Pheno) ]
  Pheno <- cbind(Pheno, REGION = region)
  for( Region in names(table(Pheno$REGION)) ){
    
    #Build Target Location for Statistics
    Template <- cbind(Template, NA, NA, NA, NA, NA )
    Names <- c('Log2FC', 'TtestPVal', 'Cor_TTest_BF', 'Wilcox_PVal', 'Cor_Wilcox_BF')
    ##Targeted Columns in final data frame for this Region
    targets <- c( First = dim(Template)[2]-4, Last = dim(Template)[2])
    colnames(Template)[ targets["First"]:targets["Last"] ] <- c( paste0( paste0(Region, "_"), Names ))
    
    ##Filter for region of interest
    Phen <- Pheno[ Pheno$REGION == Region, ]
    if( dim(Phen)[1] == dim(EXP)[2]){
      Ex <- EXP[ , colnames(EXP) %in% row.names(Phen) ]
      
      ##Number of Corrections for this dataset:
      Corrections <- dim( EXP[ complete.cases(EXP), ] )[1]
      
      #Asess if the Gene/Protein is Quantified in that dataset
      for( i in 1:dim(EXP)[1]){
        if( length(na.omit(EXP[ i ,])) == 0 ){
        }else{
          CNTRLs <- Ex[ i , row.names(Phen[ Phen$Control == 1, ]) ]
          CASEs <- Ex[ i , row.names(Phen[ Phen$AD == 1, ]) ]
          
          OUT <- c( rep(NA,5))
          OUT[1] <- median(CASEs) - median(CNTRLs)
          #Ttest:
          OUT[2] <- t.test( CNTRLs, CASEs )$p.value
          if( (Corrections * OUT[2]) > 1 ){
            OUT[3] <- 1
          }else{
            OUT[3] <- Corrections * OUT[2]
          }
          
          #Wilcox
          OUT[4] <- wilcox.test( CNTRLs, CASEs  )$p.value
          if( (Corrections * OUT[4]) > 1 ){
            OUT[5] <- 1
          }else{
            OUT[5] <- Corrections * OUT[4]
          }
          
          Template[ i, targets['First']:targets['Last'] ] <- OUT
        }
      }
    }else{
      print( paste0( "WARNING: ", Region, " Samples are either not in the expression or the sample numbers don't match up!") )
      print( paste0( "Skiping: ", Region, " Region") )
    }
    
  }
  Template <- cbind( Template[, 1:7], 
                     p.adjust( Template[,6], method = "BH", n=dim(Template)[1] ), 
                     Template[, 8:9], 
                     p.adjust( Template[,8], method = "BH", n=dim(Template)[1] )
  )
  colnames(Template)[8] <- paste0(Region, "_Cor_TTest_BH")
  colnames(Template)[11] <- paste0(Region, "_Cor_Wilcox_BH")
  return(Template)
}

#Eq for Std Error
se <- function(x) sqrt(var(x)/length(x))

CalcEr <- function(EXP, Meta){
  Out <- matrix( 0, dim(EXP)[1], 7)
  row.names(Out) <- row.names(EXP)
  colnames(Out) <- c( "CT_X", "Ct_SE", "AD_X", "AD_SE", "FC_X", "FC_lwr", "FC_upr")
  
  for( gene in row.names(Out)){
    CT <- mean( EXP[ gene, row.names(Meta[ Meta$SampleForCombat == "CT",]) ] )
    AD <- mean( EXP[ gene, row.names(Meta[ Meta$SampleForCombat == "AD",]) ] )
    
    CT_se <- se( EXP[ gene, row.names(Meta[ Meta$SampleForCombat == "CT",]) ] )
    AD_se <-  se( EXP[ gene, row.names(Meta[ Meta$SampleForCombat == "AD",]) ] )
    
    #Calc mean log2 Fold Change
    FC <- AD - CT
    LWR <- (AD-1.96*AD_se) - (CT-1.96*CT_se)
    UPR <- (AD+1.96*AD_se) - (CT+1.96*CT_se)
    if( LWR < UPR ){
      Out[ gene, ] <- c( CT, CT_se, AD, AD_se, FC, LWR, UPR)
    }else{
      Out[ gene, ] <- c( CT, CT_se, AD, AD_se, FC, UPR, LWR)
    }
  }
  return(as.data.frame(Out))
}

Filler <- function( Data, Template){
  #'Usage: Filler( AG_Banner, FinMast)
  #'@Data the INDV data frame to be integrated
  #'@Template the template to integrate into
  
  if( length(table( Data$Tissue )) > 1 ){
    print( "WARNING Multiple Tissue Types in Data" )
  }else{
    Template$Tissue <- names(table( Data$Tissue ))
  }
  
  for( name in row.names(Data) ){
    Template[ name, ]$Log2_FC <- Data[ name, ]$Log2_FC
    Template[ name, ]$CI_Upr <- Data[ name, ]$CI_Upr
    Template[ name, ]$CI_Lwr <- Data[ name, ]$CI_Lwr
    Template[ name, ]$PVal <- Data[ name, ]$PVal
    Template[ name, ]$FDR_Cor_PVal <- Data[ name, ]$FDR_Cor_PVal
  }
  Template <- Template[ complete.cases(Template$ENSG),]
  return( Template )
}

##############################
######    LOAD DATA     ######
##############################
#Load AGORA Template
Master <- syn_temp$get('syn18914952')
Master <- read.csv(Master$path, row.names = 1 )

####Load Proteomics Data:
BLSA <- syn_temp$get('syn18914920')
BLSA <- read.csv(BLSA$path)

Mayo_dat <- syn_temp$get('syn9637748')
Mayo_dat <- read.table(Mayo_dat$path, sep = '\t', header = T )

MSBB_dat <- syn_temp$get('syn6100414')
MSBB_dat <- read.table(MSBB_dat$path, sep = '\t', header = T )
##Take out second reps from MSBB
MSBB_dat <- MSBB_dat[ , colnames(MSBB_dat)[ grepl('b7r2', colnames(MSBB_dat) ) == F ] ]

BANNER <- syn_temp$get('syn18918360')
BANNER <- as.matrix(read.csv( file = BANNER$path, header=T, row.names=1))

### - Process
#MSBB_CDat <- Process( MSBB_dat, 0 )
#Mayo_CDat <- Process( Mayo_dat, 0 )

MSBB_CDat <- Process( MSBB_dat, 1 )
Mayo_CDat <- Process( Mayo_dat, 1 )
BLSA_CDat <- Process( BLSA, 1 )
Banner_CDat <- BANNER#Process( BANNER, 0 )

### - Impute
BLSA_RDY <- Logger_KNN( BLSA_CDat, "BLSA" )
Mayo_RDY <- Logger_KNN( Mayo_CDat, "Mayo" )
MSBB_RDY <- Logger_KNN( MSBB_CDat, "MSBB" )
exprMat.VAL <- Logger_KNN(Banner_CDat, "Banner")

##########Load Phenotype Meta Data
#BANNER
Combat_Info <- syn_temp$get('syn18914620')
Combat_Info <- read.csv( Combat_Info$path, header=T)
rownames(Combat_Info) <- Combat_Info$Sample

Combat_Info <- cbind( Combat_Info, Diagnosis = 0 )
Combat_Info[ Combat_Info$AD == 1 , ]$Diagnosis <- 'AD'
Combat_Info[ Combat_Info$Control == 1 , ]$Diagnosis <- 'CT'
Combat_Info[ grepl( 'bgis', Combat_Info$Sample) == T , ]$Diagnosis <- 'BGIS'
Combat_Info[ grepl( 'egis', Combat_Info$Sample) == T , ]$Diagnosis <- 'EGIS'

##BLSA
BLSA_Met <- syn_temp$get('syn18914694')
BLSA_Met <- read.csv(file = BLSA_Met$path, row.names = 1)
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
BLSA_Met <- BLSA_Met[ (row.names(BLSA_Met) %in% colnames(BLSA_RDY))==T, ]
complete.cases(BLSA_Met)

BLSA_data <- syn_temp$get('syn18918327')
BLSA_data <- read_excel(BLSA_data$path, sheet = 2)

BLSA_Fdata <- as.data.frame( cbind( 'Sample' = BLSA_data$SanpleName, 'SampleName' = BLSA_data$`Sample#`, 'BRC' = BLSA_data$`BRC#`, 'MFG' = BLSA_data$MFG, 'Precuneus' = BLSA_data$Precuneus, 'Control' = BLSA_data$CT, 'AsymAD' = BLSA_data$AsymAD, 'AD' = BLSA_data$AD, 'CERAD' = BLSA_data$CERAD, 'Braak' = BLSA_data$BRAAK, 'Age' = BLSA_data$AGE, 'Gender' = BLSA_data$SEX,  'PMI' = BLSA_data$PMI, 'ApoE' = BLSA_data$`ApoE Numeric Score` ) )
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
BLSA_Fdata <- BLSA_Fdata[ colnames(BLSA_RDY), ]
for( i in 2:15){
  BLSA_Fdata[,i] <- as.numeric( BLSA_Fdata[,i] )
}

###MAYO
Mayo_Met <- syn_temp$get('syn18914935')
Mayo_Met <- read.csv( file = Mayo_Met$path, row.names = 1)
row.names(Mayo_Met) <- Mayo_Met$SampleNames

Mayo_Met <- Mayo_Met[match(colnames(Mayo_RDY),rownames(Mayo_Met)),]
Mayo_Met <- Mayo_Met[ , c("SampleNames", "batch", "Diagnosis", "Gender", "Braak", "AgeAtDeath", "PMI" ) ]
Mayo_Met$batch <- gsub( "b", "", Mayo_Met$batch )
colnames( Mayo_Met ) <- c( "Sample", "Batch", "Diagnosis", "Gender", "Braak", "Age", "PMI" )
Mayo_Met <- cbind( Mayo_Met, AD = 0, Control = 0, PSP =0)
Mayo_Met[ (Mayo_Met$Diagnosis == 'AD') == T , ]$AD <- 1
Mayo_Met[ (Mayo_Met$Diagnosis == 'Control') == T , ]$Diagnosis <- 'CT'
Mayo_Met[ (Mayo_Met$Diagnosis == 'CT') == T , ]$Control <- 1
Mayo_Met[ (Mayo_Met$Diagnosis == 'PSP') == T , ]$PSP <- 1
Mayo_Met<- cbind(Mayo_Met, SampleForCombat=Mayo_Met$Diagnosis)
Mayo_Met <- Mayo_Met[ colnames(Mayo_RDY), ]
Mayo_Met[ grepl( 'egis', Mayo_Met$Sample) == T , ]$Diagnosis <- 'EGIS'
Mayo_Met[ grepl( 'mgis', Mayo_Met$Sample) == T , ]$Diagnosis <- 'MGIS'

###MSBB
MSBB_Pheno <- syn_temp$get('syn18914936')
MSBB_Pheno <- read.csv(file = MSBB_Pheno$path, header = T, sep = "\t")

#Get rid of the repeated sample runs
Repeats <- MSBB_Pheno[ MSBB_Pheno$Batch == 8, ]$RunName
MSBB_RDY <- MSBB_RDY[ , (colnames(MSBB_RDY) %in% Repeats) == F ]
MSBB_Pheno <- MSBB_Pheno[ (MSBB_Pheno$RunName %in% Repeats) == F , ]
row.names(MSBB_Pheno) <- MSBB_Pheno$RunName

##Add GIS Samples to the Pheno File
Adds <- colnames(MSBB_RDY)[ (colnames(MSBB_RDY) %in% MSBB_Pheno$RunName) == F ]
Addon <- as.data.frame( matrix( NA, length(Adds), dim(MSBB_Pheno)[2] ) )
row.names( Addon ) <- Adds 
colnames( Addon ) <- colnames(MSBB_Pheno)
Addon$RunName <- row.names( Addon )
Addon$CT <- 0
Addon$AFF <- 0

##Need to add on batch
SplitOut <- do.call( rbind, strsplit(Addon$RunName, '_') ) 
Addon$Batch <- gsub( "r2", "", gsub( "b", "", SplitOut[,1] ) )
#here...
MSBB_Pheno <- rbind( MSBB_Pheno, Addon )
MSBB_Pheno <- MSBB_Pheno[ colnames(MSBB_RDY), ]

colnames(MSBB_Pheno)[2] <- "Sample"
colnames(MSBB_Pheno)[7] <- "Braak"
colnames(MSBB_Pheno)[4] <- "Control"
colnames(MSBB_Pheno)[5] <- "AD"
colnames(MSBB_Pheno)[8] <- "Age"
colnames(MSBB_Pheno)[9] <- "Gender"
MSBB_Pheno <- cbind( MSBB_Pheno, Diagnosis = 0 )
MSBB_Pheno[ MSBB_Pheno$AD == 1 , ]$Diagnosis <- 'AD'
MSBB_Pheno[ MSBB_Pheno$Control == 1 , ]$Diagnosis <- 'CT'
MSBB_Pheno[ grepl( "bmgis", MSBB_Pheno$Sample) == T, ]$Diagnosis <- 'BGIS'

tPhen <- syn_temp$get('syn18914939')
tPhen <- read.csv(file=tPhen$path)
row.names(tPhen) <- tPhen$RunName
MSBB_Pheno$Age <- tPhen[ row.names(MSBB_Pheno), ]$AOD

### Batch Correct and Regress
Imputed <- Runner( exprMat.VAL, Combat_Info, "BANNER")

#BLSA:
BLSA_Fdata <- cbind(BLSA_Fdata, Batch = 1 )
BLSA_RunPData <- BLSA_Fdata[ , c("Sample", "Gender", "Age", "Control", "AD", "AsymAD", "CERAD","PMI", "Braak", "ApoE", "Batch", "Diagnosis"  )]
#Code Sex as 2 and 1 not 1 and 0
BLSA_RunPData$Gender <- BLSA_RunPData$Gender+1
Imputed_BLSA <- Runner( BLSA_RDY, BLSA_RunPData, "BLSA")

#Mayo:
colnames(Mayo_Met)[11] <- "Diagnosis"
Imputed_Mayo <- Runner( Mayo_RDY, Mayo_Met, "MAYO")

#MSBB:
#colnames(Combat_Info)[ colnames(Combat_Info) %in% colnames(MSBB_Pheno)]
MSBB_Pheno$REGION <- "AntPFC"
##_##  MSBB_Pheno <- MSBB_Pheno[ colnames(MSBB_CDat),]
MSBB_Pheno <- MSBB_Pheno[ row.names(MSBB_Pheno)[ row.names(MSBB_Pheno)  %in% colnames(MSBB_CDat)], ]
head(as.matrix(as.data.frame(MSBB_RDY))[1:10])
temp<-apply(MSBB_RDY,2,as.numeric) 
row.names(temp) <- row.names(MSBB_RDY) 
MSBB_RDY <- temp
Imputed_MSBB <- Runner( MSBB_RDY, MSBB_Pheno, "MSBB")
#Imputed_MSBB <- Runner( MSBB_RDY[,row.names(MSBB_Pheno)], MSBB_Pheno, "MSBB")

###Run Statistical Testing
#Anova
MAYO_Anov <- Anov( Imputed_Mayo$AdjEXP, Imputed_Mayo$Cov, Imputed_Mayo$NuMeta)
Banner_Anov <- Anov( Imputed$AdjEXP, Imputed$Cov, Imputed$NuMeta)
BLSA_Anov <- Anov( Imputed_BLSA$AdjEXP, Imputed_BLSA$Cov, Imputed_BLSA$NuMeta)
MSBB_Anov <- Anov( Imputed_MSBB$AdjEXP, Imputed_MSBB$Cov, Imputed_MSBB$NuMeta)
#tTest and Wilcoxon
BannerStats <- Test( Imputed$AdjEXP, Imputed$Cov, Master, "DLPFC")
MAYOStats <- Test( Imputed_Mayo$AdjEXP, Imputed_Mayo$Cov, Master, "TC")
BLSAStats <- Test( Imputed_BLSA$AdjEXP, Imputed_BLSA$Cov, Master, "MFG")
MSBBStats <- Test( Imputed_MSBB$AdjEXP, Imputed_MSBB$Cov, Master, "AntPFC")

Stats <- list( Banner = list(Anova = Banner_Anov$Stats, Regular = BannerStats), 
               Mayo = list(Anova = MAYO_Anov$Stats, Regular = MAYOStats),
               BLSA = list(Anova = BLSA_Anov$Stats, Regular = BLSAStats),
               MSBB = list(Anova = MSBB_Anov$Stats, Regular = MSBBStats)
)


##Final Agora Data Frame Check Stats:
Banner_FC <- CalcEr( Imputed$AdjEXP, Imputed$NuMeta)
Mayo_FC <- CalcEr( Imputed_Mayo$AdjEXP, Imputed_Mayo$NuMeta)
BLSA_FC <- CalcEr( Imputed_BLSA$AdjEXP, Imputed_BLSA$NuMeta)
MSBB_FC <- CalcEr( Imputed_MSBB$AdjEXP, Imputed_MSBB$NuMeta)


Stats <- list( Banner = list(Anova = Banner_Anov$Stats, Regular = BannerStats, RoughFC = Banner_FC ), 
               Mayo = list(Anova = MAYO_Anov$Stats, Regular = MAYOStats, RoughFC = Mayo_FC ),
               BLSA = list(Anova = BLSA_Anov$Stats, Regular = BLSAStats, RoughFC = BLSA_FC),
               MSBB = list(Anova = MSBB_Anov$Stats, Regular = MSBBStats, RoughFC = MSBB_FC )
)

###Final Agora Data Frame - Make Frame
AG_Banner <- cbind( Stats$Banner$Regular[ , 1:4],
                    Tissue = "DLPFC",
                    Log2_FC = -1*Stats$Banner$Anova$`diff CT-AD`,
                    CI_Upr = -1*Stats$Banner$Anova$`lwr CT-AD`,
                    CI_Lwr = -1*Stats$Banner$Anova$`upr CT-AD`,
                    PVal = Stats$Banner$Anova$`CT-AD`,
                    FDR_Cor_PVal = Stats$Banner$Anova$`CT-AD_Cor_FDR`
)

AG_BLSA <- cbind( Stats$BLSA$Regular[ , 1:4],
                  Tissue = "MFG",
                  Log2_FC = -1*Stats$BLSA$Anova$`diff CT-AD`,
                  CI_Upr = -1*Stats$BLSA$Anova$`lwr CT-AD`,
                  CI_Lwr = -1*Stats$BLSA$Anova$`upr CT-AD`,
                  PVal = Stats$BLSA$Anova$`CT-AD`,
                  FDR_Cor_PVal = Stats$BLSA$Anova$`CT-AD_Cor_FDR`
)

AG_Mayo <- cbind( Stats$Mayo$Regular[ , 1:4],
                  Tissue = "TCX",
                  Log2_FC = -1*Stats$Mayo$Anova$`diff CT-AD`,
                  CI_Upr = -1*Stats$Mayo$Anova$`lwr CT-AD`,
                  CI_Lwr = -1*Stats$Mayo$Anova$`upr CT-AD`,
                  PVal = Stats$Mayo$Anova$`CT-AD`,
                  FDR_Cor_PVal = Stats$Mayo$Anova$`CT-AD_Cor_FDR`
)

AG_MSBB <- cbind( Stats$MSBB$Regular[ , 1:4 ],
                  Tissue = rep.int( "AntPFC", dim(Stats$MSBB$Regular)[1] ),
                  Log2_FC = -1*Stats$MSBB$Anova$`diff CT-AD`,
                  CI_Upr = -1*Stats$MSBB$Anova$`lwr CT-AD`,
                  CI_Lwr = -1*Stats$MSBB$Anova$`upr CT-AD`,
                  PVal = Stats$MSBB$Anova$`CT-AD`,
                  FDR_Cor_PVal = Stats$MSBB$Anova$`CT-AD_Cor_FDR`
)

FinMast <- Master[ (row.names(Master) %in% c( row.names(AG_MSBB), row.names(AG_Banner), row.names(AG_BLSA), row.names(AG_Mayo) )) == T, ]
FinMast <- cbind(FinMast, Tissue = NA, Log2_FC=NA, CI_Upr=NA, CI_Lwr=NA, PVal=NA, FDR_Cor_PVal=NA)

Banner_Mast <- Filler( AG_Banner, FinMast)
BLSA_Mast<- Filler( AG_BLSA, FinMast)
Mayo_Mast <- Filler( AG_Mayo, FinMast)
MSBB_Mast <- Filler( AG_MSBB, FinMast)

Final_Agora <- rbind( Banner_Mast, BLSA_Mast, Mayo_Mast, MSBB_Mast  )

WD <- getwd()
colnames(Final_Agora)[ 10 ] <- "Cor_PVal"
#_# write.csv( Final_Agora, paste0("/Users/jgockley/Desktop/Projects/AMP-AD/Proteomics/Compare_To_Emory/", "/Final_Agora_DE.csv") , row.names = F)
#_# saveRDS(Stats, paste0(WD, "Final_DE_Stats.RDS") )   

#############################################################################################
#############################################################################################
# - Meta Analysis - Data Clean

#Mayo - TCX
#Imputed_Mayo$ScalWins

#BLSA - MFG
#Imputed_BLSA$ScalWins

#Banner - DLPFC
#Imputed$ScalWins

#MSBB - AntPFC
#Imputed_MSBB$ScalWins

Cols <- dim(Imputed_Mayo$ScalWins)[2]+dim(Imputed_BLSA$ScalWins)[2]+dim(Imputed$ScalWins)[2]+dim(Imputed_MSBB$ScalWins)[2]

COLS <- c( colnames(Imputed_Mayo$ScalWins),
           colnames(Imputed_BLSA$ScalWins),
           colnames(Imputed$ScalWins),
           colnames(Imputed_MSBB$ScalWins))

ROWS <- c(row.names(Imputed_Mayo$ScalWins),
  row.names(Imputed_BLSA$ScalWins),
  row.names(Imputed$ScalWins),
  row.names(Imputed_MSBB$ScalWins))
ROWS <- ROWS[ !duplicated(ROWS) ]

expr <- as.data.frame(matrix(NA, length(ROWS), Cols ))

colnames(expr) <- COLS
row.names(expr) <- ROWS

expr[ row.names(Imputed_Mayo$ScalWins),colnames(Imputed_Mayo$ScalWins) ] <- Imputed_Mayo$ScalWins[ row.names(Imputed_Mayo$ScalWins),colnames(Imputed_Mayo$ScalWins) ] 
expr[ row.names(Imputed_BLSA$ScalWins),colnames(Imputed_BLSA$ScalWins) ] <- Imputed_BLSA$ScalWins[ row.names(Imputed_BLSA$ScalWins),colnames(Imputed_BLSA$ScalWins) ] 
expr[ row.names(Imputed$ScalWins),colnames(Imputed$ScalWins) ] <- Imputed$ScalWins[ row.names(Imputed$ScalWins),colnames(Imputed$ScalWins) ] 
expr[ row.names(Imputed_MSBB$ScalWins),colnames(Imputed_MSBB$ScalWins) ] <- Imputed_MSBB$ScalWins[ row.names(Imputed_MSBB$ScalWins),colnames(Imputed_MSBB$ScalWins) ] 


####Get Meta Data
Mayo_Met <- syn_temp$get('syn18914935')
Mayo_Met <- read.csv( file = Mayo_Met$path, row.names = 1)
row.names(Mayo_Met) <- Mayo_Met$SampleNames

Mayo_Met <- Mayo_Met[match(colnames(Mayo_RDY),rownames(Mayo_Met)),]
Mayo_Met$batch <- gsub( "b", "", Mayo_Met$batch )
Mayo_Met <- Mayo_Met[ colnames(Mayo_RDY), ]
Mayo_Met$Diagnosis <- as.character(Mayo_Met$Diagnosis)
Mayo_Met[ grepl( 'egis', Mayo_Met$SampleNames) == T , ]$Diagnosis <- 'EGIS'
Mayo_Met[ grepl( 'mgis', Mayo_Met$SampleNames) == T , ]$Diagnosis <- 'MGIS'

#Format for Metanalysis:
Mayo_Met <- Mayo_Met[ Mayo_Met$Diagnosis %in% c('AD','Control'), ]
Mayo_Met$Study<-'Mayo'
colnames(Mayo_Met)[ colnames(Mayo_Met) == 'SampleNames' ] <- 'SampleID'
#APOE4
Mayo_Met$APOE4 <-  0 
Mayo_Met[ Mayo_Met$ApoE %in% c('e24', 'e34'), ]$APOE4 <- 1
Mayo_Met[ Mayo_Met$ApoE %in% 'e44', ]$APOE4 <- 2
Mayo_Met$Tissue <- 'TCX'

#SEX
Mayo_Met[ Mayo_Met$Gender == 1,]$Gender <- 'MALE'
Mayo_Met[ Mayo_Met$Gender == 2,]$Gender <- 'FEMALE'

#Diagnosis
Mayo_Met[ Mayo_Met$Diagnosis == 'Control',]$Diagnosis <- 'CONTROL'
Mayo_Met[ Mayo_Met$Diagnosis == 'AD',]$Diagnosis <- 'AD'
Mayo_Met <- Mayo_Met[ , c('Study', 'SampleID', 'Gender', 'Diagnosis', 'APOE4', 'Tissue' ) ]
  #head(Mayo_Met[ , c('Study', 'SampleID', 'Gender', 'Diagnosis', 'APOE4', 'Tissue' ) ])


##BLSA
BLSA_Met <- syn_temp$get('syn18914694')
BLSA_Met <- read.csv(file = BLSA_Met$path, row.names = 1)
BLSA_Met$Tissue <- 'MFG'
BLSA_Met$Study <- 'BLSA'
colnames(BLSA_Met)[ colnames(BLSA_Met) =='ApoE'] <- 'APOE4'
colnames(BLSA_Met)[ colnames(BLSA_Met) =='SEX'] <- 'Gender'

#SEX
BLSA_Met[ BLSA_Met$Gender == 0,]$Gender <- 'MALE'
BLSA_Met[ BLSA_Met$Gender == 1,]$Gender <- 'FEMALE'

#Diagnosis
BLSA_Met[BLSA_Met$Affected==0,]$Affected <- 'CONTROL'
BLSA_Met[BLSA_Met$Affected==1,]$Affected <- 'ASYMP_AD'
BLSA_Met[BLSA_Met$Affected==2,]$Affected <- 'AD'

colnames(BLSA_Met)[ colnames(BLSA_Met) == 'SampleName' ] <- 'SampleID'
colnames(BLSA_Met)[ colnames(BLSA_Met) == 'Affected' ] <- 'Diagnosis'

#Filter
BLSA_Met <- BLSA_Met[ BLSA_Met$Diagnosis %in% c('CONTROL','AD'),]
BLSA_Met <- BLSA_Met[ , c('Study', 'SampleID', 'Gender', 'Diagnosis', 'APOE4', 'Tissue' ) ]
#head(BLSA_Met[ , c('Study', 'SampleID', 'Gender', 'Diagnosis', 'APOE4', 'Tissue' ) ])

#BANNER
Combat_Info <- syn_temp$get('syn18914620')
Combat_Info <- read.csv( Combat_Info$path, header=T)
rownames(Combat_Info) <- Combat_Info$Sample

Combat_Info$Study <- 'BANNER'
Combat_Info$Tissue <- 'DLPFC'
colnames(Combat_Info)[ colnames(Combat_Info) == 'Sample' ] <- 'SampleID'
colnames(Combat_Info)[ colnames(Combat_Info) == 'ApoE' ] <- 'APOE4'

Combat_Info <- Combat_Info[ Combat_Info$Control == 1 | Combat_Info$AD == 1, ]
Combat_Info[ Combat_Info$Gender == 1,]$Gender <- 'MALE'
Combat_Info[ Combat_Info$Gender == 2,]$Gender <- 'FEMALE'

Combat_Info[ Combat_Info$APOE4 == -1 & (is.na(Combat_Info$APOE4)==F), ]$APOE4 <- 0
Combat_Info[ Combat_Info$APOE4 == -2 & (is.na(Combat_Info$APOE4)==F), ]$APOE4 <- 0

Combat_Info$Diagnosis <- NA
Combat_Info[ Combat_Info$Control == 1, ]$Diagnosis <- 'CONTROL'
Combat_Info[ Combat_Info$AD == 1, ]$Diagnosis <- 'AD'

Banner_Met <- Combat_Info[ , c('Study', 'SampleID', 'Gender', 'Diagnosis', 'APOE4', 'Tissue' ) ]
  #head(Banner_Met)

###MSBB
MSBB_Pheno <- syn_temp$get('syn18914936')
MSBB_Pheno <- read.csv(file = MSBB_Pheno$path, header = T, sep = "\t")

#Get rid of the repeated sample runs
Repeats <- MSBB_Pheno[ MSBB_Pheno$Batch == 8, ]$RunName
MSBB_RDY <- MSBB_RDY[ , (colnames(MSBB_RDY) %in% Repeats) == F ]
MSBB_Pheno <- MSBB_Pheno[ (MSBB_Pheno$RunName %in% Repeats) == F , ]
row.names(MSBB_Pheno) <- MSBB_Pheno$RunName

##Add GIS Samples to the Pheno File
Adds <- colnames(MSBB_RDY)[ (colnames(MSBB_RDY) %in% MSBB_Pheno$RunName) == F ]
Addon <- as.data.frame( matrix( NA, length(Adds), dim(MSBB_Pheno)[2] ) )
row.names( Addon ) <- Adds 
colnames( Addon ) <- colnames(MSBB_Pheno)
Addon$RunName <- row.names( Addon )
Addon$CT <- 0
Addon$AFF <- 0

##Need to add on batch
SplitOut <- do.call( rbind, strsplit(Addon$RunName, '_') ) 
Addon$Batch <- gsub( "r2", "", gsub( "b", "", SplitOut[,1] ) )
#here...
MSBB_Pheno <- rbind( MSBB_Pheno, Addon )
MSBB_Pheno <- MSBB_Pheno[ colnames(MSBB_RDY), ]

colnames(MSBB_Pheno)[2] <- "Sample"
colnames(MSBB_Pheno)[7] <- "Braak"
colnames(MSBB_Pheno)[4] <- "Control"
colnames(MSBB_Pheno)[5] <- "AD"
colnames(MSBB_Pheno)[8] <- "Age"
colnames(MSBB_Pheno)[9] <- "Gender"
MSBB_Pheno <- cbind( MSBB_Pheno, Diagnosis = 0 )
MSBB_Pheno[ MSBB_Pheno$AD == 1 , ]$Diagnosis <- 'AD'
MSBB_Pheno[ MSBB_Pheno$Control == 1 , ]$Diagnosis <- 'CT'
MSBB_Pheno[ grepl( "bmgis", MSBB_Pheno$Sample) == T, ]$Diagnosis <- 'BGIS'

tPhen <- syn_temp$get('syn18914939')
tPhen <- read.csv(file=tPhen$path)
row.names(tPhen) <- tPhen$RunName
MSBB_Pheno$Age <- tPhen[ row.names(MSBB_Pheno), ]$AOD

MSBB_Pheno <- MSBB_Pheno[ MSBB_Pheno$Diagnosis == 'AD' | MSBB_Pheno$Diagnosis == 'CT', ]
MSBB_Pheno[ MSBB_Pheno$Diagnosis == 'CT', ]$Diagnosis <- "CONTROL"

colnames(MSBB_Pheno)[ colnames(MSBB_Pheno) == 'Sample' ] <- 'SampleID'
colnames(MSBB_Pheno)[ colnames(MSBB_Pheno) == 'REGION' ] <- 'Tissue'

MSBB_Pheno[ MSBB_Pheno$Gender == 1, ]$Gender <- 'MALE'
MSBB_Pheno[ MSBB_Pheno$Gender == 2, ]$Gender <- 'FEMALE'

MSBB_Pheno$APOEt1 <- NA
MSBB_Pheno$APOEt2 <- NA

MSBB_Pheno[ MSBB_Pheno$Apo1 < 4 & is.na(MSBB_Pheno$Apo1)==F, ]$APOEt1 <- 0
MSBB_Pheno[ MSBB_Pheno$Apo1 == 4 & is.na(MSBB_Pheno$Apo1)==F, ]$APOEt1 <- 1
MSBB_Pheno[ MSBB_Pheno$Apo2 < 4 & is.na(MSBB_Pheno$Apo2)==F, ]$APOEt2 <- 0
MSBB_Pheno[ MSBB_Pheno$Apo2 == 4 & is.na(MSBB_Pheno$Apo2)==F, ]$APOEt2 <- 1

MSBB_Pheno$APOE4 <- MSBB_Pheno$APOEt1 + MSBB_Pheno$APOEt2
MSBB_Pheno$Study <- 'MSBB'
MSBB_Pheno <- MSBB_Pheno[ , c('Study', 'SampleID', 'Gender', 'Diagnosis', 'APOE4', 'Tissue' ) ]
  #head(MSBB_Pheno)

covar <- rbind(MSBB_Pheno, Banner_Met, Mayo_Met, BLSA_Met)

Keeps <- as.character( covar$SampleID[ as.character( covar$SampleID ) %in% colnames(expr)] )
covar <- covar[ Keeps, ]
covar$SampleID <- as.character(covar$SampleID)
expr <- expr[ , Keeps ]

#head(expr[,1:4])
expr$peptide_id <- row.names(expr)
expr <- expr[ , c('peptide_id', colnames(expr)[ colnames(expr) != 'peptide_id']) ]


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# - Meta Analysis - Model

tissue.dx.summary = plyr::ddply(covar, .(Tissue, Diagnosis), .fun = function(x, y){
  data.frame(peptide_id = y$peptide_id,
             n = dim(x)[1],
             mn = rowMeans(y[,x$SampleID], na.rm = T),
             sd = apply(y[,x$SampleID], 1, sd, na.rm = T))
}, expr)
# Perform meta-analysis for AD-CONTROL comparison
#[(tissue.dx.summary$Tissue%in%'CBE')==F,]
library(meta)
meta.anlz.ad_cntrl = plyr::ddply(tissue.dx.summary[complete.cases(tissue.dx.summary),], .(peptide_id), .fun = function(x){
  exp.effect = dplyr::filter(x, Diagnosis == 'AD')
  rownames(exp.effect) = exp.effect$Tissue
  cntrl.effect = dplyr::filter(x, Diagnosis == 'CONTROL')
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

# Save DataFrame and push to Synapse
write.csv( meta.anlz.ad_cntrl, file = 'Raw_Proteomics_MetaAnalysis.csv')




