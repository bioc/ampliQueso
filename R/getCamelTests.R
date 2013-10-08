camelTest<-function(iBedFile, iCovdesc='covdesc',iT1, iT2, iNorm=c('none','density','minMax'),  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iSizes = NULL, iNPerm=100,iParallel=TRUE)
{
  if(iParallel == TRUE)
    .camelTestParReg(iBedFile, iCovdesc,iT1, iT2, iNorm,  iMeasure, iSizes = NULL, iNPerm)
  else
    .camelTest(iBedFile, iCovdesc,iT1, iT2, iNorm,  iMeasure, iSizes = NULL, iNPerm)
}



.camelTest <- function(iBedFile, iCovdesc='covdesc',iT1, iT2, iNorm=c('none','density','minMax'),  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iSizes = NULL, iNPerm=100)
{
  
  #mwiewior: reading BED
  annot<-read.table(iBedFile)
  numgenes <- dim(annot)[1]
  cvd <- read.table(iCovdesc)
  n.samples <- dim(cvd)[1]
  n1 <- which(cvd== iT1) 
  n2 <- which(cvd== iT2)
  
  #mwiewior:initiate output matrix
  #testResults <- matrix(nrow=numgenes, ncol= (length(iMeasure)* length(iNorm)) )
  testResults <- data.frame(t(rep(NA, (length(iMeasure)* length(iNorm)) ) ) )
  measureColNames <-c()
  for(n in iNorm){
    for(m in iMeasure){
      measureColNames<-cbind(measureColNames,paste(m,n,sep="_"))
    }
  }
  colnames(testResults)<-measureColNames
  testResults <- testResults[-1,]
  for (i in 1:numgenes)
  {
    chr <- as.character(annot[i,1]);
    st <- as.numeric(annot[i,2]);
    end <- as.numeric(annot[i,3]);
    regName <-annot[i,6];
    strand <- as.character(annot[i,4]);
    #print(strand)
    strand <- .strand2number(str)
    
    rs <- newSeqReads(chr,st,end,strand)
    rs <- getBamData(rs,1:n.samples, covdesc.file = iCovdesc)
    nd <- getCoverageFromRS(rs,1:n.samples)
    
    
    ###normalizacja
    if (!is.null(iSizes) )
      for (j in 1:n.samples)
        setSAXPYData(nd,iSizes[j,2],j) 
       # nd@data[[j]]<-nd@data[[j]]*iSizes[j,2] # tutaj jest normalizacja globalna wzgledem rozmiarow prob
    
    
    ###mwiewior:iterate first over normalization methods then over measures
    for(n in iNorm){
      
      if (n=='none') 
        nd.norm <- nd
      else if (n=='density') 
        nd.norm <- densityNormalize(nd)
      else if (n=='minMax') 
        nd.norm <- min_maxNormalize(nd)
      else
        nd.norm <- nd
      
      nd.norm <- RleList2matrix(getData(nd.norm) )
     # nd.norm <- RleList2matrix(nd.norm@data)
      combos<-.perm.samples(n.samples=n.samples, n.subsamples=length(n1), n.perm=iNPerm)
      dd_all <- .aveSimpleSamples(nd.norm, n1, n2)
      
      for(m in iMeasure)
      {
        if (m=="DA") 
          out <- c(diff_area(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="diff_area")))
        else if (m=="QQ") 
          out <- c(qq_plot(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="qq_plot")))
        else if (m=="PP") 
          out <- c(pp_plot(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="pp_plot")))
        else if (m=="HD1") 
          out <- c(hump_diff1(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="hump_diff1")))
        else if (m=="HD2") 
          out <- c(hump_diff2(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="hump_diff2")))
        
        M<-length(which(out>out[1]))
        pVal<-permp(M,nperm=iNPerm,length(n1),length(n2),twosided=T)
        
        
        testResults[i,paste(m,n,sep="_")]<- pVal
      
        
        
      }
      
    }   
    
  }
  rownames(testResults)<-as.character(annot[1:numgenes,6])
  return(testResults)
}  


.camelTestParReg <- function(iBedFile, iCovdesc='covdesc',iT1, iT2, iNorm=c('none','density','minMax'),  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iSizes = NULL, iNPerm=100)
{
  
  #mwiewior: reading BED
  annot<-read.table(iBedFile)
  numgenes <- dim(annot)[1]
  cvd <- read.table(iCovdesc)
  n.samples <- dim(cvd)[1]
  n1 <- which(cvd== iT1) 
  n2 <- which(cvd== iT2)
  
  #mwiewior:initiate output matrix
  #testResults <- matrix(nrow=numgenes, ncol= (length(iMeasure)* length(iNorm)) )
  testResults <- data.frame(t(rep(NA, (length(iMeasure)* length(iNorm)) ) ) )
  measureColNames <-c()
  for(n in iNorm){
    for(m in iMeasure){
      measureColNames<-cbind(measureColNames,paste(m,n,sep="_"))
    }
  }
  colnames(testResults)<-measureColNames
  testResults <- testResults[-1,]
  .loadLibs("parallel")
  .loadLibs("doParallel")
  .loadLibs("foreach")
  #automatic no of cores detection
  cl <- makeCluster(detectCores())
  #manual detection
  #cl <- makeCluster(6)
  registerDoParallel(cl)
  i=NULL
  testResults2<- foreach( i=1:nrow(annot), .export=c(".callCamelTestReg",".strand2number",".perm.samples",".shuffle",".aveSimpleSamples",
                                                    ".diff_stat",".ave_varSamples",".averageSamples",".rs2Counts",".localMax"), 
                         .packages=c("statmod","ampliQueso","rnaSeqMap"),.combine=rbind) %dopar% {
    .callCamelTestReg(annot, i,iCovdesc,n.samples, iSizes, iNorm, n1, iNPerm, n2, iMeasure,measureColNames,testResults)
  }
  
  stopCluster(cl)
  return(testResults2)
}  





.callCamelTestReg <- function (annot, i,iCovdesc, n.samples, iSizes, iNorm, n1, iNPerm, n2, iMeasure,measureColNames,testResults) {

    chr <- as.character(annot[i,1]);
    st <- as.numeric(annot[i,2]);
    end <- as.numeric(annot[i,3]);
    regName <-as.character(annot[i,6]);
    strand <- as.character(annot[i,4]);
    strand <- .strand2number(strand)
    #print("ch1")
    rs <- newSeqReads(chr,st,end,strand)
    rs <- getBamData(rs,1:n.samples, covdesc.file = iCovdesc)
    nd <- getCoverageFromRS(rs,1:n.samples)
    
    ###normalizacja
    if (!is.null(iSizes) )
      for (j in 1:n.samples) 
        setSAXPYData(nd,iSizes[j,2],j)
	#nd@data[[j]]<-nd@data[[j]]*iSizes[j,2] # tutaj jest normalizacja globalna wzgledem rozmiarow prob
    
    
    ###mwiewior:iterate first over normalization methods then over measures
    for(n in iNorm){
      
      if (n=='none') 
        nd.norm <- nd
      else if (n=='density') 
        nd.norm <- densityNormalize(nd)
      else if (n=='minMax') 
        nd.norm <- min_maxNormalize(nd)
      else
        nd.norm <- nd
      
      nd.norm <- RleList2matrix(getData(nd.norm) )
      #nd.norm <- RleList2matrix(nd.norm@data)
      combos<-.perm.samples(n.samples=n.samples, n.subsamples=length(n1), n.perm=iNPerm)
      dd_all <- .aveSimpleSamples(nd.norm, n1, n2)
      #print("ch3")
      #testResults=matrix(ncol=length(iMeasure)* length(iNorm),nrow=1)
      #testResults <- data.frame(t(rep(NA, (length(iMeasure)* length(iNorm)) ) ) )
      #colnames(testResults)<-measureColNames
      for(m in iMeasure)
      {
        if (m=="DA") 
          out <- c(diff_area(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="diff_area")))
        else if (m=="QQ") 
          out <- c(qq_plot(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="qq_plot")))
        else if (m=="PP") 
          out <- c(pp_plot(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="pp_plot")))
        else if (m=="HD1") 
          out <- c(hump_diff1(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="hump_diff1")))
        else if (m=="HD2") 
          out <- c(hump_diff2(dd_all), apply(combos,1,function(x) .shuffle(nd.norm,group1=x,measure="hump_diff2")))
        
        M<-length(which(out>out[1]))
        pVal<-permp(M,nperm=iNPerm,length(n1),length(n2),twosided=T)
        
        
        testResults[i,paste(m,n,sep="_")]<- pVal
        
        
      }
      
    }
    rownames(testResults)<-as.character(annot[i,6])
    return(testResults[i,])
}



.perm.samples<-function(n.samples=12,n.subsamples=6, n.perm=100)
{
  a <- 1:n.samples
  combos<-t(combn(n.samples, n.subsamples))
  B<-dim(combos)[1]
  
  if (n.perm>B) combos<-combos[sample(1:B,n.perm,replace=T),]
  
  if (n.perm<B) combos<-combos[sample(1:B,n.perm,replace=F),]
  
  combos
}

.shuffle<-function(nd,group1=1:6,measure="diff_area")
{
  n<-dim(nd)[2]
  group2<-c(1:n)[-group1]
  dd_all <- .aveSimpleSamples(nd,group1,group2)
  
  if (measure=="diff_area") out <- diff_area(dd_all)
  if (measure=="qq_plot") out <- qq_plot(dd_all)
  if (measure=="pp_plot") out <- pp_plot(dd_all)
  if (measure=="hump_diff1") out <- hump_diff1(dd_all)
  if (measure=="hump_diff2") out <- hump_diff2(dd_all)
  
  out
}


.aveSimpleSamples <- function(nd, idx1=1:4, idx2=5:8)
{
  avg1 <- rowMeans(nd[,idx1])
  avg2 <- rowMeans(nd[,idx2])
  dd1 <- cbind(avg1,avg2)
  dd1
}


.diff_stat<-function(nd,group1,group2)
{
  n1<-length(group1)
  n2<-length(group2)
  dd_all <- .ave_varSamples(nd,group1,group2)
  
  #  if (is.null(cconst))
  #  cconst <- max(abs(dd[, 1:2]))
  #  out <- out/(cconst * n)
  m <- dim(dd_all[[1]])[1]
  
  out <- sqrt((n1+n2-2)/(1/n1+1/n2))/m*sum((abs(dd_all$average[, 1] - dd_all$average[, 2]))/sqrt(dd_all$variance[, 1]+dd_all$variance[, 1]))
  
  out
}


.ave_varSamples <- function(nd, group1=group1, group2=group2)
{
  out<-list()
  dd_all <- RleList2matrix(getData(nd) )
  #dd_all <- RleList2matrix(nd@data)
  avg1 <- apply(dd_all[,group1],1,mean)
  avg2 <- apply(dd_all[,group2],1,mean)
  dd1 <- matrix(avg1,length(avg1),1)
  dd1 <- cbind(dd1,avg2)
  var1 <- apply(dd_all[,group1],1,function(x) sum((x-mean(x))^2))
  var2 <- apply(dd_all[,group2],1,function(x) sum((x-mean(x))^2))
  dd2 <- matrix(var1,length(var1),1)
  dd2 <- cbind(dd2,var2)
  
  out[["average"]]<-dd1
  out[["variance"]]<-dd2
  out
}


.averageSamples <- function(nd, idx1=1:4, idx2=5:8)
{
  dd_all <- RleList2matrix(getData(nd) )
  #dd_all <- RleList2matrix(nd@data)
  avg1 <- apply(dd_all[,idx1],1,mean)
  avg2 <- apply(dd_all[,idx2],1,mean)
  dd1 <- matrix(avg1,length(avg1),1)
  dd1 <- cbind(dd1,avg2)
  dd1
}

# x<-1:4
# varianceSamples <- function(nd, idx1=1:4, idx2=5:8)
# {
#   dd_all <- RleList2matrix(nd@data)
#   avg1 <- apply(dd_all[,idx1],1,mean)
#   avg2 <- apply(dd_all[,idx2],1,mean)
#   dd1 <- matrix(avg1,length(avg1),1)
#   dd1 <- cbind(dd1,avg2)
#   dd1
# }


.rs2Counts <- function(rs)
{
  out <- NULL
  rsData <- getData(rs)
  for (i in 1:length(rsData))
  {
    if (rsData[[i]][1,1]==0) out <- c(out, 0)
    else out <- c(out, dim(rsData[[i]])[1])
  }
  out
}

.localMax <- function(v){
  #finds local maxima and returns them as index
  vl <- length(v)
  idx <-NULL
  if (vl>2){
    for (i in 2:(vl-1)){
      if (v[i]>v[i-1] & v[i]>=v[i+1]) idx <- c(idx, i)
    }
  }
  idx
}
