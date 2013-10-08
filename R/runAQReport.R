runAQReport <- function(iCovdesc, iBedFile, iRefSeqFile=NULL,iGroup, iT1, iT2,iTopN=5,iMinQual, iReportFormat="pdf", iReportType="article",iReportPath=NULL, iVerbose=TRUE,  iVCF=NULL, iParallel=TRUE) {
    
    #constant package name
    PACKAGE_NAME="ampliQueso"
    
    #load libs needed for the report prepration
    .loadLibs(PACKAGE_NAME)
    .loadLibs("knitr")
    .loadLibs("rgl")
    .loadLibs("genefilter")
    .loadLibs("rnaSeqMap")
    .loadLibs("statmod")
    #package home directory
    gvWorkdir = getwd()
    gvPackagePath<-path.package(PACKAGE_NAME)
    #setwd(gvPackagePath)
    
    
  
    #set the global variables to be sure that results will be consistent
    gvCovdesc=iCovdesc
    gvCovdescT<-read.table(gvCovdesc,header=TRUE)
    gvBedFile=iBedFile
    gvRefSeqFile=iRefSeqFile
    gvGroup=iGroup
    gvT1=iT1
    gvT2=iT2
    gvReportFormat=iReportFormat
    gvReportType=iReportType
    #show commands used for the report generation
    gvVerbose=iVerbose
    
    
    gvMinQual=iMinQual
    #Camel Test measure used for computing p-value and sorting
    gvCamelMeasurePVal="DA_density"
    
    reportParams<-t(data.frame(.escSpecialChars(iCovdesc),.escSpecialChars(iBedFile),.escSpecialChars(iRefSeqFile),.escSpecialChars(iGroup),.escSpecialChars(iT1),.escSpecialChars(iT2),iReportFormat,iReportType))
    rownames(reportParams)<-c("Covdesc file","BED File","Ref. seq. file","Group var.", "Group 1","Group 2", "Report format", "Report type")
    
    #fc thresholds vector
    gvAbsFoldChangeThr<-c(0.8,1,1.5,2,5)
    #p-value thresholds vector
    gvPValueThr <-c(1/10^(3:0))
    #report columns names
    #####old T-test
    #repColNames <-c("p.value.thr","p.value.adj.pct\\\\%","fc0.8","fc1","fc1.5","fc2","fc3")
    repColNames <-c("p.value.thr","p.value.adj.pct\\%","fc0.8","fc1","fc1.5","fc2","fc5")
    
    #topN genes cut-off
    gvTopN=iTopN
    
    #######RNA AmpliSeq analysis
    ###Counts
    print("Calculating counts...")
    countTable<-getCountTable(covdesc=gvCovdesc, bedFile=gvBedFile)
    if( ncol(countTable) > 2 ){
        logMeanCountT1 <- rowMeans(log2(countTable[,1:ncol(countTable)/2]+0.0001))
        logMeanCountT2 <- rowMeans(log2(countTable[,(ncol(countTable)/2 + 1):ncol(countTable)]+0.0001))
    }else if( ncol(countTable) == 2){
        logMeanCountT1 <- log2(countTable[,1:ncol(countTable)/2]+0.0001)
        logMeanCountT2 <- log2(countTable[,(ncol(countTable)/2 + 1):ncol(countTable)]+0.0001)  
    }else 
      stop("Minimun sample number should >=2")
    cc <- cbind(logMeanCountT1,logMeanCountT2)
    foldChange <- cc[,1] - cc[,2]
    fact<-factor(read.table(gvCovdesc)[,1])
    cleanCountTable<-apply(countTable,c(1,2), as.numeric)
    ###Cameltests
    if(ncol(countTable) > 2){
       print("Calculating camel tests...")
       if(iParallel==TRUE)
          camelTestTable <- camelTest(iBedFile=gvBedFile,iCovdesc=gvCovdesc, iT1=gvT1, iT2=gvT2,iParallel=TRUE)
       else
          camelTestTable <- camelTest(iBedFile=gvBedFile,iCovdesc=gvCovdesc, iT1=gvT1, iT2=gvT2,iParallel=FALSE)
    }
    else if( ncol(countTable) == 2)
    {
      #quick an dirty fix for 2 sample
      tTestTable<-rowttests(cleanCountTable,fact)
      pAdjVector <- p.adjust(tTestTable[,3])
      tTestTable <-cbind(tTestTable,pAdjVector)
      colnames(tTestTable)=c(colnames(tTestTable[,1:3]),"p.value.adj")
      camelTestTable<-tTestTable
      colnames(camelTestTable)[colnames(camelTestTable) =="p.value"]=gvCamelMeasurePVal
    }
    else
      stop("Minimun sample number should >=2")
      
    #####old T-test
    #tTestTable<-rowttests(cleanCountTable,fact)
    #pAdjVector <- p.adjust(tTestTable[,3])
    #tTestTable <-cbind(tTestTable,pAdjVector)
    #colnames(tTestTable)=c(colnames(tTestTable[1:3]),"p.value.adj")
    
    
    absFoldChange<-abs(foldChange)
    #####old T-test
    #countTableTran<-cbind(tTestTable,absFoldChange)
    countTableTran<-cbind(absFoldChange,camelTestTable)
    #sort by p.value.adj
    #####old T-test
    #countTableTranSort<-countTableTran[order(countTableTran$p.value),]
    
    ####choose only one Camel test measure for p-value sorting and computation
    colnames(countTableTran)[colnames(countTableTran) == gvCamelMeasurePVal] = "p.value"
    #countTableTranSort<-countTableTran[order(countTableTran[["p.value"]]),]
    countTableTranSort<-countTableTran[order(countTableTran[,"p.value"]),]
    
    #define output counts table
    countRepTable <-data.frame(matrix(NA, nrow = length(gvPValueThr), ncol = length(gvAbsFoldChangeThr)+2 ) )
    colnames(countRepTable)<-repColNames
    countRepTable[,1]<-gvPValueThr
    pLo<-0
    fcLo<-0
    for(p in gvPValueThr) {
      for (fc in gvAbsFoldChangeThr){
        val<-nrow(countTableTranSort[(countTableTranSort[,"p.value"]> pLo & countTableTranSort[,"p.value"] <= p 
                                     & countTableTranSort[,"absFoldChange"] >fcLo & countTableTranSort[,"absFoldChange"] <= fc), ,drop = FALSE ])
        #print(val)
        if(!is.null(val))
          countRepTable[countRepTable$p.value.thr==p , colnames(countRepTable)==paste("fc",as.character(fc),sep="" ) ]<-val
        else
          countRepTable[countRepTable$p.value.thr==p , colnames(countRepTable)==paste("fc",as.character(fc),sep="" ) ]<-0
        
        #fcLo<-fc
      }
      fcLo<-0
      #pLo<-p
    }
    pPcts<-rowSums(countRepTable[,3:ncol(countRepTable)] )/sum(countRepTable[3:ncol(countRepTable)]) *100
    countRepTable[,2]<-pPcts
    colnames(countRepTable)<-repColNames
    
    
    #
    ###Volcano
    #######################################3moved to the report template
    #plot(foldChange, -log(tTestTable$p.value))
    ###Heatmap
    #top N-genes with the lowest p-value 
    #moved to the report template
    #heatmap.2(countTable[rownames(countTableTranSort[1:gvTopN,]),],trace="none")
    ####Counts distributions
    binsNum=10
    #get fullpaths
    #moved to the report template
#     legendPaths<-strsplit( colnames(countTable),"/")
#     legendSampl<-rep("",ncol(countTable))
#     y<-rep(0,binsNum)
#     for(s in 1:ncol(countTable)){   
#         #remove outliers by taking 1st and 99th percentile of counts
#         countPercentiles<-quantile(countTable[,s],probs=c(0.01,0.99))
#         #50+1 for seq function due to its formula: by = ((to - from)/(length.out - 1)
#         bins <-seq(countPercentiles[1],countPercentiles[2],length.out=binsNum+1)
#         for(b in 2:length(bins)){         
#             y[b]<-length(countTable[,s][countTable[,s]>=bins[b-1] & countTable[,s]<bins[b]])   
#         }
#         if(s == 1){
#           plot(bins,y,type="l", col=(rainbow(s)[s]),ylim<-c(0,countPercentiles[2]),ylab="counts")
#         }
#         else{
#           lines(bins,y, col=(rainbow(s)[s])   )     
#         }
#         points(bins,y,col=(rainbow(s)[s]), cex=0.5, pch=1)
#         #get filenames from fullpaths
#         legendSampl[s]<-unlist(legendPaths[s])[length(unlist(legendPaths[s]))]
#     }
#       legend(legend=legendSampl,x=2000,y=150,lty=c(1,1),col=rainbow(ncol(countTable)) )
    ###check if colours are generated correctly
    
    
    #############Camels
    #call the parallel version
    print("Calculating camel measures...")
    if(iParallel==TRUE)
        camelsTable <-compareCoveragesReg(iBedFile=gvBedFile, iGroup=gvGroup, iT1=gvT1, iT2=gvT2, iCovdesc=gvCovdesc,iParallel=TRUE)
    else
        camelsTable <-compareCoveragesReg(iBedFile=gvBedFile, iGroup=gvGroup, iT1=gvT1, iT2=gvT2, iCovdesc=gvCovdesc,iParallel=FALSE)
    rownames(camelsTable)<-.escSpecialChars(rownames(countTable))
    camelsRankTable<-cbind(rank(camelsTable[,1]),rank(camelsTable[,2]), rank(camelsTable[,3]),rank(camelsTable[,4]),rank(camelsTable[,5]))
    camelsRankTable<-cbind(camelsRankTable, rowSums(camelsRankTable))
    camelsRankTableSort<-camelsRankTable[order(camelsRankTable[,6],decreasing=FALSE),]
    
    colnames(camelsRankTableSort)<-c("DA","QQ", "PP", "HD1", "HD2","Rank sum")
    camelsRankTableSort<-round(camelsRankTableSort)
    
    
    ################SNPs
    #call the parallel version
    #with BED file if refSeqFile not NULL
    if(!is.null(gvRefSeqFile))
    {
      print("Calling SNPs...")
        legendPaths<-strsplit( colnames(countTable),"/")
        if(iParallel==TRUE)
            aqSNPL<-getSNP(covdesc=gvCovdesc,minQual=gvMinQual,refSeqFile=iRefSeqFile,bedFile=iBedFile,iParallel=TRUE)
        else
          aqSNPL<-getSNP(covdesc=gvCovdesc,minQual=gvMinQual,refSeqFile=iRefSeqFile,bedFile=iBedFile,iParallel=FALSE)
        if(!is.null(aqSNPL) && length(aqSNPL) > 0) {
            #prepare SNP summary
            #add BAM filename as SNP lists names
            names(aqSNPL)<-foreach( i=1:nrow(gvCovdescT), .export=c("gvCovdescT","legendPaths")) %do%(unlist(legendPaths[i])[length(unlist(legendPaths[i]))])
            #merge SNP lists from the same group
            #get the id of the last element from the group column
            groupLengths<-rle(as.character(gvCovdescT[,1]))
            #merge SNP into groups
            snpGroups <- list()
            snpGroupsAggr <-list()
            snpGroupsAggrSort <- list()
            shift=0
            for(i in 1:length(groupLengths$lengths) ){
              #define an empty data frame with proper column names
              snpGroups[i][[1]]<- data.frame(t(rep(NA,length(colnames(aqSNPL[1][[1]])) )))
              colnames( snpGroups[i][[1]]) <- colnames(aqSNPL[1][[1]])
              snpGroups[i][[1]] <- snpGroups[i][[1]][-1,]
              
              for(j in (1+shift):(groupLengths$lengths[i]+shift) ){  
                snpGroups[i][[1]] <- merge( snpGroups[i][[1]],aqSNPL[j][[1]], all=TRUE)
                if(j == shift+groupLengths$lengths[i])
                  shift=j
              }
              if( nrow(snpGroups[i][[1]])> 0 ) 
                {
                  snpGroupsAggr[i][[1]] <-aggregate(snpGroups[i][[1]]$sequence, list(snpGroups[i][[1]]$sequence,snpGroups[i][[1]]$start,snpGroups[i][[1]]$end), FUN="length")
                  snpGroupsAggrSort[i][[1]] <-snpGroupsAggr[i][[1]][order(snpGroupsAggr[i][[1]][,4],decreasing=TRUE),]
                  colnames(snpGroupsAggrSort[i][[1]]) <- c("chr","start","end","count")
              }
            }
            if ( length(snpGroupsAggrSort)>1 && !is.null(snpGroupsAggrSort[[2]]) && !is.null(snpGroupsAggrSort[[1 ]]) )
            {
                snpDiffSick<-.setdiff.data.frame(snpGroupsAggrSort[[2]][,1:3],snpGroupsAggrSort[[1]][,1:3])
                snpDiffSickCount<-merge(snpDiffSick,snpGroupsAggrSort[[2]])
                snpDiffSickCountSort<-snpDiffSickCount[order(snpDiffSickCount[,4],decreasing=TRUE),]
                snpDiffSickCountFilter <- snpDiffSickCountSort[snpDiffSickCountSort[,4] >1, ]
                #find BAM files containing the SNPs to be reported
                snpDiffSickCountFilter<- .getSNPSamples(aqSNPL,snpDiffSickCountFilter)
                
                snpDiffHealthy<-.setdiff.data.frame(snpGroupsAggrSort[[1]][,1:3],snpGroupsAggrSort[[2]][,1:3])
                snpDiffHealthyCount<-merge(snpDiffHealthy,snpGroupsAggrSort[[1]])
                snpDiffHealthyCountSort<-snpDiffHealthyCount[order(snpDiffHealthyCount[,4],decreasing=TRUE),]
                snpDiffHealthyCountFilter <- snpDiffHealthyCountSort[snpDiffHealthyCountSort[,4] >1, ]
                snpDiffHealthyCountFilter<- .getSNPSamples(aqSNPL,snpDiffHealthyCountFilter)
            }
            else
            {
              snpDiffSickCountFilter <-data.frame()
              snpDiffHealthyCountFilter <- data.frame()
              
            }
        
        }
    
    
    }
    
    if(is.null(iReportPath))
      iReportPath<-getwd()
    #run the report with the proper template
    if (tolower(gvReportFormat) == "pdf"  && tolower(gvReportType) == "article"){
        #use aq_article.Rnw template
        knit(paste(gvPackagePath,"/reports/aq_article.Rnw",sep="")) 
        tools::texi2dvi("aq_article.tex", pdf=TRUE)
        file.rename("aq_article.pdf",paste(iReportPath,"/","aq_article.pdf",sep=""))
        #setwd(gvWorkdir)
    }
    else if(tolower(gvReportFormat) == "pdf"  && tolower(gvReportType) == "beamer"){
        #use aq_beamer.Rnw template
        knit(paste(gvPackagePath,"/reports/aq_beamer.Rnw",sep="")) 
        tools::texi2dvi("aq_beamer.tex", pdf=TRUE)
        file.rename("aq_beamer.pdf",paste(iReportPath,"/","aq_beamer.pdf",sep=""))
        #setwd(gvWorkdir)
    }
    else if (tolower(gvReportFormat) == "html"){
        #use aq_html.Rhtml template
        knit(paste(gvPackagePath,"/reports/aq_html.Rhtml",sep=""),output=paste(iReportPath,"/","aq_html.Rhtml",sep=""))
    }
    else 
      stop("Report format or type not supported.")
  
}
.getSNPSamples <- function(snpDF,sampleDF){
  if(nrow(sampleDF) > 0){
    
    sampleDF["files"]<-""
    for(i in 1:nrow(sampleDF)){
      chr <- as.character(sampleDF[i,][[1]])
      st <- sampleDF[i,][[2]]
      en <- sampleDF[i,][[3]]
      for(j in 1:length(snpDF)){
        sampleSNP <-snpDF[j][[1]]  
        if( nrow(sampleSNP[ as.character(sampleSNP$sequence) ==chr & sampleSNP$start==st & sampleSNP$end == en,,drop = FALSE]) == 1 )
          sampleDF[i,5]=paste(sampleDF[i,5],names(snpDF[j]),sep="," )
      }
    }
    sampleDF[,5]<-substring(.escSpecialChars(sampleDF[,5]),2)
    return(sampleDF)
  }
}