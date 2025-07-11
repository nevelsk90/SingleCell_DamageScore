tempo.out.list <- readRDS( file="tempo.out_listSCSN.rda")

### plot PDS vs maxP
tempo.out <- tempo.out.list[[1]]
tempo.out.maxP <- Reduce( rbind , lapply( seq(tempo.out), function(ii){
  
  sampleName <- sub("__.*" , "",  names(tempo.out)[ii])
  dataSet <-  sub(".*__" , "",  names(tempo.out)[ii])
  XX <- data.frame( tempo.maxP= apply( tempo.out[[ii]], 1, max),
                    tempo.maxPtime= apply( tempo.out[[ii]], 1, which.max ),
                    sample= sampleName ,
                    datSet=  dataSet )
}) )
tempo.out.maxP$PDS <- listSCSN_merged$PDS[ rownames(tempo.out.maxP) ]
tempo.out.maxP$gtype <- listSCSN_merged$gtypeDE[ rownames(tempo.out.maxP) ]
tempo.out.maxP$group <- listSCSN_merged$group[ rownames(tempo.out.maxP) ]
tempo.out.maxP$sampleLbl <- paste0( tempo.out.maxP$sample , "..", tempo.out.maxP$group)

## deviation from the likely time
toPlot <- Reduce( rbind, lapply( unique(tempo.out.maxP$sample), function(ss){
  print(ss)
  datt <- tempo.out.maxP[tempo.out.maxP$sample==ss,]
  dd <- density( datt$tempo.maxPtime , window = "optcosine",na.rm = T)
  # plot( dd)
  X<-which.max( dd$y)
  timeMost <- dd$x[X]
  
  datt$timeMost.diff <-  abs( unlist(
    lapply( datt$tempo.maxPtime, circadian_time_difference , ct1=timeMost)) )
  
  # df  <- data.frame( timeMost.diff=datt$timeMost.diff , PDS=datt$PDS)
  # pca <- prcomp(~PDS + timeMost.diff, data = df)
  # slp <- with(pca, rotation[2,1] / rotation[1,1])
  # int <- with(pca, center[2] - slp*center[1])
  # datt$slp<- slp
  # datt$int<- int
  return(datt)
}) )

### load  processed tempo results
llgg <- readRDS( file="tempo.out.list_maxP.Ptime.rda")

# separate data and plots
ppll <- lapply( llgg ,"[[",1)
ddll <-  lapply( llgg ,"[[",2)
names(ppll) <- names(ddll) <- names(llgg)

ddll  <- Reduce(rbind, ddll)
ddll$maxPvsPDS.rho <- as.numeric(ddll$maxPvsPDS.rho)
ddll$maxPvsPDS <- as.numeric( ddll$maxPvsPDS)
ddll$maxPtimevsPDS.rho <- as.numeric(ddll$maxPtimevsPDS.rho)
ddll$maxPtimevsPDS <- as.numeric( ddll$maxPtimevsPDS)



### plot summary stat for a selected methods
{
  # toPlot <- ddll[ !is.na(ddll$maxPvsPDS.rho.sig) | !is.na(ddll$maxPtimevsPDS.rho.sig), ]
  toPlot <- ddll
  
  toPlot$sigStat1 <- ifelse( toPlot$maxPvsPDS < 0.05, "sig", "nonsig")
  toPlot$sigStat2 <- ifelse( toPlot$maxPtimevsPDS < 0.05, "sig", "nonsig")
  toPlot.indv <- toPlot[ toPlot$testRun=="counts_core.clock.Clock2",] 
  # toPlot.indv <- toPlot[ toPlot$testRun=="counts_core.clock.Arntl",] 
  
  
  toPlot.indv$sampleVis1 <- ifelse( toPlot.indv$maxPvsPDS<0.05 , 
                                    paste0( toPlot.indv$datSet , "_._", toPlot.indv$group , toPlot.indv$sample) , NA)
  toPlot.indv$sampleVis2 <- ifelse( toPlot.indv$maxPtimevsPDS<0.05 , 
                                    paste0( toPlot.indv$datSet , "_._", toPlot.indv$group , toPlot.indv$sample) , NA)
  
  gg1<- ggplot( data = toPlot.indv, aes( x= gtype, y= maxPvsPDS.rho ,  shape=type, color=sigStat1))  + 
    geom_label_repel(aes(label = sampleVis1), color="grey30", size=3, label.size = NA)+
    geom_jitter(width = 0.1, alpha=0.5, size=5)+ theme_bw()+ 
    scale_color_tableau()+ theme( text = element_text( size = 20))
  gg2<- ggplot( data = toPlot.indv, aes( x= gtype, y= maxPtimevsPDS.rho ,  shape=type, color=sigStat2))  + 
    geom_label_repel(aes(label = sampleVis2), color="grey30", size=3,label.size = NA )+
    geom_jitter(width = 0.1, alpha=0.5, size=5)+ theme_bw( )+ 
    scale_color_tableau() + theme( text = element_text( size = 20))
  
  ccppl <- cowplot::plot_grid( plotlist = list(gg1,gg2))  
  ccppl
  
  pdf(width = 12 , height = 8, file = "/data/user/tpadvits/PROJECTS/PDS/GRN/tempo/tempo.out_listSCSN.counts_core.Clock.pdf")
  ccppl
  dev.off()
}