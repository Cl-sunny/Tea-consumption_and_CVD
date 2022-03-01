rm(list=ls())
library(forestplot)
###Reading in data for Figure 3###
data <- read.csv ('data for Figure 3.csv') 

tabletext<-cbind( 
  c("Outcomes",as.character(data[,1])),
  c("Data source",as.character(data[,2])),
  c("SNPs",as.character(data[,3])),
  c("OR (95% CI)",as.character(data[,4])),
  c("P-value",as.character(data[,5])))
tabletext

cochrane_from_rmeta <- structure(list(
  mean =c(NA, data[,6]),
  lower =c(NA, data[,7]),
  upper =c(NA, data[,8])),
  .Names =c("mean", "lower", "upper"),
  row.names =c(NA, -13L),
  class ="data.frame")

forestplot(tabletext,cochrane_from_rmeta,new_page =T,
           graph.pos = 4, zero=1, 
           hrzl_lines=list("2" = gpar(lwd=2.5, col="#272727"),
                           "4" =gpar(lwd=2.5,lty=2, col ="#5B5B5B"),
                           "5" = gpar(lwd=2.5, col="#272727"),
                           "7" =gpar(lwd=2.5,lty=2, col ="#5B5B5B"),
                           "8" = gpar(lwd=2.5, col="#272727"),      
                           "10" =gpar(lwd=2.5,lty=2, col ="#5B5B5B"),
                           "11" = gpar(lwd=2.5, col="#272727"),
                           "13" =gpar(lwd=2.5,lty=2, col ="#5B5B5B"),
                           "14" = gpar(lwd=2.5, col="black")
                           ),
           txt_gp=fpTxtGp(label=gpar(cex=1.8),
                          ticks=gpar(cex=1.8),
                          xlab=gpar(cex=1.8),
                          title=gpar(cex=1.8)),
           is.summary=c(T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,F),
           xlog=T,
           xticks =c(0.85,0.9,0.95,1.0,1.05,1.1),
           clip=c(0.8,1.4),
           graphwidth = unit(0.2,"npc"), lineheight = unit(.17,"npc"), 
           boxsize = 0.2, 
           fn.ci_norm="fpDrawNormalCI",
           lwd.ci=3, ci.vertices=T, ci.vertices.height = 0.1,
           lwd.xaxis=2.5,
           lwd.zero=2.5,
           title="",
           col=fpColors(box="#FDC453",line="#FDC453",zero="#ADADAD",
                       summary="#9ADBC5",hrz_lines ="#808080"),
           align=c("l","l","c","c","c","c"))

