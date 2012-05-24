Lin.reg.distance <- function(dist.m, emp.dist, pairs, tab.pop.pop, offhalf.list, offfull.list, relate.off.non.Mxy.mean, file.output, directory.name, out.name, inputdata, object, value, iteration)
{
  

  # Function caluclates CI intervals and exports plots and matrices

  dist.m <- as.dist(dist.m)
  emp.dist <- as.dist(emp.dist)
  
  # Bootstrap for dist.m and emp number determined by pairs
  bt.sample <- sample(1:length(as.dist(dist.m)), pairs, replace=TRUE)
  
  reg <- lm(emp.dist[bt.sample]~dist.m[bt.sample])
  
  sum.reg <- summary(reg)
  r.sum <- data.frame(rbind(sum.reg[[8]], sum.reg[[9]]),row.names=c("R-square","R-square-adjusted"))
  names(r.sum) <- " "
  reg.out <- list(as.data.frame(t(as.matrix(summary(sum.reg[[3]])))), as.data.frame(sum.reg[[4]]), r.sum)
  names(reg.out) <- c("Residuals of regression", "Coefficients of estimate", "Coefficients of Correlation")
  
  reg.pred <- predict(reg, interval="confidence")
  new.dist <- data.frame(dist.m[bt.sample],reg.pred[,2],reg.pred[,3])
  
  offh <- offhalf.list[!is.na(as.vector(offhalf.list))]
  offs <- offfull.list[!is.na(as.vector(offfull.list))]
  non <- relate.off.non.Mxy.mean[!is.na(as.vector(relate.off.non.Mxy.mean))]

if (file.output==TRUE)
{  
  pdf(paste(".","/",directory.name,"/","Total-Regression",tab.pop.pop[1,2],".pdf",sep=""))
  plot(dist.m, emp.dist, main="Regression of pairwise relatedness with geographic distance", xlab="Geographic distance", ylab="Pairwise relatedness", col=colors()[228])
  abline(lm(emp.dist~dist.m),lty=50)
  
  lines(new.dist[order(new.dist[,1]),][,1], new.dist[order(new.dist[,1]),][,2], lty=40)
  lines(new.dist[order(new.dist[,1]),][,1], new.dist[order(new.dist[,1]),][,3], lty=40)
  
  polygon(c(range(dist.m)[1]-10,range(dist.m)[2]+10,range(dist.m)[2]+10,range(dist.m)[1]-10),c(mean(non)-var(non),mean(non)-var(non),mean(non)+var(non),mean(non)+var(non)),border=NA,col=colors()[228],density=100)
  polygon(c(range(dist.m)[1]-10,range(dist.m)[2]+10,range(dist.m)[2]+10,range(dist.m)[1]-10),c(mean(offs)-var(offs),mean(offs)-var(offs),mean(offs)+var(offs),mean(offs)+var(offs)),border=NA,col=colors()[307],density=100)
  polygon(c(range(dist.m)[1]-10,range(dist.m)[2]+10,range(dist.m)[2]+10,range(dist.m)[1]-10),c(mean(offh)-var(offh),mean(offh)-var(offh),mean(offh)+var(offh),mean(offh)+var(offh)),border=NA,col=colors()[274],density=100)
  
  lines(seq(range(dist.m)[1],range(dist.m)[2],0.1),rep(mean(offh),length(seq(range(dist.m)[1],range(dist.m)[2],0.1))),lty=40)
  lines(seq(range(dist.m)[1],range(dist.m)[2],0.1),rep(mean(offs),length(seq(range(dist.m)[1],range(dist.m)[2],0.1))),lty=40)
  lines(seq(range(dist.m)[1],range(dist.m)[2],0.1),rep(mean(non),length(seq(range(dist.m)[1],range(dist.m)[2],0.1))),lty=40)                                                       
  
  dev.off()

out.file <- file(paste(".","/",directory.name,"/","Total.Regression",tab.pop.pop[1,2],out.name,".txt",sep=""),"w")
  
writeLines(paste(
  "Demerelate - v.0.8", "---","\n",
  "Relatedness outputfile on file:", inputdata,"\n",
  "Analysis had been made using", iteration,"iterations","and",pairs,"pairs","using the",value,"estimator.","\n",
  if (value=="Bxy"){paste("Calculations are based on Li and Horvitz 1953. The values represent an indication on relatedness based on allele sharing.","\n", sep=" ")},
  if (value=="Mxy"){paste("Calculations are based on Bluoin et al. 1996. The values represent relatedness assessment based on genotype sharing.","\n", sep=" ")},
  if (value=="rxy"){paste("Calculations are based on Queller and Goodnight 1989. The values represent relatedness value corrected for total allele diversity.","\n", sep=" ")},
  "For further information mind References at the end of this file.","\n","\n",
  "Calculations had been made for population:", as.character(tab.pop.pop[1,2]),"\n",
  "\n",
  "Summary on relatedness regression","\n",
  "---","\n",sep=" "),con=out.file)
  
  write.table(reg.out[[1]],file=out.file,append=T,sep="\t", quote=F)
  write.table(reg.out[[2]],file=out.file,append=T,sep="\t", quote=F)
  write.table(reg.out[[3]],file=out.file,append=T,sep="\t", quote=F, col.names=F)
  

  writeLines(paste(
  "---","\n","\n",
  "\n","\n","\n","References","\n",
  "Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.","\n",
  "Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.","\n",
  "Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.","\n",sep=" "),
  con=out.file
           )
  
close(out.file)

}
  return(reg.out)
  
}
