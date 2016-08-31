Lin.reg.distance <- function(dist.m, emp.dist, pairs, tab.pop.pop, offhalf.list, offfull.list, relate.off.non.Mxy.mean, file.output, directory.name, out.name, inputdata, object, value, iteration)
{
  

  # Function caluclates CI intervals and exports plots and matrices

  dist.m <- as.dist(dist.m)
  emp.dist <- as.dist(1-emp.dist)
  
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
  
  man.out <- mantel(dist.m,emp.dist,permutations=pairs)
  
  offh <- 1-offhalf.list[!is.na(as.vector(offhalf.list))]
  offs <- 1-offfull.list[!is.na(as.vector(offfull.list))]
  non <- 1-relate.off.non.Mxy.mean[!is.na(as.vector(relate.off.non.Mxy.mean))]

if (file.output==TRUE)
{  
  pdf(paste(".","/",directory.name,"/","Total-Regression",tab.pop.pop[1,2],".pdf",sep=""))
  plot(dist.m, emp.dist, main="Regression of pairwise relatedness with geographic distance", xlab="Geographic distance", ylab=paste("1-Pairwise relatedness [",value,"]",sep=""), col=colors()[228])
  abline(lm(emp.dist~dist.m),lty="solid")
  
  lines(seq(range(dist.m)[1],range(dist.m)[2],0.1),rep(mean(offh),length(seq(range(dist.m)[1],range(dist.m)[2],0.1))),col="red",lty="dotdash")
  lines(seq(range(dist.m)[1],range(dist.m)[2],0.1),rep(mean(offs),length(seq(range(dist.m)[1],range(dist.m)[2],0.1))),col="blue",lty="dashed")
  lines(seq(range(dist.m)[1],range(dist.m)[2],0.1),rep(mean(non),length(seq(range(dist.m)[1],range(dist.m)[2],0.1))),lty="dotted")

  dev.off()

out.file <- file(paste(".","/",directory.name,"/","Total.Regression",tab.pop.pop[1,2],out.name,".txt",sep=""),"w")
  
writeLines(paste(
  "Demerelate - v.0.9-1", " ---","\n","Relatedness outputfile on file:", inputdata,"\n","Analysis had been made based on ",pairs," pairs using the ",value," estimator.","\n",
  if (value=="Bxy"){paste("Calculations are based on Li and Horvitz 1953. The values represent an indication on relatedness based on allele sharing.","\n", sep="")},
  if (value=="Mxy"){paste("Calculations are based on Bluoin et al. 1996. The values represent relatedness assessment based on genotype sharing.","\n", sep="")},
  if (value=="rxy"){paste("Calculations are based on Queller and Goodnight 1989. The values represent relatedness value corrected for total allele diversity.","\n", sep="")},
  if (value=="Sxy"){paste("Calculations are based on Lynch 1988. The values represent an indication on relatedness based on proportion of shared alleles.","\n", sep="")},
  if (value=="lxy"){paste("Calculations are based on Lynch and Ritland 1999. Final values represent relatedness estimates and are averaged via RE-RAT procedure.","\n", sep="")},
  if (value=="ritland"){paste("Calculations are based on Ritland 1996 with slight modifications based on Lynch and Ritland 1999.","\n", sep="")},
  if (value=="wang.fin"){paste("Calculations are based on Wang 2002. The estimates should only be used for finite populations.","\n", sep="")},
  if (value=="wang"){paste("Calculations are based on Wang 2002. The estimates is corrected for sample size and for expected frequency of unrelated individuals according to Li et al 1993.","\n", sep="")},
  if (value=="morans"){paste("Calculations are based on Hardy and Vekemans 1999 which describes morans I as estimator for genetic relatedness.","\n", sep="")},
  if (value=="morans.fin"){paste("Calculations are based on Hardy and Vekemans 1999 but ommiting sample size corrections. This value is only applicable for finite populations.","\n", sep="")},
  if (value=="loiselle"){paste("Calculations are based on Loiselle et al. 1995 ans are implemented as described by Hardy and
Vekemans 2015.","\n", sep="")},
  if (value=="Li"){paste("The estimator is calculated according to Li et al. 1993","\n", sep="")},
  "\n",
  "Calculations had been made for population:", as.character(tab.pop.pop[1,2]),"\n",
  "\n",
  "Summary on Mantel statistics","\n",
  "---","\n",
  man.out[[2]],"\n",
  "Statistic               :",round(man.out[[3]],3),"\n",
  "Number of Permutations  :",man.out[[6]],"\n",
  "Significance            :",round(man.out[[4]],3),sep=""),con=out.file)  

  writeLines(paste(
  "---","\n","\n",
  "\n","\n","\n","References","\n",
  "Blouin, M., Parsons, M., Lacaille, V. and Lotz, S. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.","\n",
  "Hardy, O.J. and Vekemans, X. (1999) Isolation by distance in a contiuous population: reconciliation between spatial autocorrelation analysis and population genetics models. Heredity, 83, 145-154.","\n",
  "Li, C.C., Weeks, D.E. and Chakravarti, A. (1993) Similarity of DNA fingerprints due to chance and relatedness. Human Heredity, 43, 45-52.","\n",
  "Li, C.C. and Horvitz, D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics, 5, 107-17.","\n",
  "Loiselle, B.A., Sork, V.L., Nason, J. and Graham, C. (1995) Spatial genetic structure of a tropical understory shrub, Psychotria officinalis (Rubiaceae). American Journal of Botany, 82, 1420-1425.","\n",
  "Lynch, M. (1988) Estimation of relatedness by DNA fingerprinting. Molecular Biology and Evolution, 5(5), 584-599.","\n",
  "Lynch, M. and Ritland, K. (1999) Estimation of pairwise relatedness with molecular markers. Genetics, 152, 1753-1766.","\n",
  "Mantel, N. (1967) The detection of disease clustering and a generalized regression approach. Cancer Research, 27, 209-220.","\n",
  "Oksanen, J. et al. (2013) vegan: Community Ecology Package. R package version 2.0-8.","\n",
  "Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.","\n",
  "Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.","\n",
  "Ritland, K. (1999) Estimators for pairwise relatedness and individual inbreeding coefficients. Genetics Research, 67, 175-185.","\n",
  "Wang, J. (2002) An estimator for pairwise relatedness using molecular markers. Genetics, 160, 1203-1215.","\n",
  sep=""),
  con=out.file
           )
  
close(out.file)

}
  return(man.out)
  
}
