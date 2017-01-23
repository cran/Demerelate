# Changed 23.5.2016
# Changed 16.8.2016 added pm.e as allelefreq for emp

stat.pops <- function(Thresholds, tab.pop.pop,  pairs, p.correct, directory.name, out.name, file.output, inputdata, object, value, iteration, pm, genotype.ref)

    {

    number.loci <- (ncol(tab.pop.pop)-2)/2
    empirical.share.ls <- vector("list",number.loci)
    empirical.share.ls.pm <- empirical.share.ls
    relate.non.X <- vector("list",number.loci)
    if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="wang" | value=="wang.fin" | value=="ritland")
    {
      empirical.share.ls.w <- vector("list",number.loci)
      relate.non.X.w <- vector("list",number.loci)
    }
    
    # Preparing vectors/dfs for applying
    data <- expand.grid(seq(1:length(tab.pop.pop[,1])),seq(1:length(tab.pop.pop[,1])))
    data <- data[as.numeric(data[,1]) <= as.numeric(data[,2]),]
    data <- data[data[,1]!=data[,2],]
  							   
    # Calculation of value for each locus in population tab.pop.pop
    for (i in 1:number.loci)
				{
				message(paste("---","Calculations for empirical values are performed for Locus",i,"----",Sys.time()),"\n")

				# Empirical share calculated for each locus
        
        empirical.share.ls[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,data, value, pm[i])
				if (value=="lxy" | value=="loiselle"){empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,data, value=paste(value,".w",sep=""), pm[[i]])}
        if (value=="ritland"){empirical.share.ls.w[[i]] <- length(pm[[i]])-1}
        if (value=="loiselle"){empirical.share.ls.w[[i]] <- sum(pm[[i]]*(1-pm[[i]]))}
        if (value=="morans" | value=="morans.fin"){empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,data, value="morans.w", pm[[i]])}
        if (value=="wang") {empirical.share.ls.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
        if (value=="wang.fin") {empirical.share.ls.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm[[i]])}
				names(empirical.share.ls)[i] <- paste(names(tab.pop.pop)[(i*2)+1],names(tab.pop.pop)[(i*2)+2],sep="-")
     							
				}
				
    			# Mean empirical
    if (value=="rxy" | value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="Li")
    {
			empirical.list <- rowMeans(do.call("cbind", empirical.share.ls), na.rm=TRUE)
			
    }
    
    if (value=="loiselle" | value=="morans" | value=="morans.fin"| value=="ritland")
    {
      empirical.list <- rowSums(do.call("cbind", empirical.share.ls), na.rm=TRUE)
      empirical.list.w <- rowSums(do.call("cbind", empirical.share.ls.w), na.rm=TRUE)
      empirical.list <- empirical.list/empirical.list.w
    }
    
    if (value=="lxy")
    {
      # RE
      empirical.list <- do.call("cbind",sapply(empirical.share.ls,function(x)(x[1])))
      empirical.list.w <- do.call("cbind",sapply(empirical.share.ls.w,function(x)(x[1])))
      empirical.list <- rowSums(empirical.list,na.rm=TRUE)
      empirical.list.w <- rowSums(empirical.list.w,na.rm=TRUE)
      empirical.list <- empirical.list/empirical.list.w
      
      # RAT
      empirical.list.rat <- do.call("cbind",sapply(empirical.share.ls,function(x)(x[2])))
      empirical.list.w.rat <- do.call("cbind",sapply(empirical.share.ls.w,function(x)(x[2]))) 
      empirical.list.rat <- rowSums(empirical.list.rat,na.rm=TRUE)
      empirical.list.w.rat <- rowSums(empirical.list.w.rat,na.rm=TRUE)
      empirical.list.rat <- empirical.list.rat/empirical.list.w.rat
      
      # AVERAGE
      empirical.list <- (empirical.list+empirical.list.rat)/2
      
    }
    
    if (value=="wang" | value=="wang.fin")
    {
    # weight for loci
    # According to frotran code of related, b-g and Pi are corrected for ul and average Pi and average b-g are corrected for 1/sum(1/ul)
    # Strangely average means here the sum f Pi and b-g ...?
    # Calculation is made for finite samples omitting equation 12-14 in wang2002 in wang.fin
    # Option wang takes the bias correction for sampling bias into account
    u <- unlist(lapply(seq(1:length(empirical.share.ls.w)),function(x){u<-empirical.share.ls.w[[x]][7]}))
    empirical.share.ls <- lapply(seq(1:length(empirical.share.ls)),function(x){empirical.share.ls[[x]] * 1/u[x]})
    empirical.share.ls <- Reduce("+",empirical.share.ls)
    empirical.share.ls <- empirical.share.ls*(1/(sum(1/u)))
    empirical.share.ls.w <- Reduce("+",empirical.share.ls.w)
    empirical.share.ls.w <- empirical.share.ls.w*(1/(sum(1/u)))
    
    empirical.list <- rowMeans(do.call("rbind",lapply(seq(1:length(empirical.share.ls[,1])), function(x){wang.compose(as=empirical.share.ls.w,Ps=empirical.share.ls[x,])})))
    
    
    }
    
    #Naming indivdual pairings
    
    names(empirical.list) <- apply(expand.grid(tab.pop.pop[,1],tab.pop.pop[,1]),1,paste,collapse="_")[as.numeric(row.names(data))]
    
# Testing on NAs
			if (length(table(is.na((empirical.list))))==2)
				{

warning(" ############## ERROR FULL STOP ################ ")
warning(" ############## ERROR FULL STOP ################ ")
warning(" ############## ERROR FULL STOP ################ ")
warning(" ## DON'T PANIK JUST READ FURTHER INSTRUCTIONS ## ")
warning("---","\n","\n")
warning("---","\n","\n")
warning(paste("An error occured in calculating allele sharings for population ",tab.pop.pop[1,2],"----",Sys.time(),sep=" "))
warning("---","\n","\n")
warning("Due to too much missing values in inputdata some individual sharings could not be calculated ----")
warning("Check on output for further information")
warning("---","\n","\n")
warning("NAs mark invalid combinations for individual allele sharings from individuals indicated in row and column names of 'error.individuals' ----")
warning("Please remove at least one individual indicated in 'error.individuals' or NAs from at least one individual from your input.txt and retry analysis")
warning("---","\n","\n")
warning("error.pairings (output=position in emp):")
warning("---","\n","\n")
error.individuals <- which(is.na(empirical.list),arr.ind=TRUE)
error.individuals
warning("---","\n","\n")
warning("If NAs not removable try starting analysis with mode cluster=FALSE; instead of cluster=TRUE (default)")
return(error.individuals)

if (file.output==TRUE)
{ write.table(file="error.individuals.txt",x=error.individuals)
}
stop()
}

    # Clusteranalysis
if (file.output==TRUE)
{
    # as.dist conversion
    empirical.dist <- data.frame(do.call("rbind",strsplit(names(empirical.list),"_")),empirical.list)
    empirical.dist <- with(empirical.dist,
                         structure(
                           empirical.dist[,3],
                           Size = length(unique(c(as.character(empirical.dist[,1]),as.character(empirical.dist[,2])))),
                           Labels = unique(c(as.character(empirical.dist[,1]),as.character(empirical.dist[,2]))),
                           Diag = FALSE,
                           Upper = FALSE,
                           method = "user",
                           class = "dist"))

		pdf(paste(".","/",directory.name,"/","Cluster",tab.pop.pop[1,2],out.name,".pdf",sep=""))


				par(cex=0.7,font=3)
				dis.plot <- plot(hclust(as.dist(1-empirical.dist),method="complete"),xlab="",ylab=paste("1-",value,sep=""),sub=paste("Inter-individual relatedness in population", tab.pop.pop[1,2],sep=" "), main=paste(value,"normalized dissimilarity"))
				abline(h=1-Thresholds[[1]], col="red", lty="dotdash")
				abline(h=1-Thresholds[[2]], col="blue", lty="dashed")

		hist(empirical.list, main = paste("Histogram of normalized mean allele sharing in", tab.pop.pop[1,2], sep=" "), xlab=paste("Empirical relatedness [",value,"]"))
		abline(v=Thresholds[[1]], col="red", lty="dotdash")
		abline(v=Thresholds[[2]], col="blue", lty="dashed")

		dev.off()
}

	    
    # length of reference population == empirical population
    data1 <- data.frame(seq(1:length(empirical.list[empirical.list!=NA])),seq(1:length(empirical.list[empirical.list!=NA])))
		# Random non-related for X-square
    for (i in 1:number.loci) 
    {
    
    if (genotype.ref==TRUE)
    {
    random.pairs.non.ls.X.1 <- data.frame(paste(rep("NON",length(empirical.list[empirical.list!=NA])),seq(1:length(empirical.list[empirical.list!=NA])),sep="-"),"NON",matrix(sapply(pm[i],function(x){sample(x=as.numeric(names(x)),prob=x,size=2*length(empirical.list[empirical.list!=NA]),replace=T)}),ncol=2))
    random.pairs.non.ls.X.2 <- data.frame(paste(rep("NON",length(empirical.list[empirical.list!=NA])),seq(1:length(empirical.list[empirical.list!=NA])),sep="-"),"NON",matrix(sapply(pm[i],function(x){sample(x=as.numeric(names(x)),prob=x,size=2*length(empirical.list[empirical.list!=NA]),replace=T)}),ncol=2))
    }
    
    if (genotype.ref==FALSE)
    {
    random.pairs.non.ls.X.1 <- data.frame(paste(rep("NON",length(empirical.list[empirical.list!=NA])),seq(1:length(empirical.list[empirical.list!=NA])),sep="-"),"NON",sapply(pm[[length(pm)]][,(i*2+1):(i*2+2)],function(x){sample(x,size=length(empirical.list[empirical.list!=NA]),replace=T)}))
    random.pairs.non.ls.X.2 <- data.frame(paste(rep("NON",length(empirical.list[empirical.list!=NA])),seq(1:length(empirical.list[empirical.list!=NA])),sep="-"),"NON",sapply(pm[[length(pm)]][,(i*2+1):(i*2+2)],function(x){sample(x,size=length(empirical.list[empirical.list!=NA]),replace=T)}))
    }
    
    relate.non.X[[i]] <- allele.sharing(random.pairs.non.ls.X.1, random.pairs.non.ls.X.2, 1, data1, value, pm[i])
    if (value=="lxy") {relate.non.X.w[[i]] <- allele.sharing(random.pairs.non.ls.X.1,random.pairs.non.ls.X.2,1,data1, value=paste(value,".w",sep=""), pm[[i]])}
    if (value=="loiselle") {relate.non.X.w[[i]] <- sum(pm[[i]]*(1-pm[[i]]))}
    if (value=="ritland") {relate.non.X.w[[i]] <- length(pm[[i]])-1}
    if (value=="morans" | value=="morans.fin") {relate.non.X.w[[i]] <- allele.sharing(random.pairs.non.ls.X.1,random.pairs.non.ls.X.2,1,data1, value="morans.w", pm[[i]])}
    if (value=="wang") {relate.non.X.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
    if (value=="wang.fin") {relate.non.X.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm[[i]])}
    }

		if (value=="Mxy" | value=="Sxy" | value=="Bxy" | value=="rxy" | value=="Li")
		{
				relate.non.X.mean <- rowMeans(do.call("cbind",relate.non.X),na.rm=TRUE)
		}
		
		if (value=="loiselle" | value=="morans.fin" | value=="morans" | value=="ritland")
		{
		  relate.non.X.mean <- rowSums(do.call("cbind",relate.non.X),na.rm=TRUE)
		  relate.non.X.mean.w <- rowSums(do.call("cbind",relate.non.X.w),na.rm=TRUE)
		  relate.non.X.mean <- relate.non.X.mean/relate.non.X.mean.w
		}

    if (value=="lxy")
    {
      ## NON
      relate.non.X.mean <- do.call("cbind",sapply(relate.non.X,function(x)(x[1])))
      relate.non.X.mean.w <- do.call("cbind",sapply(relate.non.X.w,function(x)(x[1])))
      relate.non.X.mean <- rowSums(relate.non.X.mean,na.rm=TRUE)
      relate.non.X.mean.w <- rowSums(relate.non.X.mean.w,na.rm=TRUE)
      relate.non.X.mean <- relate.non.X.mean/relate.non.X.mean.w
      
      # RAT
      relate.non.X.mean.rat <- do.call("cbind",sapply(relate.non.X,function(x)(x[2])))
      relate.non.X.mean.w.rat <- do.call("cbind",sapply(relate.non.X.w,function(x)(x[2]))) 
      relate.non.X.mean.rat <- rowSums(relate.non.X.mean.rat,na.rm=TRUE)
      relate.non.X.mean.w.rat <- rowSums(relate.non.X.mean.w.rat,na.rm=TRUE)
      relate.non.X.mean.rat <- relate.non.X.mean.rat/relate.non.X.mean.w.rat
      
      # AVERAGE
      relate.non.X.mean <- (relate.non.X.mean+relate.non.X.mean.rat)/2
 
    }

    if (value=="wang" | value=="wang.fin")
    {
    # weight for loci
    # According to frotran code of related, b-g and Pi are corrected for ul and average Pi and average b-g are corrected for 1/sum(1/ul)
    # Strangely average means here the sum f Pi and b-g ...?
    # Calculation is made for finite samples omitting equation 12-14 in wang2002 in wang.fin
    # Option wang takes the bias correction for sampling bias into account
    u <- unlist(lapply(seq(1:length(relate.non.X.w)),function(x){u<-relate.non.X.w[[x]][7]}))
    relate.non.X <- lapply(seq(1:length(relate.non.X)),function(x){relate.non.X[[x]] * 1/u[x]})
    relate.non.X <- Reduce("+",relate.non.X)#/length(empirical.share.ls)
    relate.non.X <- relate.non.X*(1/(sum(1/u)))
    relate.non.X.w <- Reduce("+",relate.non.X.w)#/length(empirical.share.ls.w)
    relate.non.X.w <- relate.non.X.w*(1/(sum(1/u)))
  
    relate.non.X.mean <- rowMeans(do.call("rbind",lapply(seq(1:length(relate.non.X[,1])), function(x){wang.compose(as=relate.non.X.w,Ps=relate.non.X[x,])})))
    
    }

    relate.non.X.mean.values <- relate.non.X.mean

		relate.non.X.mean[relate.non.X.mean.values>=Thresholds[1,2]] <- "FS"
		relate.non.X.mean[relate.non.X.mean.values>=Thresholds[1,1] & relate.non.X.mean.values<Thresholds[1,2]] <- "HS"
		relate.non.X.mean[relate.non.X.mean.values<Thresholds[1,1]] <- "NON"
    
    empirical.list.values <- empirical.list

    # Matrix conversion
    empirical.list[empirical.list.values>=Thresholds[1,2]] <- "FS"
    empirical.list[empirical.list.values>=Thresholds[1,1] & empirical.list.values<Thresholds[1,2]] <- "HS"
    empirical.list[empirical.list.values<Thresholds[1,1]] <- "NON"
		
    
		# P Statistics
		# only FS
		emp <- 0
		non <- 0
		if (!is.na(as.numeric(table(empirical.list)["FS"]))==TRUE){emp <- as.numeric(table(empirical.list)["FS"])}
		if (!is.na(as.numeric(table(relate.non.X.mean)["FS"]))==TRUE){non <- as.numeric(table(relate.non.X.mean)["FS"])}
		stat.p <- prop.test(c(emp,non), c(sum(table(empirical.list)), sum(table(relate.non.X.mean))), conf.level=0.95,correct=p.correct)

		# FS + HS
		emp <- 0
		non <- 0
		if (!is.na(sum(as.numeric(table(empirical.list)["HS"])==TRUE, as.numeric(table(empirical.list)["FS"])==TRUE, na.rm=TRUE))){emp <- sum(as.numeric(table (empirical.list)["FS"]), as.numeric(table(empirical.list)["HS"]),na.rm=TRUE)}
		if (!is.na(sum(as.numeric(table(relate.non.X.mean)["HS"])==TRUE, as.numeric(table(relate.non.X.mean)["FS"])==TRUE, na.rm=TRUE))){non <- sum(as.numeric(table(relate.non.X.mean)["FS"]), as.numeric(table(relate.non.X.mean)["HS"]),na.rm=TRUE)}
		stat.p.HS <- prop.test(c(emp,non),c(sum(table(empirical.list)),sum(table(relate.non.X.mean))),conf.level=0.95,correct=p.correct)
        
		f.p <- as.data.frame(rbind(c(stat.p[[4]][1], stat.p[[4]][2],stat.p[[1]],stat.p[[2]],stat.p[[3]],stat.p[[6]][1],stat.p[[6]][2]),c(stat.p.HS[[4]][1],stat.p.HS[[4]][2],stat.p.HS[[1]],stat.p.HS[[2]],stat.p.HS[[3]],stat.p.HS[[6]][1],stat.p.HS[[6]][2])))
		f.p <- round(f.p,3)		
		row.names(f.p) <- c("Full Siblings","Full+Half Siblings")
		colnames(f.p) <- c("Observed","Expected","Chi^2","d.f.","p-value","0.95 Lower CI","0.95 Upper CI")

if (file.output==TRUE)
{
  out.file <- file(paste(".","/",directory.name,"/","Relate.mean",tab.pop.pop[1,2],out.name,".txt",sep=""),"w")
  writeLines(
    paste(
      "Demerelate - v.0.9-2", " ---","\n","Relatedness outputfile on file: ", inputdata,"\n","Analysis had been made based on ",pairs," pairs using the ",value," estimator.","\n",
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
      if (value=="loiselle"){paste("Calculations are based on Loiselle et al. 1995 and are implemented as described by Hardy and
Vekemans 2015.","\n", sep="")},
      if (value=="Li"){paste("The estimator is calculated according to Li et al. 1993","\n", sep="")},
  "For further information mind References at the end of this file.","\n","\n",
  "Calculations had been made for population:", as.character(tab.pop.pop[1,2]),"\n",
  "\n",
  "Relatedness Thresholds","\n",
  "---", sep=""),con=out.file)
  write.table(round(Thresholds,3), file=out.file, append=T, sep="\t", quote=F, row.names=F) 
              writeLines(paste("---","\n","\n",
  "Relatedness calculations","\n",
  "---","\n",
	"Observed frequencies of full siblings (FS), half siblings (HS) and non related pairs (NON)","\n",sep=""),con=out.file)
              write.table(table(empirical.list),file=out.file,append=T,sep="\t", quote=F, row.names=F, col.names=F)
              
  writeLines(paste(
  "---","\n","\n",
  "---","\n",
	"Expected frequencies of full siblings (FS), half siblings (HS) and non related pairs (NON)","\n",sep=""),con=out.file)
              
  write.table(table(relate.non.X.mean),file=out.file,append=T,sep="\t", quote=F, row.names=F, col.names=F)
  
  writeLines(paste("---","\n","\n","Chisquare Statistics","\n","---","\n","\n",sep=""),con=out.file)
    
  write.table(f.p,file=out.file,append=T,sep="\t", quote=F, col.names = NA)
              
  writeLines(paste("\n",

  "----------------------------------------------------------------------------------------------------","\n",

      
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
  con=out.file)
    
    close(out.file)

}

		
    
    
      out.stat <- list(empirical.list.values, f.p)
      names(out.stat) <- c("Empirical_List", "Chi-square statistics")
      return(out.stat)
    
  
    
  }
