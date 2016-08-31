# Changed 23.5.2016

relate.calc <- function(tab.pop, pairs, file.output, value, directory.name, ref.pop)
	{	

		number.loci <- (ncol(tab.pop)-2)/2
   	empirical.share.ls <- vector("list",number.loci)
		relate.full.mean <- vector("list",number.loci)
		relate.half.mean <- vector("list",number.loci)
		relate.off.non <- vector("list",number.loci)
		random.pairs.fsib.ls <- vector("list",2)
		random.pairs.hsib1.ls	<- vector("list",2)
		random.pairs.hsib2.ls <- vector("list",2)
    
      x <- seq(3,ncol(ref.pop),2)
      pm <- lapply(x,function(x){table(c(ref.pop[,x],ref.pop[,x+1]))/length(c(ref.pop[,x],ref.pop[,x+1]))})
      names(pm) <- lapply(x,function(x){sum(complete.cases(cbind(ref.pop[,x],ref.pop[,x+1])))})
        
    
    if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="wang" | value=="wang.fin" | value=="ritland")
    {
      empirical.share.ls.w <- vector("list",number.loci)
      relate.full.mean.w <- vector("list",number.loci)
      relate.half.mean.w <- vector("list",number.loci)
      relate.off.non.w <- vector("list",number.loci)
    }
    
		# Random pairs for fullsibs overall
		random.pairs.fsib.ls[[1]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.fsib.ls[[1]][,1]<-paste("FS-",random.pairs.fsib.ls[[1]][,1],"-",seq(1:length(random.pairs.fsib.ls[[1]][,1])),sep="")
		random.pairs.fsib.ls[[2]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
		random.pairs.fsib.ls[[2]][,1]<-paste("FS-",random.pairs.fsib.ls[[2]][,1],"-",seq((length(random.pairs.fsib.ls[[1]][,1])+1):(length(random.pairs.fsib.ls[[2]][,1])+length(random.pairs.fsib.ls[[1]][,1]))),sep="")
		# Random pairs for halfsibs overall 
		random.pairs.hsib1.ls[[1]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.hsib1.ls[[1]][,1]<-paste("HS-",random.pairs.hsib1.ls[[1]][,1],"-",seq(1:length(random.pairs.hsib1.ls[[1]][,1])),sep="")
		random.pairs.hsib2.ls[[1]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.hsib2.ls[[1]][,1]<-paste("HS-",random.pairs.hsib2.ls[[1]][,1],"-",seq((length(random.pairs.hsib1.ls[[1]][,1])+1):(length(random.pairs.hsib1.ls[[1]][,1])+length(random.pairs.hsib2.ls[[1]][,1]))),sep="")
		random.pairs.hsib2.ls[[2]] <- ref.pop[sample(1:nrow(ref.pop),pairs,replace=TRUE),]
    random.pairs.hsib2.ls[[2]][,1]<-paste("HS-",random.pairs.fsib.ls[[1]][,1],"-",seq((length(random.pairs.hsib1.ls[[1]][,1])+length(random.pairs.hsib2.ls[[1]][,1])+1):(length(random.pairs.hsib1.ls[[1]][,1])+length(random.pairs.hsib2.ls[[1]][,1])+length(random.pairs.hsib2.ls[[2]][,1]))),sep="")
		
       for (i in 1:number.loci)	
				{
		
		
		message(paste("---","Calculations are performed for Locus",i,"----",Sys.time(),"\n"))
		
			# 1. Empirisches sharing for each locus calculated
				empirical.share.ls[[i]] <- allele.sharing(tab.pop,tab.pop, i, FALSE, value, pm)
				names(empirical.share.ls)[i] <- names(tab.pop)[(i*2)+1]
		    if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="ritland") {empirical.share.ls.w[[i]] <- allele.sharing(tab.pop,tab.pop, i, FALSE, value=paste(value,".w",sep=""), pm)}
				if (value=="wang") {empirical.share.ls.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
				if (value=="wang.fin") {empirical.share.ls.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm)}
				
		 	# 2. Random  pairs calculated
				  random.pairs.non.ls <- random.pairs(ref.pop,(i*2)+1,pairs)
    
		      off.full.ls.1 <- offspring(random.pairs.fsib.ls[[1]],random.pairs.fsib.ls[[2]],(i*2)+1, pairs)
		      off.full.ls.2 <- offspring(random.pairs.fsib.ls[[1]],random.pairs.fsib.ls[[2]],(i*2)+1, pairs)
		
		      off.half.ls.1 <- offspring(random.pairs.hsib1.ls[[1]],random.pairs.hsib2.ls[[1]],(i*2)+1, pairs)
		      off.half.ls.2 <- offspring(random.pairs.hsib1.ls[[1]],random.pairs.hsib2.ls[[2]],(i*2)+1, pairs)
    		  		
			# 3. Offsprings for reference are calculated
    
		      relate.off.non[[i]] <- allele.sharing(random.pairs.non.ls[[1]],random.pairs.non.ls[[2]],1,TRUE, value, pm[i])
		      if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="ritland") {relate.off.non.w[[i]] <- allele.sharing(random.pairs.non.ls[[1]],random.pairs.non.ls[[2]],1,TRUE, value=paste(value,".w",sep=""), pm[i])}
		      if (value=="wang") {relate.off.non.w[[i]] <- wang.w(allele.column=1, ref.pop=pm[i])}
		      if (value=="wang.fin") {relate.off.non.w[[i]] <- wang.fin.w(allele.column=1, ref.pop=pm[i])}
    
		      relate.full.mean[[i]] <- allele.sharing(off.full.ls.1,off.full.ls.2,1,TRUE, value, pm[i])
		      if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="ritland") {relate.full.mean.w[[i]] <- allele.sharing(off.full.ls.1,off.full.ls.2,1,TRUE, value=paste(value,".w",sep=""), pm[i])}
		      if (value=="wang") {relate.full.mean.w[[i]] <- wang.w(allele.column=1, ref.pop=pm[i])}
		      if (value=="wang.fin") {relate.full.mean.w[[i]] <- wang.fin.w(allele.column=1, ref.pop=pm[i])}
    
					relate.half.mean[[i]] <- allele.sharing(off.half.ls.1,off.half.ls.2,1,TRUE, value, pm[i])
					if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="ritland") {relate.half.mean.w[[i]] <- allele.sharing(off.half.ls.1,off.half.ls.2,1,TRUE, value=paste(value,".w",sep=""), pm[i])}
					if (value=="wang") {relate.half.mean.w[[i]] <- wang.w(allele.column=1, ref.pop=pm[i])}
					if (value=="wang.fin") {relate.half.mean.w[[i]] <- wang.fin.w(allele.column=1, ref.pop=pm[i])}
    
		   		}
    
 
			# Empirical
      if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li")
      {
        
      empirical.Mxy.mean <- rowMeans(do.call("cbind",empirical.share.ls),na.rm=TRUE)  
			relate.off.full.Mxy.mean <- rowMeans(do.call("cbind",relate.full.mean),na.rm=TRUE)
			relate.off.non.Mxy.mean <- rowMeans(do.call("cbind",relate.off.non),na.rm=TRUE)	
			relate.off.half.Mxy.mean <- rowMeans(do.call("cbind",relate.half.mean),na.rm=TRUE)	

      }
 
		  if (value=="loiselle"  | value=="morans.fin" | value=="morans" | value=="ritland")
		  {
		  
		  empirical.Mxy.mean <- rowSums(do.call("cbind",empirical.share.ls),na.rm=TRUE)
		  empirical.Mxy.mean.w <- rowSums(do.call("cbind",empirical.share.ls.w),na.rm=TRUE)
		  empirical.Mxy.mean <- empirical.Mxy.mean/empirical.Mxy.mean.w

		  relate.off.full.Mxy.mean <- rowSums(do.call("cbind",relate.full.mean),na.rm=TRUE)
		  relate.off.full.Mxy.mean.w <- rowSums(do.call("cbind",relate.full.mean.w),na.rm=TRUE)
		  relate.off.full.Mxy.mean <- relate.off.full.Mxy.mean/relate.off.full.Mxy.mean.w
      
		  relate.off.non.Mxy.mean <- rowSums(do.call("cbind",relate.off.non),na.rm=TRUE)	
		  relate.off.non.Mxy.mean.w <- rowSums(do.call("cbind",relate.off.non.w),na.rm=TRUE)  
		  relate.off.non.Mxy.mean <- relate.off.non.Mxy.mean/relate.off.non.Mxy.mean.w
      
      relate.off.half.Mxy.mean <- rowSums(do.call("cbind",relate.half.mean),na.rm=TRUE)
		  relate.off.half.Mxy.mean.w <- rowSums(do.call("cbind",relate.half.mean.w),na.rm=TRUE)  
		  relate.off.half.Mxy.mean <- relate.off.half.Mxy.mean/relate.off.half.Mxy.mean.w
      
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
		  empirical.Mxy.mean <- (empirical.list+empirical.list.rat)/2
      
      ## FS
		  relate.off.full.Mxy.mean <- do.call("cbind",sapply(relate.full.mean,function(x)(x[1])))
		  relate.off.full.Mxy.mean.w <- do.call("cbind",sapply(relate.full.mean.w,function(x)(x[1])))
		  relate.off.full.Mxy.mean <- rowSums(relate.off.full.Mxy.mean,na.rm=TRUE)
		  relate.off.full.Mxy.mean.w <- rowSums(relate.off.full.Mxy.mean.w,na.rm=TRUE)
		  relate.off.full.Mxy.mean <- relate.off.full.Mxy.mean/relate.off.full.Mxy.mean.w
		  
		  # RAT
		  relate.off.full.Mxy.mean.rat <- do.call("cbind",sapply(relate.full.mean,function(x)(x[2])))
		  relate.off.full.Mxy.mean.w.rat <- do.call("cbind",sapply(relate.full.mean.w,function(x)(x[2]))) 
		  relate.off.full.Mxy.mean.rat <- rowSums(relate.off.full.Mxy.mean.rat,na.rm=TRUE)
		  relate.off.full.Mxy.mean.w.rat <- rowSums(relate.off.full.Mxy.mean.w.rat,na.rm=TRUE)
		  relate.off.full.Mxy.mean.rat <- relate.off.full.Mxy.mean.rat/relate.off.full.Mxy.mean.w.rat
		  
		  # AVERAGE
		  relate.off.full.Mxy.mean <- (relate.off.full.Mxy.mean+relate.off.full.Mxy.mean.rat)/2
      
      ## NON
		  relate.off.non.Mxy.mean <- do.call("cbind",sapply(relate.off.non,function(x)(x[1])))
		  relate.off.non.Mxy.mean.w <- do.call("cbind",sapply(relate.off.non.w,function(x)(x[1])))
		  relate.off.non.Mxy.mean <- rowSums(relate.off.non.Mxy.mean,na.rm=TRUE)
		  relate.off.non.Mxy.mean.w <- rowSums(relate.off.non.Mxy.mean.w,na.rm=TRUE)
		  relate.off.non.Mxy.mean <- relate.off.non.Mxy.mean/relate.off.non.Mxy.mean.w
		  
		  # RAT
		  relate.off.non.Mxy.mean.rat <- do.call("cbind",sapply(relate.off.non,function(x)(x[2])))
		  relate.off.non.Mxy.mean.w.rat <- do.call("cbind",sapply(relate.off.non.w,function(x)(x[2]))) 
		  relate.off.non.Mxy.mean.rat <- rowSums(relate.off.non.Mxy.mean.rat,na.rm=TRUE)
		  relate.off.non.Mxy.mean.w.rat <- rowSums(relate.off.non.Mxy.mean.w.rat,na.rm=TRUE)
		  relate.off.non.Mxy.mean.rat <- relate.off.non.Mxy.mean.rat/relate.off.non.Mxy.mean.w.rat
		  
		  # AVERAGE
		  relate.off.non.Mxy.mean <- (relate.off.non.Mxy.mean+relate.off.non.Mxy.mean.rat)/2

      ## NON
		  relate.off.half.Mxy.mean <- do.call("cbind",sapply(relate.half.mean,function(x)(x[1])))
		  relate.off.half.Mxy.mean.w <- do.call("cbind",sapply(relate.half.mean.w,function(x)(x[1])))
		  relate.off.half.Mxy.mean <- rowSums(relate.off.half.Mxy.mean,na.rm=TRUE)
		  relate.off.half.Mxy.mean.w <- rowSums(relate.off.half.Mxy.mean.w,na.rm=TRUE)
		  relate.off.half.Mxy.mean <- relate.off.half.Mxy.mean/relate.off.half.Mxy.mean.w
		  
		  # RAT
		  relate.off.half.Mxy.mean.rat <- do.call("cbind",sapply(relate.half.mean,function(x)(x[2])))
		  relate.off.half.Mxy.mean.w.rat <- do.call("cbind",sapply(relate.half.mean.w,function(x)(x[2]))) 
		  relate.off.half.Mxy.mean.rat <- rowSums(relate.off.half.Mxy.mean.rat,na.rm=TRUE)
		  relate.off.half.Mxy.mean.w.rat <- rowSums(relate.off.half.Mxy.mean.w.rat,na.rm=TRUE)
		  relate.off.half.Mxy.mean.rat <- relate.off.half.Mxy.mean.rat/relate.off.half.Mxy.mean.w.rat
		  
		  # AVERAGE
		  relate.off.half.Mxy.mean <- (relate.off.half.Mxy.mean+relate.off.half.Mxy.mean.rat)/2
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
	    
	    empirical.Mxy.mean <- rowMeans(do.call("rbind",lapply(seq(1:length(empirical.share.ls[,1])), function(x){wang.compose(as=empirical.share.ls.w,Ps=empirical.share.ls[x,])})))
	    
      
      # Full Sib randomized
	    u <- unlist(lapply(seq(1:length(relate.full.mean.w)),function(x){u<-relate.full.mean.w[[x]][7]}))
	    relate.full.mean <- lapply(seq(1:length(relate.full.mean)),function(x){relate.full.mean[[x]] * 1/u[x]})
	    relate.full.mean <- Reduce("+",relate.full.mean)
	    relate.full.mean <- relate.full.mean*(1/(sum(1/u)))
	    relate.full.mean.w <- Reduce("+",relate.full.mean.w)
	    relate.full.mean.w <- relate.full.mean.w*(1/(sum(1/u)))
	    
	    relate.off.full.Mxy.mean <- rowMeans(do.call("rbind",lapply(seq(1:length(relate.full.mean[,1])), function(x){wang.compose(as=relate.full.mean.w,Ps=relate.full.mean[x,])})))
	  
      # Non Sib randomized
	    u <- unlist(lapply(seq(1:length(relate.off.non.w)),function(x){u<-relate.off.non.w[[x]][7]}))
	    relate.off.non <- lapply(seq(1:length(relate.off.non)),function(x){relate.off.non[[x]] * 1/u[x]})
	    relate.off.non <- Reduce("+",relate.off.non)
	    relate.off.non <- relate.off.non*(1/(sum(1/u)))
	    relate.off.non.w <- Reduce("+",relate.off.non.w)
	    relate.off.non.w <- relate.off.non.w*(1/(sum(1/u)))
	    
	    relate.off.non.Mxy.mean <- rowMeans(do.call("rbind",lapply(seq(1:length(relate.off.non[,1])), function(x){wang.compose(as=relate.off.non.w,Ps=relate.off.non[x,])})))
	   
      # Half sib randomized
	    u <- unlist(lapply(seq(1:length(relate.half.mean.w)),function(x){u<-relate.half.mean.w[[x]][7]}))
	    relate.half.mean <- lapply(seq(1:length(relate.half.mean)),function(x){relate.half.mean[[x]] * 1/u[x]})
	    relate.half.mean <- Reduce("+",relate.half.mean)
	    relate.half.mean <- relate.half.mean*(1/(sum(1/u)))
	    relate.half.mean.w <- Reduce("+",relate.half.mean.w)
	    relate.half.mean.w <- relate.half.mean.w*(1/(sum(1/u)))
	    
	    relate.off.half.Mxy.mean <- rowMeans(do.call("rbind",lapply(seq(1:length(relate.half.mean[,1])), function(x){wang.compose(as=relate.half.mean.w,Ps=relate.half.mean[x,])})))
	  
	    }

		
		# Calculating multiple logistic regression
		return.glm <- glm.prep(empirical.Mxy.mean, relate.off.full.Mxy.mean, relate.off.half.Mxy.mean, relate.off.non.Mxy.mean)
    half <- return.glm[[1]]
    sumlrm <- return.glm[[2]]
		  
	  # Threshold calculation
    full <- (sumlrm[[1]][1]-sumlrm[[1]][2])/(sumlrm[[1]][4]-sumlrm[[1]][3])		
		Thres <- data.frame(half,full)
    		    
if (file.output==TRUE)
{
    write.table(file=paste(".","/",directory.name,"/","Random.Fullsib.relatedness.overall.txt",sep=""),x=relate.off.full.Mxy.mean, quote=FALSE, sep=" ", col.names=value)
    write.table(file=paste(".","/",directory.name,"/","Random.Halfsib.relatedness.overall.txt",sep=""),x=relate.off.half.Mxy.mean, quote=FALSE, sep=" ", col.names=value)
    write.table(file=paste(".","/",directory.name,"/","Random.NonRelated.relatedness.overall.txt",sep=""),x=relate.off.non.Mxy.mean, quote=FALSE, sep=" ", col.names=value)
}
		# Assigning Output
		relate.return <- list(empirical.Mxy.mean, relate.off.full.Mxy.mean, relate.off.half.Mxy.mean, relate.off.non.Mxy.mean, Thres)
    names(relate.return) <- c("Relatedness_Empirical"," Randomized_Fullssibs","Randomized_Halfsibs","Ranodmized_Non","Thresholds")
		return(relate.return)
			
	}	

