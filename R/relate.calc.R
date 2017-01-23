# Changed 23.5.2016

relate.calc <- function(tab.pop, pairs, file.output, value, directory.name, pm)
	{	

		number.loci <- (ncol(tab.pop)-2)/2
		relate.full.mean <- vector("list",number.loci)
		relate.half.mean <- vector("list",number.loci)
		relate.off.non <- vector("list",number.loci)
		random.pairs.fsib.ls <- vector("list",2)
		random.pairs.hsib1.ls	<- vector("list",2)
		random.pairs.hsib2.ls <- vector("list",2)
        
    
    if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="wang" | value=="wang.fin" | value=="ritland")
    {
      relate.full.mean.w <- vector("list",number.loci)
      relate.half.mean.w <- vector("list",number.loci)
      relate.off.non.w <- vector("list",number.loci)
    }

		
		### sample from pm...
		# Random pairs for fullsibs overall
		random.pairs.fsib.ls[[1]] <- data.frame(paste(rep("FSIB",pairs),seq(1:pairs),sep="-"),"FSIB",matrix(sapply(pm[-((length(pm)-1):length(pm))],function(x){sample(as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2*number.loci))
		random.pairs.fsib.ls[[2]] <- data.frame(paste(rep("FSIB",pairs),seq(1:pairs),sep="-"),"FSIB",matrix(sapply(pm[-((length(pm)-1):length(pm))],function(x){sample(as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2*number.loci))
		# Random pairs for halfsibs overall 
		random.pairs.hsib1.ls[[1]] <- data.frame(paste(rep("HSIB",pairs),seq(1:pairs),sep="-"),"HSIB",matrix(sapply(pm[-((length(pm)-1):length(pm))],function(x){sample(as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2*number.loci))
		random.pairs.hsib2.ls[[1]] <- data.frame(paste(rep("HSIB",pairs),seq(1:pairs),sep="-"),"HSIB",matrix(sapply(pm[-((length(pm)-1):length(pm))],function(x){sample(as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2*number.loci))
		random.pairs.hsib2.ls[[2]] <- data.frame(paste(rep("HSIB",pairs),seq(1:pairs),sep="-"),"HSIB",matrix(sapply(pm[-((length(pm)-1):length(pm))],function(x){sample(as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2*number.loci))
		
		
		
       for (i in 1:number.loci)	
				{
		
		
		message(paste("---","Calculations are performed for Locus",i,"----",Sys.time(),"\n"))
         
          random.pairs.non.ls.1 <- data.frame(paste(rep("NON",pairs),seq(1:pairs),sep="-"),"NON",matrix(sapply(pm[i],function(x){sample(x=as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2))
          random.pairs.non.ls.2 <- data.frame(paste(rep("NON",pairs),seq(1:pairs),sep="-"),"NON",matrix(sapply(pm[i],function(x){sample(x=as.numeric(names(x)),prob=x,size=pairs,replace=T)}),ncol=2))
         
		      off.full.ls.1 <- offspring(random.pairs.fsib.ls[[1]],random.pairs.fsib.ls[[2]],(i*2)+1, pairs)
		      off.full.ls.2 <- offspring(random.pairs.fsib.ls[[1]],random.pairs.fsib.ls[[2]],(i*2)+1, pairs)
		
		      off.half.ls.1 <- offspring(random.pairs.hsib1.ls[[1]],random.pairs.hsib2.ls[[1]],(i*2)+1, pairs)
		      off.half.ls.2 <- offspring(random.pairs.hsib1.ls[[1]],random.pairs.hsib2.ls[[2]],(i*2)+1, pairs)
    		  		
			# 3. Offsprings for reference are calculated
 ### for all allele column 1 instead of 1 and pm[i] !!Problem da pm aus overall refpop        
		      data1 <- data.frame(seq(1:pairs), seq(1:pairs))
		      relate.off.non[[i]] <- allele.sharing(random.pairs.non.ls.1,random.pairs.non.ls.2,1,data1, value, pm[i])
		      if (value=="lxy") {relate.off.non.w[[i]] <- allele.sharing(random.pairs.non.ls.1,random.pairs.non.ls.2,1,data1, value=paste(value,".w",sep=""), pm[[i]])}
		      if (value=="loiselle") {relate.off.non.w[[i]] <- sum(pm[[i]]*(1-pm[[i]]))}
		      if (value=="ritland") {relate.off.non.w[[i]] <- length(pm[[i]])-1}
		      if (value=="morans" | value=="morans.fin") {relate.off.non.w[[i]] <- allele.sharing(random.pairs.non.ls.1,random.pairs.non.ls.2,1,data1, value="morans.w", pm[[i]])}
		      if (value=="wang") {relate.off.non.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
		      if (value=="wang.fin") {relate.off.non.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm[[i]])}
    
		      relate.full.mean[[i]] <- allele.sharing(off.full.ls.1, off.full.ls.2, 1, data1, value, pm[i])
		      if (value=="lxy") {relate.full.mean.w[[i]] <- allele.sharing(off.full.ls.1,off.full.ls.2,1,data1, value=paste(value,".w",sep=""), pm[[i]])}
		      if (value=="loiselle") {relate.full.mean.w[[i]] <- sum(pm[[i]]*(1-pm[[i]]))}
		      if (value=="ritland") {relate.full.mean.w[[i]] <- length(pm[[i]])-1}
		      if (value=="morans" | value=="morans.fin") {relate.full.mean.w[[i]] <- allele.sharing(off.full.ls.1,off.full.ls.2,1,data1, value="morans.w", pm[[i]])}
		      if (value=="wang") {relate.full.mean.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
		      if (value=="wang.fin") {relate.full.mean.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm[[i]])}
    
					relate.half.mean[[i]] <- allele.sharing(off.half.ls.1, off.half.ls.2, 1, data1, value, pm[i])
					if (value=="lxy") {relate.half.mean.w[[i]] <- allele.sharing(off.half.ls.1,off.half.ls.2,1, data1, value=paste(value,".w",sep=""), pm[[i]])}
					if (value=="loiselle") {relate.half.mean.w[[i]] <- sum(pm[[i]]*(1-pm[[i]]))}
					if (value=="ritland") {relate.half.mean.w[[i]] <- length(pm[[i]])-1}
					if (value=="morans" | value=="morans.fin") {relate.half.mean.w[[i]] <- allele.sharing(off.half.ls.1,off.half.ls.2,1, data1, value="morans.w", pm[[i]])}
					if (value=="wang") {relate.half.mean.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
					if (value=="wang.fin") {relate.half.mean.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm[[i]])}
    
		   		}
    
 
			# Empirical
      if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li")
      {
      
			relate.off.full.Mxy.mean <- rowMeans(do.call("cbind",relate.full.mean),na.rm=TRUE)
			relate.off.non.Mxy.mean <- rowMeans(do.call("cbind",relate.off.non),na.rm=TRUE)	
			relate.off.half.Mxy.mean <- rowMeans(do.call("cbind",relate.half.mean),na.rm=TRUE)	

      }
 
		  if (value=="loiselle"  | value=="morans.fin" | value=="morans" | value=="ritland")
		  {

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
      ## FS
		  relate.off.full.Mxy.mean <- do.call("cbind",lapply(relate.full.mean,function(x)(.subset2(x,1))))
		  relate.off.full.Mxy.mean.w <- do.call("cbind",lapply(relate.full.mean.w,function(x)(.subset2(x,1))))
		  relate.off.full.Mxy.mean <- rowSums(relate.off.full.Mxy.mean,na.rm=TRUE)
		  relate.off.full.Mxy.mean.w <- rowSums(relate.off.full.Mxy.mean.w,na.rm=TRUE)
		  relate.off.full.Mxy.mean <- relate.off.full.Mxy.mean/relate.off.full.Mxy.mean.w
		  
		  # RAT
		  relate.off.full.Mxy.mean.rat <- do.call("cbind",lapply(relate.full.mean,function(x)(.subset2(x,2))))
		  relate.off.full.Mxy.mean.w.rat <- do.call("cbind",lapply(relate.full.mean.w,function(x)(.subset2(x,2))))
		  relate.off.full.Mxy.mean.rat <- rowSums(relate.off.full.Mxy.mean.rat,na.rm=TRUE)
		  relate.off.full.Mxy.mean.w.rat <- rowSums(relate.off.full.Mxy.mean.w.rat,na.rm=TRUE)
		  relate.off.full.Mxy.mean.rat <- relate.off.full.Mxy.mean.rat/relate.off.full.Mxy.mean.w.rat
		  
		  # AVERAGE
		  relate.off.full.Mxy.mean <- (relate.off.full.Mxy.mean+relate.off.full.Mxy.mean.rat)/2
      
      ## NON
		  relate.off.non.Mxy.mean <- do.call("cbind",lapply(relate.off.non,function(x)(.subset2(x,1))))
		  relate.off.non.Mxy.mean.w <- do.call("cbind",lapply(relate.off.non.w,function(x)(.subset2(x,1))))
		  relate.off.non.Mxy.mean <- rowSums(relate.off.non.Mxy.mean,na.rm=TRUE)
		  relate.off.non.Mxy.mean.w <- rowSums(relate.off.non.Mxy.mean.w,na.rm=TRUE)
		  relate.off.non.Mxy.mean <- relate.off.non.Mxy.mean/relate.off.non.Mxy.mean.w
		  
		  # RAT
		  relate.off.non.Mxy.mean.rat <- do.call("cbind",lapply(relate.off.non,function(x)(.subset2(x,2))))
		  relate.off.non.Mxy.mean.w.rat <- do.call("cbind",lapply(relate.off.non.w,function(x)(.subset2(x,2)))) 
		  relate.off.non.Mxy.mean.rat <- rowSums(relate.off.non.Mxy.mean.rat,na.rm=TRUE)
		  relate.off.non.Mxy.mean.w.rat <- rowSums(relate.off.non.Mxy.mean.w.rat,na.rm=TRUE)
		  relate.off.non.Mxy.mean.rat <- relate.off.non.Mxy.mean.rat/relate.off.non.Mxy.mean.w.rat
		  
		  # AVERAGE
		  relate.off.non.Mxy.mean <- (relate.off.non.Mxy.mean+relate.off.non.Mxy.mean.rat)/2

      ## NON
		  relate.off.half.Mxy.mean <- do.call("cbind",lapply(relate.half.mean,function(x)(.subset2(x,1))))
		  relate.off.half.Mxy.mean.w <- do.call("cbind",lapply(relate.half.mean.w,function(x)(.subset2(x,1))))
		  relate.off.half.Mxy.mean <- rowSums(relate.off.half.Mxy.mean,na.rm=TRUE)
		  relate.off.half.Mxy.mean.w <- rowSums(relate.off.half.Mxy.mean.w,na.rm=TRUE)
		  relate.off.half.Mxy.mean <- relate.off.half.Mxy.mean/relate.off.half.Mxy.mean.w
		  
		  # RAT
		  relate.off.half.Mxy.mean.rat <- do.call("cbind",lapply(relate.half.mean,function(x)(.subset2(x,2))))
		  relate.off.half.Mxy.mean.w.rat <- do.call("cbind",lapply(relate.half.mean.w,function(x)(.subset2(x,2))))
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
		return.glm <- glm.prep(relate.off.full.Mxy.mean, relate.off.half.Mxy.mean, relate.off.non.Mxy.mean)
    half <- .subset2(return.glm,1)
    sumlrm <- .subset2(return.glm,2)
		  
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
		relate.return <- list(relate.off.full.Mxy.mean, relate.off.half.Mxy.mean, relate.off.non.Mxy.mean, Thres)
    names(relate.return) <- c("Randomized_Fullssibs", "Randomized_Halfsibs", "Randomized_Non", "Thresholds")
		return(relate.return)
			
	}	

