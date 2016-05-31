# Written by P. Kraemer
# Last modified 07.10.2015

Emp.calc <- function(tab.pop.pop, value="NA", ref.pop="NA")
{
   
  if (is(ref.pop)[1]!="data.frame") {ref.pop <- tab.pop.pop}
  
  
single.pop <- function(tab.pop.pop, value, ref.pop)    
  {
  
    number.loci <- (ncol(tab.pop.pop)-2)/2
    
    # List of all reference allele frequencies
    x <- seq(3,length(ref.pop[1,]),2)
    pm <- lapply(x,function(x){table(c(ref.pop[,x],ref.pop[,x+1]))/length(c(ref.pop[,x],ref.pop[,x+1]))})
    names(pm) <- lapply(x,function(x){sum(complete.cases(cbind(ref.pop[,x],ref.pop[,x+1])))})    
    
    empirical.share.ls <- vector("list",number.loci)
    if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="wang" | value=="wang.fin" | value=="ritland") {empirical.share.ls.w <- vector("list",number.loci)}
    
    # Calculation of value for each locus in population tab.pop.pop
    for (i in 1:number.loci)
    {
      # Empirical weights calculated for each locus
      if (value=="lxy") {empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,FALSE, value="lxy.w", pm)}
      if (value=="ritland") {empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,FALSE, value="ritland.w", pm)}
      if (value=="loiselle") {empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,FALSE, value="loiselle.w", pm)}
      if (value=="morans") {empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,FALSE, value="morans.w", pm)}
      if (value=="morans.fin") {empirical.share.ls.w[[i]] <- allele.sharing(tab.pop.pop,tab.pop.pop,i,FALSE, value="morans.fin.w", pm)}
      if (value=="wang") {empirical.share.ls.w[[i]] <- wang.w(allele.column=i, ref.pop=pm)}
      if (value=="wang.fin") {empirical.share.ls.w[[i]] <- wang.fin.w(allele.column=i, ref.pop=pm)}
      
      # Empirical share calculated for each locus
      message(paste("---","Calculations for empirical values are performed for Locus",i,"----",Sys.time()),"\n")
      empirical.share.ls[[i]] <- allele.sharing(tab.pop.pop, tab.pop.pop, i, FALSE, value, pm)
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
      empirical.share.ls <- Reduce("+",empirical.share.ls)#/length(empirical.share.ls)
      empirical.share.ls <- empirical.share.ls*(1/(sum(1/u)))
      empirical.share.ls.w <- Reduce("+",empirical.share.ls.w)#/length(empirical.share.ls.w)
      empirical.share.ls.w <- empirical.share.ls.w*(1/(sum(1/u)))
         
      # compose average
      
     empirical.list <- rowMeans(do.call("rbind",lapply(seq(1:length(empirical.share.ls[,1])), function(x){wang.compose(as=empirical.share.ls.w,Ps=empirical.share.ls[x,])})))
    } 
    
    if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li")
    {
    empirical.list <- do.call("cbind",empirical.share.ls)	
    empirical.list <- matrix(rowMeans(empirical.list,na.rm=TRUE),nrow(empirical.share.ls[[1]]))
    row.names(empirical.list) <- row.names(empirical.share.ls[[1]])
    colnames(empirical.list) <- colnames(empirical.share.ls[[1]])
    }
    
    if (value=="loiselle" | value=="morans.fin" | value=="morans" | value=="ritland")
    {
      empirical.list <- do.call("cbind",empirical.share.ls)
      empirical.list.w <- do.call("cbind",empirical.share.ls.w) 
      empirical.list <- matrix(rowSums(empirical.list,na.rm=TRUE),nrow(empirical.share.ls[[1]]))
      empirical.list.w <- matrix(rowSums(empirical.list.w,na.rm=TRUE),nrow(empirical.share.ls.w[[1]]))
      empirical.list <- empirical.list/empirical.list.w
      row.names(empirical.list) <- row.names(empirical.share.ls[[1]])
      colnames(empirical.list) <- colnames(empirical.share.ls[[1]])
    }
    
    if (value=="lxy")
    {
      # RE
      empirical.list <- do.call("cbind",sapply(empirical.share.ls,function(x)(x[[1]])))
      empirical.list.w <- do.call("cbind",sapply(empirical.share.ls.w,function(x)(x[[1]])))
      empirical.list <- matrix(rowSums(empirical.list,na.rm=TRUE),nrow(empirical.share.ls[[1]][[1]]))
      empirical.list.w <- matrix(rowSums(empirical.list.w,na.rm=TRUE),nrow(empirical.share.ls.w[[1]][[1]]))
      empirical.list <- empirical.list/empirical.list.w
      
      # RAT
      empirical.list.rat <- do.call("cbind",sapply(empirical.share.ls,function(x)(x[[2]])))
      empirical.list.w.rat <- do.call("cbind",sapply(empirical.share.ls.w,function(x)(x[[2]]))) 
      empirical.list.rat <- matrix(rowSums(empirical.list.rat,na.rm=TRUE),nrow(empirical.share.ls[[1]][[2]]))
      empirical.list.w.rat <- matrix(rowSums(empirical.list.w.rat,na.rm=TRUE),nrow(empirical.share.ls.w[[1]][[2]]))
      empirical.list.rat <- empirical.list.rat/empirical.list.w.rat
      
      # AVERAGE
      empirical.list <- (empirical.list+empirical.list.rat)/2
      row.names(empirical.list) <- row.names(empirical.share.ls[[1]][[1]])
      colnames(empirical.list) <- colnames(empirical.share.ls[[1]][[1]])
    }

    return(empirical.list)
    
  }

    if (value=="NA")
      {
      message(paste("---","Relatedness calculations are performed for Bxy","----", Sys.time()),"\n")
      empirical.list.out <- single.pop(tab.pop.pop, "Bxy", ref.pop)
      message(paste("---","Relatedness calculations are performed for Sxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "Sxy", ref.pop))
      message(paste("---","Relatedness calculations are performed for Mxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "Mxy", ref.pop))
      message(paste("---","Relatedness calculations are performed for Li","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "Li", ref.pop))
      message(paste("---","Relatedness calculations are performed for rxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "rxy", ref.pop))
      message(paste("---","Relatedness calculations are performed for lxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "lxy", ref.pop))
      message(paste("---","Relatedness calculations are performed for loiselle","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "loiselle", ref.pop))
      message(paste("---","Relatedness calculations are performed for Wang (Finite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "wang.fin", ref.pop))
      message(paste("---","Relatedness calculations are performed for Wang (Infinite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "wang", ref.pop))
      message(paste("---","Relatedness calculations are performed for Ritland","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "ritland", ref.pop))
      message(paste("---","Relatedness calculations are performed for Morans I (Finite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "morans.fin", ref.pop))
      message(paste("---","Relatedness calculations are performed for Morans I (Infinite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "morans", ref.pop))
      
      colnames(empirical.list.out) <- c("Bxy", "Sxy", "Mxy", "Li", "rxy", "lxy", "loiselle", "wang.fin", "wang", "ritland", "morans_I.fin", "morans_I")
      }
    else
      {
      empirical.list.out <- single.pop(tab.pop.pop, value, ref.pop)      
      }
    
 return(empirical.list.out)
}
