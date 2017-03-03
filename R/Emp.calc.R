Emp.calc <- function(tab.pop.pop, value="NA", ref.pop="NA")
{
  
  if (is(ref.pop)[1]=="vector" | is(ref.pop)[1]=="list") {pm <- ref.pop} 
  if (is(ref.pop)[1]=="character") {ref.pop <- tab.pop.pop}
  if (is(ref.pop)[1]=="data.frame") {
    x <- seq(3,length(ref.pop[1,]),2)
    pm <- lapply(x,function(x){table(c(ref.pop[,x],ref.pop[,x+1]))/length(c(ref.pop[,x],ref.pop[,x+1]))})
    names(pm) <- (sapply(x,function(x){sum(complete.cases(cbind(ref.pop[,x],ref.pop[,x+1])))}))
    pm <- c(pm,list(sapply(x,function(x){sum(complete.cases(cbind(ref.pop[,x],ref.pop[,x+1])))})))
    
  }
  
  if (any(lapply(pm,length)==2,TRUE)){warning("Careful, bi-allelic markers are detected! 
  Especially, rxy and ritland estimator are not defined when bi-allelic estimates are used with allele frequencies of 0.5.
  You should consider removing bi-allelics which tend to have very evenly distributed alleles or swich to another estimator. 
  Be careful even if allele frequencies are not perfectly 0.5, during randomizations problems may occur due to producing randomly such populations.")}
  
single.pop <- function(tab.pop.pop, value, ref.pop)    
  {
  
    number.loci <- (ncol(tab.pop.pop)-2)/2
    
    # Preparing vectors/dfs for applying
    data <- expand.grid(seq(1:length(tab.pop.pop[,1])),seq(1:length(tab.pop.pop[,1])))
    data <- data[as.numeric(data[,1]) <= as.numeric(data[,2]),]
    data <- data[data[,1]!=data[,2],]
    
    empirical.share.ls <- vector("list",number.loci)
    if (value=="lxy" | value=="loiselle" | value=="morans" | value=="morans.fin" | value=="wang" | value=="wang.fin" | value=="ritland") {empirical.share.ls.w <- vector("list",number.loci)}
    
    loci <- seq(1:number.loci)
    
    # Calculate each locus in list
    empirical.share.ls <- lapply(loci,function(x){allele.sharing(tab.pop.pop, tab.pop.pop, x, data, value, pm[x])})
    if (value=="lxy") {empirical.share.ls.w <- lapply(loci,function(x){allele.sharing(tab.pop.pop, tab.pop.pop, x, data, value="lxy.w", pm[[x]])})}
    if (value=="ritland") {empirical.share.ls.w <- as.list(sapply(pm[-length(pm)],function(x){length(x)-1}))}
    if (value=="loiselle")  {empirical.share.ls.w<- as.list(sapply(pm[-length(pm)],function(x){sum(x*(1-x))}))}
    if (value=="morans" | value=="morans.fin") {empirical.share.ls.w <- lapply(loci,function(x){allele.sharing(tab.pop.pop, tab.pop.pop, x, data, value="morans.w", pm[[x]])})}
    if (value=="wang") {empirical.share.ls.w <- lapply(loci,function(x){wang.w(allele.column=x, ref.pop=pm)})}
    if (value=="wang.fin") {empirical.share.ls.w <- lapply(loci,function(x){wang.fin.w(allele.column=x, ref.pop=pm[[x]])})}
  
    
    if (value=="wang" | value=="wang.fin")
    {
      u <- unlist(lapply(seq(1:length(empirical.share.ls.w)),function(x){u<-empirical.share.ls.w[[x]][7]}))
      empirical.share.ls <- lapply(seq(1:length(empirical.share.ls)),function(x){empirical.share.ls[[x]] * 1/u[x]})
      empirical.share.ls <- Reduce("+",empirical.share.ls)
      empirical.share.ls <- empirical.share.ls*(1/sum(1/u))
      empirical.share.ls.w <- Reduce("+",empirical.share.ls.w)
      empirical.share.ls.w <- empirical.share.ls.w*(1/sum(1/u))
         
      # compose average
      
     empirical.list <- rowMeans(do.call("rbind",lapply(seq(1:length(empirical.share.ls[,1])), function(x){wang.compose(as=empirical.share.ls.w,Ps=empirical.share.ls[x,])})))
    } 
    
    if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li")
    {
    empirical.list <- do.call("cbind",empirical.share.ls)	
    empirical.list <- rowMeans(empirical.list,na.rm=TRUE)
    }
    
    if (value=="loiselle" | value=="morans.fin" | value=="morans" | value=="ritland")
    {
      empirical.list <- do.call("cbind",empirical.share.ls)
      empirical.list.w <- do.call("cbind",empirical.share.ls.w) 
      empirical.list <- rowSums(empirical.list,na.rm=TRUE)
      empirical.list.w <- rowSums(empirical.list.w,na.rm=TRUE)
      empirical.list <- empirical.list/empirical.list.w
    }
    
    if (value=="lxy")
    {
      # RE
      empirical.list <- do.call("cbind", lapply(empirical.share.ls, function(x)(.subset2(x,1))))
      empirical.list.w <- do.call("cbind", lapply(empirical.share.ls.w, function(x)(.subset2(x,1))))
      empirical.list <- rowSums(empirical.list,na.rm=TRUE)
      empirical.list.w <- rowSums(empirical.list.w,na.rm=TRUE)
      empirical.list <- empirical.list/empirical.list.w
      
      # RAT
      empirical.list.rat <- do.call("cbind", lapply(empirical.share.ls, function(x)(.subset2(x,2))))
      empirical.list.w.rat <- do.call("cbind", lapply(empirical.share.ls.w, function(x)(.subset2(x,2))))
      empirical.list.rat <- rowSums(empirical.list.rat,na.rm=TRUE)
      empirical.list.w.rat <- rowSums(empirical.list.w.rat,na.rm=TRUE)
      empirical.list.rat <- empirical.list.rat/empirical.list.w.rat
      
      # AVERAGE
      empirical.list <- (empirical.list+empirical.list.rat)/2
    }
    

    names(empirical.list) <- apply(expand.grid(tab.pop.pop[,1],tab.pop.pop[,1]),1,paste,collapse="_")[as.numeric(row.names(data))]

    
    return(empirical.list)
    
  }

    if (value=="NA")
      {
      message(paste("---","Relatedness calculations are performed for Bxy","----", Sys.time()),"\n")
      empirical.list.out <- single.pop(tab.pop.pop, "Bxy", )
      message(paste("---","Relatedness calculations are performed for Sxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "Sxy", pm))
      message(paste("---","Relatedness calculations are performed for Mxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "Mxy", pm))
      message(paste("---","Relatedness calculations are performed for Li","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "Li", pm))
      message(paste("---","Relatedness calculations are performed for rxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "rxy", pm))
      message(paste("---","Relatedness calculations are performed for lxy","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "lxy", pm))
      message(paste("---","Relatedness calculations are performed for loiselle","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "loiselle", pm))
      message(paste("---","Relatedness calculations are performed for Wang (Finite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "wang.fin", pm))
      message(paste("---","Relatedness calculations are performed for Wang (Infinite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "wang", pm))
      message(paste("---","Relatedness calculations are performed for Ritland","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "ritland", pm))
      message(paste("---","Relatedness calculations are performed for Morans I (Finite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "morans.fin", pm))
      message(paste("---","Relatedness calculations are performed for Morans I (Infinite)","----", Sys.time()),"\n")
      empirical.list.out <- cbind(empirical.list.out,single.pop(tab.pop.pop, "morans", pm))
      
      colnames(empirical.list.out) <- c("Bxy", "Sxy", "Mxy", "Li", "rxy", "lxy", "loiselle", "wang.fin", "wang", "ritland", "morans_I.fin", "morans_I")
      
      }
    else
      {
      empirical.list.out <- single.pop(tab.pop.pop, value, pm)      
      }
    
 return(empirical.list.out)
}
