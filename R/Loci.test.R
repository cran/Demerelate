Loci.test <- function(tab.pop, ref.pop="NA", object=FALSE, value=NA, bt=1000, file.output=FALSE)
  
{

  # Function tests on differences in mean relatedness based on number of loci used for calculation
  
  # Data are loaded for input
    if (object==FALSE) {
                          tab.pop <- input.txt(tab.pop, "pop")
                          if (ref.pop=="NA") {ref.pop <- tab.pop} else {ref.pop <- input.txt(ref.pop, "ref.pop")}
                        }
    
    if (object==TRUE) {
                          if (is(ref.pop)[1]!="data.frame") {ref.pop <- tab.pop}
    }
  
  if (is(ref.pop)[1]=="vector" | is(ref.pop)[1]=="list") {pm <- ref.pop} 
  if (is(ref.pop)[1]=="character") {ref.pop <- tab.pop}
  if (is(ref.pop)[1]=="data.frame") {
    x <- seq(3,length(ref.pop[1,]),2)
    pm <- lapply(x,function(x){table(c(ref.pop[,x],ref.pop[,x+1]))/length(c(ref.pop[,x],ref.pop[,x+1]))})
    names(pm) <- (sapply(x,function(x){sum(complete.cases(cbind(ref.pop[,x],ref.pop[,x+1])))}))
    pm <- c(pm,list(sapply(x,function(x){sum(complete.cases(cbind(ref.pop[,x],ref.pop[,x+1])))})))
  }
    
    
    # Preparing vectors/dfs for applying
    data <- expand.grid(seq(1:length(tab.pop[,1])),seq(1:length(tab.pop[,1])))
    data <- data[as.numeric(data[,1]) <= as.numeric(data[,2]),]
    data <- data[data[,1]!=data[,2],]
  
  
  loci <- seq(1,(length(tab.pop[1,])/2)-1)
  boots <- vector("list",bt)
  loc <- vector("list", length(loci))
  
  # Calculate each locus in list
  lis.boot <- lapply(loci,function(x){allele.sharing(tab.pop, tab.pop, x, data, value, pm[x])})
  if (value=="lxy") {lis.boot.w <- lapply(loci,function(x){allele.sharing(tab.pop, tab.pop, x, data, value="lxy.w", pm[[x]])})}
  if (value=="ritland") {lis.boot.w <- as.list(sapply(pm[-length(pm)],function(x){length(x)-1}))}
  if (value=="loiselle")  {lis.boot.w<- as.list(sapply(pm[-length(pm)],function(x){sum(x*(1-x))}))}
  if (value=="morans" | value=="morans.fin") {lis.boot.w <- lapply(loci,function(x){allele.sharing(tab.pop, tab.pop, x, data, value="morans.w", pm[[x]])})}
  if (value=="wang") {lis.boot.w <- lapply(loci,function(x){wang.w(allele.column=x, ref.pop=pm)})}
  if (value=="wang.fin") {lis.boot.w <- lapply(loci,function(x){wang.fin.w(allele.column=x, ref.pop=pm[[x]])})}
  

  for (i in 1:length(lis.boot)){names(lis.boot[[i]])<-as.character(i)}
  
    
  # Draw loci random for i loci from bt samples
  for (i in 1:length(loc))
  {
  # Draw random bt samples
  if (value=="lxy")
    {rsamp <- lapply(boots, function(x){sample(seq(1:length(lis.boot[[1]][[1]])),replace=T)})}
    else
    {if (value=="wang" | value=="wang.fin") {rsamp <- lapply(boots, function(x){sample(seq(1:nrow(lis.boot[[1]])),replace=T)})}
      else {rsamp <- lapply(boots, function(x){sample(seq(1:length(lis.boot[[1]])),replace=T)})}}
    
  for (z in 1:length(rsamp))
    {
    r.loc <- lapply(vector("list",i),function(x){sample(1:length(loc),1)})
    
    if (value=="lxy"){list.bt.loc.re <- lapply(seq(1:i), function(x){lis.boot[[r.loc[[x]]]][[1]][rsamp[[z]]]})
                      list.bt.loc.rat <- lapply(seq(1:i), function(x){lis.boot[[r.loc[[x]]]][[2]][rsamp[[z]]]})}
    else
    {if (value=="wang" | value=="wang.fin") {list.bt.loc <- lapply(seq(1:i), function(x){lis.boot[[r.loc[[x]]]][rsamp[[z]],]})}
      else {list.bt.loc <- lapply(seq(1:i), function(x){lis.boot[[r.loc[[x]]]][rsamp[[z]]]})}}
    
      if (value=="lxy")
      {
        list.bt.loc.w.re <- lapply(seq(1:i), function(x){lis.boot.w[[r.loc[[x]]]][[1]][rsamp[[z]]]})
        list.bt.loc.w.rat <- lapply(seq(1:i), function(x){lis.boot.w[[r.loc[[x]]]][[2]][rsamp[[z]]]})
      }
      
      if (value=="loiselle" | value=="ritland" | value=="morans" | value=="morans.fin")
      {
        list.bt.loc.w <- lapply(seq(1:i), function(x){rep(lis.boot.w[[r.loc[[x]]]],length(rsamp[[z]]))})
      }
    

      if (value=="wang" | value=="wang.fin")
      {
        list.bt.loc.w <- lapply(seq(1:i), function(x){lis.boot.w[[r.loc[[x]]]]})
      }
        
  
    # Mean over loci
    if (value=="loiselle" | value=="morans" | value=="morans.fin" | value=="ritland")
    {
    empirical.list <- do.call("cbind",list.bt.loc)
    empirical.list.w <- do.call("cbind",list.bt.loc.w) 
    empirical.list <- rowSums(empirical.list,na.rm=TRUE)
    empirical.list.w <- rowSums(empirical.list.w,na.rm=TRUE)
    empirical.list <- empirical.list/empirical.list.w
    row.names(empirical.list) <- row.names(list.bt.loc[[1]])
    colnames(empirical.list) <- colnames(list.bt.loc[[1]])
    }
  
    if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li")
    {
    empirical.list <- do.call("cbind",list.bt.loc)  
    empirical.list <- rowMeans(empirical.list,na.rm=TRUE)
    row.names(empirical.list) <- row.names(list.bt.loc[[1]])
    colnames(empirical.list) <- colnames(list.bt.loc[[1]])
    }
  
    if (value=="wang" | value=="wang.fin")
    {
    u <- unlist(lapply(seq(1:length(list.bt.loc.w)),function(x){u<-list.bt.loc.w[[x]][7]}))
    list.bt.loc <- lapply(seq(1:length(list.bt.loc)),function(x){list.bt.loc[[x]] * 1/u[x]})
    list.bt.loc <- Reduce("+",list.bt.loc)
    list.bt.loc <- list.bt.loc*(1/sum(1/u))
    list.bt.loc.w <- Reduce("+",list.bt.loc.w)
    list.bt.loc.w <- list.bt.loc.w*(1/sum(1/u))
    
    # compose average
    
    empirical.list <- rowMeans(do.call("rbind",lapply(seq(1:length(list.bt.loc[,1])), function(x){wang.compose(as=list.bt.loc.w,Ps=list.bt.loc[x,])})))
    } 
    
    if (value=="lxy")
    {
    # RE
    empirical.list <- do.call("cbind",list.bt.loc.re)
    empirical.list.w <- do.call("cbind",list.bt.loc.w.re)
    empirical.list <- rowSums(empirical.list,na.rm=TRUE)
    empirical.list.w <- rowSums(empirical.list.w,na.rm=TRUE)
    empirical.list <- empirical.list/empirical.list.w
    
    # RAT
    empirical.list.rat <- do.call("cbind",list.bt.loc.rat)
    empirical.list.w.rat <- do.call("cbind",list.bt.loc.w.rat) 
    empirical.list.rat <- rowSums(empirical.list.rat,na.rm=TRUE)
    empirical.list.w.rat <- rowSums(empirical.list.w.rat,na.rm=TRUE)
    empirical.list.rat <- empirical.list.rat/empirical.list.w.rat
    
    # AVERAGE
    empirical.list <- (empirical.list+empirical.list.rat)/2
    row.names(empirical.list) <- row.names(list.bt.loc.re[[1]])
    colnames(empirical.list) <- colnames(list.bt.loc.re[[1]])
    }
    
   # names(empirical.list) <- apply(expand.grid(tab.pop[,1],tab.pop[,1]),1,paste,collapse="_")[as.numeric(row.names(data))]
    
    boots[[z]] <- empirical.list
    
  }
   
  loc[[i]] <- dist(lapply(boots, function(x){
    mean(x, na.rm = T)}))
  }
  
  
if (file.output==TRUE)
{  
  pdf(paste("Loci.test",value,Sys.Date(),".pdf"))
  
  y <- unlist(lapply(loc,mean))
  errbar(x=c(1:length(loc)),y,yplus=y+unlist(lapply(loc,sd)),yminus=y-unlist(lapply(loc,sd)),xlab="Number of Random Loci", ylab=paste("Mean Difference in Relatedness [",value,"]"))
  lines(unlist(lapply(loc,mean))~c(1:length(loc)))
  
  dev.off()
  
}
  return(loc)
  
}
