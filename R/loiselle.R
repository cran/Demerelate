loiselle <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  
  n <- as.numeric(names(ref.pop))
  ref.pop <- .subset2(ref.pop,1)
  p <- ref.pop
  n.ref.pop <- names(p)
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  ai<-n.ref.pop==ai
  aj<-n.ref.pop==aj
  bi<-n.ref.pop==bi
  bj<-n.ref.pop==bj
  
  r.return <- sum((ai/2+aj/2-p)*(bi/2+bj/2-p))+sum(p*(1-p)/(2*n-1))
  
  return(r.return)
  

}
