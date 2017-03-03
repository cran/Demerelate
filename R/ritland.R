ritland <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  
  ref.pop <-.subset2(ref.pop,1)
  p <- ref.pop
  n.p <- names(p)
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  
  # Relative frequency for specific alleles in ref.pop
  pxma <- p[which(n.p==ai)]
  pxmb <- p[which(n.p==aj)]
  pyma <- p[which(n.p==bi)]
  pymb <- p[which(n.p==bj)]
  
  p<-c(pxma,pxmb,pyma,pymb)[unique(names(c(pxma,pxmb,pyma,pymb)))]
  n.p <- names(p)
  
  ai<-n.p==ai
  aj<-n.p==aj
  bi<-n.p==bi
  bj<-n.p==bj
  
  return(2*((sum((ai*bi)/p+(ai*bj)/p+(aj*bi)/p+(aj*bj)/p)/4)-1))

  ## 2* is calculated due to Lynch & Ritland 1999 in order to be compareable to loiselle
}
