### P. Kraemer 13.5.2015
# Updated 10.12.2015 

# Morans I is calculated based on Hardy1999 Equation (4), localities is considered as different alleles of i+j i.e individuals, this could be any other deme here as well.

morans.fin <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1

  n.ref.pop <- names(ref.pop)
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  ai<-n.ref.pop==ai
  aj<-n.ref.pop==aj
  bi<-n.ref.pop==bi
  bj<-n.ref.pop==bj
  
  r.return<-sum((ai/2+aj/2-ref.pop)*(bi/2+bj/2-ref.pop))
 
  ## Do I need to fill p instead of ref.pop?? - It makes no difference, ref.pop is better..
  ## The problem is every ref.pop frequency adds to (-ref.pop[[x]])^2 and adds 0 to morans.w
  ## Solution: ref.pop is calculated from empirical data even for randomized data this is the way to go since the variance in allele freq in individuals needs to be empirical and cannot be tagen from a better panmicitc population.
  ## 

  
  return(r.return)
  
}
