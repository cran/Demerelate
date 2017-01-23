Mxy <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
  
{
  
  
  # Estimate compare Blouin 1996
  
  # Similarity is given by:
  # 1,1 - 1,2 : 0.5
  # 1,2 - 1,2 : 1
  # 1,2 - 1,3 : 0.5
  # 1,1 - 1,1 : 1
  # Rest      : 0
  
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  sum(c(ai,aj)==c(bi,bj))/2
  
}