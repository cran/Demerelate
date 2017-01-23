# Estimate compare Li and Horvitz 1953
# Compare calculation in Wang 2002 equation (15). It is based on Li et al 1993 and Lynch 1988 and is denoted as LL in Coancestry

Bxy <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
  
{  
  a <- allele.column*2+1
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  
  # Estimate compare Li and Horvitz 1953
  
  # Similarity is given by:
  # 1,1 - 1,2 : 0.5
  # 1,2 - 1,2 : 0.5
  # 1,2 - 1,3 : 0.25
  # 1,1 - 1,1 : 1
  # Rest      : 0
  
  ((ai==bi)+(ai==bj)+(aj==bi)+(aj==bj))/4
  
}
