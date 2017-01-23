Li <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
  
{
  ## Calculations are not equivalent to Spagedi, since other bias corrections are applied. For single locus estimates 
  # It is based on Li et al 1993 and Lynch 1988 and is denoted as LL in Coancestry
  
  # Similarity is given by:
  # 1,1 - 1,2 : 0.75
  # 1,2 - 1,2 : 1
  # 1,2 - 1,3 : 0.5
  # 1,1 - 1,1 : 1
  # Rest      : 0
  
  
  N <- as.numeric(names(ref.pop))
  p <- .subset2(ref.pop,1)
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  aij <- c(ai,aj)
  bij <- c(bi,bj)
  
  
  a2 <- (N*(sum(p^2))-1)/(N-1)
  
  # eq. 13  
  a3 <- (N^2*sum(p^3)-3*(N-1)*a2-1)/((N-1)*(N-2))
  
  u <- 2*a2-a3
    
  # Equation (9) in Li et al 1993
  
 
 return(((
            sum(aij%in%bij,bij%in%aij)/4)-u)/(1-u)
        )
  
}
