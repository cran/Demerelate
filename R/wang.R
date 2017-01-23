### P. Kraemer 13.5.2015
# Updated 28.7.2015
## Basic functions from Wang 2002
## N in wang.w needs to be N from reference population ie. the poipulation ps are caulculated from
## This is the case for morans and loiselle as well

wang <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # P = vector of classes of observed similarities with 0/1 for each category
  #  one of 4 categories is possible:
  # 1: AiAi - AiAi or Aij - Aij; 1 -> P[1]
  # 2: AiAi - AiAj; 1 -> P[2]
  # 3: AiAj - AiAk; 1 -> P[3]
  # 4: all S=0    ; 1 -> P[4] 
  P <- rep(0,3)
  
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]


# Calculate the code for P  
if (sum((c(ai,aj)%in%c(bi,bj)),(c(bi,bj)%in%c(ai,aj)))==4){  P[1] <- 1 }
if (sum((c(ai,aj)%in%c(bi,bj)),(c(bi,bj)%in%c(ai,aj)))==3){  P[2] <- 1 }
if (sum((c(ai,aj)%in%c(bi,bj)),(c(bi,bj)%in%c(ai,aj)))==2){  P[3] <- 1 }

return(P)
}