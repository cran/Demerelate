rxy <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
	{ 
  
  ref.pop <- .subset2(ref.pop,1)
  n.ref.pop <- names(ref.pop)
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
  
  # Calculations after Oliehoek et al 2006 (15) + (16) 
  # Share.RE for colnames/pop2
  pxma <- .subset2(ref.pop,which(n.ref.pop==ai))
  pymb <- .subset2(ref.pop,which(n.ref.pop==aj))
  het <-  ai==aj
  
  share.RE <- (((ai==bi)+(ai==bj)+(aj==bi)+(aj==bj))*0.5-pxma-pymb)/(1+het-pxma-pymb)
  
  # Share.RAT for pop2
  # het == 1 if Alleles identical otherwise het == 0
  pxma <- .subset2(ref.pop,which(n.ref.pop==bi))
  pymb <- .subset2(ref.pop,which(n.ref.pop==bj))
  het <-  bi==bj
    
  share.RAT <- (((bi==ai)+(bi==aj)+(bj==ai)+(bj==aj))*0.5-pxma-pymb)/(1+het-pxma-pymb)

  
  share <- (share.RE+share.RAT)/2
  if (length(share)>0){return(share)}else{return(NA)}
  
	}
