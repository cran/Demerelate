lxy.w <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{ 
  
  n.ref.pop <- names(ref.pop)
  re <- .subset2(data,1)[row]
  rat <- .subset2(data,2)[row]
  a <- allele.column*2+1
  
  ai <- .subset2(pop1,a)[re]
  aj <- .subset2(pop1,a+1)[re]
  bi <- .subset2(pop2,a)[rat]
  bj <- .subset2(pop2,a+1)[rat]
   
  # Share.RE for colnames/pop2
  pxma <- .subset2(ref.pop,which(n.ref.pop==ai))
  pymb <- .subset2(ref.pop,which(n.ref.pop==aj))
  het <-  ai==aj
  
  # weighted for Lynch 1999 equ 7a
  
  weight.RE <- ((1+het)*(pxma+pymb)-4*pxma*pymb)/(2*pxma*pymb)
  
  # Share.RAT for pop2
  # het == 1 if Alleles identical otherwise het == 0
  pxma <- .subset2(ref.pop,which(n.ref.pop==bi))
  pymb <- .subset2(ref.pop,which(n.ref.pop==bj))
  het <-  bi==bj
  
  # weighted for Lynch 1999 equ 7a
  
  weight.RAT <- ((1+het)*(pxma+pymb)-4*pxma*pymb)/(2*pxma*pymb)
  
  return(c(weight.RE,weight.RAT))
  

}
