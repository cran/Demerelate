lxy.w <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{ 
  ref.pop <- ref.pop[[allele.column]]
   
  # Share.RE for colnames/pop2
  pxma <- as.numeric(ref.pop[which(names(ref.pop)==pop1[data[row,1],(allele.column*2+1)])])
  pymb <- as.numeric(ref.pop[which(names(ref.pop)==pop1[data[row,1],(allele.column*2+2)])])
  het <-  as.numeric(pop1[data[row,1],(allele.column*2+1)]==pop1[data[row,1],(allele.column*2+2)])
  
  # weighted for Lynch 1999 equ 7a
  
  weight.RE <- ((1+het)*(pxma+pymb)-4*pxma*pymb)/(2*pxma*pymb)
  
  # Share.RAT for pop2
  # het == 1 if Alleles identical otherwise het == 0
  pxma <- as.numeric(ref.pop[which(names(ref.pop)==pop2[data[row,2],(allele.column*2+1)])])
  pymb <- as.numeric(ref.pop[which(names(ref.pop)==pop2[data[row,2],(allele.column*2+2)])])
  het <-  as.numeric(pop2[data[row,2],(allele.column*2+1)]==pop2[data[row,2],(allele.column*2+2)])
  
  # weighted for Lynch 1999 equ 7a
  
  weight.RAT <- ((1+het)*(pxma+pymb)-4*pxma*pymb)/(2*pxma*pymb)
  
  weight <- c(weight.RE,weight.RAT)		
  if (length(weight)>0){return(weight)}else{return(NA)}

}
