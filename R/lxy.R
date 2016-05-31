lxy <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
	{ 

  ref.pop <- ref.pop[[allele.column]]

  # Calculations after Lynch and Ritland 1999 (5a) 
  # Share.RE for colnames/pop2
  pxma <- as.numeric(ref.pop[which(names(ref.pop)==pop1[data[row,1],(allele.column*2+1)])])
  pymb <- as.numeric(ref.pop[which(names(ref.pop)==pop1[data[row,1],(allele.column*2+2)])])
  het <-  as.numeric(pop1[data[row,1],(allele.column*2+1)]==pop1[data[row,1],(allele.column*2+2)])
  
  share.RE <-
    (
      sum(
        as.numeric(pop1[data[row,1],(allele.column*2+1)]==pop2[data[row,2],(allele.column*2+1)]),
        as.numeric(pop1[data[row,1],(allele.column*2+1)]==pop2[data[row,2],(allele.column*2+2)]))*pymb
          +
        sum(
          as.numeric(pop1[data[row,1],(allele.column*2+2)]==pop2[data[row,2],(allele.column*2+1)]),
          as.numeric(pop1[data[row,1],(allele.column*2+2)]==pop2[data[row,2],(allele.column*2+2)]))*pxma-4*pxma*pymb    
    )/((1+het)*(pxma+pymb)-4*pxma*pymb)
  
  # weighted for Lynch 1999 equ 7a
  
 share.RE <- share.RE*((1+het)*(pxma+pymb)-4*pxma*pymb)/(2*pxma*pymb)
  
  # Share.RAT for pop2
  # het == 1 if Alleles identical otherwise het == 0
  pxma <- as.numeric(ref.pop[which(names(ref.pop)==pop2[data[row,2],(allele.column*2+1)])])
  pymb <- as.numeric(ref.pop[which(names(ref.pop)==pop2[data[row,2],(allele.column*2+2)])])
  het <-  as.numeric(pop2[data[row,2],(allele.column*2+1)]==pop2[data[row,2],(allele.column*2+2)])
  
  share.RAT <- 
    (
      sum(as.numeric(pop2[data[row,2],(allele.column*2+1)]==pop1[data[row,1],(allele.column*2+1)]),
          as.numeric(pop2[data[row,2],(allele.column*2+1)]==pop1[data[row,1],(allele.column*2+2)]))*pymb
            +
      sum(as.numeric(pop2[data[row,2],(allele.column*2+2)]==pop1[data[row,1],(allele.column*2+1)]),
          as.numeric(pop2[data[row,2],(allele.column*2+2)]==pop1[data[row,1],(allele.column*2+2)]))*pxma-4*pxma*pymb
      
    )/((1+het)*(pxma+pymb)-4*pxma*pymb)

  # weighted for Lynch 1999 equ 7a
  
  share.RAT <- share.RAT*((1+het)*(pxma+pymb)-4*pxma*pymb)/(2*pxma*pymb)
  
  share <- c(share.RE,share.RAT)	
  if (length(share)>0){return(share)}else{return(NA)}

	}
