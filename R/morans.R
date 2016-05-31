### P. Kraemer 19.5.2016
# Updated 10.12.2015 

# Calculations are still not correct
# Morans I is calculated based on Hardy1999 Equation (4), localities is considered as different alleles of i+j i.e individuals, this could be any other deme here as well.

morans <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  
  # n is it used from ref.pop or the sum of pop1 and pop2 ?? 
  # problematic if randomized pops are used i.e. offspring
  ## var is caluclated correctly
  
  n <- as.numeric(names(ref.pop))[allele.column]
  #n <- length(rbind(pop1,pop2)[,1])
  ref.pop <- ref.pop[[allele.column]]
  
  
  r.return <- sum(sapply(seq(1:length(ref.pop)),function(x){
    (mean(names(ref.pop[x])==pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])*
    (mean(names(ref.pop[x])==pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])
                                                     }))

  var.p <- sapply(seq(1:length(ref.pop)),function(x){rbind(
     apply(names(ref.pop[x])==pop1[,(allele.column*2+1):(allele.column*2+2)],1,mean),
     apply(names(ref.pop[x])==pop2[,(allele.column*2+1):(allele.column*2+2)],1,mean))
                                                    })

  var.p <-apply((t(apply(var.p,1,'-',apply(var.p,2,mean))))^2,2,mean)
 
  return(r.return+sum(var.p)/(n-1))
  
}
