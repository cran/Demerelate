### P. Kraemer 13.5.2015
# Updated 10.12.2015 

# Calculations are still not correct
# Morans I is calculated based on Hardy1999 Equation (4), localities is considered as different alleles of i+j i.e individuals, this could be any other deme here as well.

morans.fin <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  ref.pop <- ref.pop[[allele.column]]
 
  r.return <- sum(sapply(seq(1:length(ref.pop)),function(x){
    (mean(names(ref.pop[x])==pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])*
    (mean(names(ref.pop[x])==pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])
                                                     }))
  

  
  return(r.return)
  
}