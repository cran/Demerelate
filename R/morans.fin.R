### P. Kraemer 13.5.2015
# Updated 10.12.2015 

# Morans I is calculated based on Hardy1999 Equation (4), localities is considered as different alleles of i+j i.e individuals, this could be any other deme here as well.

morans.fin <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
 # ref.pop <- ref.pop[[allele.column]]  
  
  pop <- rbind(pop1,pop2)
  ref.pop <- table(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))/length(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))
 
  ## Do I need to fill p instead of ref.pop?? - It makes no difference, ref.pop is better..
  ## The problem is every ref.pop frequency adds to (-ref.pop[[x]])^2 and adds 0 to morans.w
  ## Solution: ref.pop is calculated from empirical data even for randomized data this is the way to go since the variance in allele freq in individuals needs to be empirical and cannot be tagen from a better panmicitc population.
  ## 
  
  r.return <- sum(sapply(seq(1:length(ref.pop)),function(x){
    (mean(names(ref.pop[x])==pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])*
    (mean(names(ref.pop[x])==pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])
                                                     }))
  
  return(r.return)
  
}