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
  ## N is calculated from population allele frequencies are based on
  ## artificial panmictic populations can only be used for the offspring randomization, allelefreq and var can only be calculated from each empirical population
  
 #n <- as.numeric(names(ref.pop))[allele.column]
 #ref.pop <- ref.pop[[allele.column]]
  
 pop <- rbind(pop1,pop2)
 ref.pop <- table(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))/length(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))
  
 if (identical(pop1,pop2)==TRUE){n<-nrow(pop1)}else{n<-nrow(pop)}
  
  
  r.return <- sum(sapply(seq(1:length(ref.pop)),function(x){
    (mean(names(ref.pop[x])==pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])*
    (mean(names(ref.pop[x])==pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])
                                                     }))

  var.p <- sapply(seq(1:length(ref.pop)),function(x){rbind(
    rowMeans(names(ref.pop[x])==pop1[,(allele.column*2+1):(allele.column*2+2)]),
    rowMeans(names(ref.pop[x])==pop2[,(allele.column*2+1):(allele.column*2+2)]))
  })
  
  var.p <-rowMeans((t(var.p)-colMeans(var.p))^2)
 
  return(r.return+sum(var.p)/(n-1))
  
}
