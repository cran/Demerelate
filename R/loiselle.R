##
## 18.05.2016 loiselle is changed according to spagedi manual, equal results
## all alleles of ref.pop are considered for comparing frequencies i.e. indivdidual frequencies of alleles are in most cases 0

loiselle <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  
  n <- as.numeric(names(ref.pop)[allele.column])
  ref.pop <- ref.pop[[allele.column]]
  
  r.return <- sum(sapply(seq(1:length(ref.pop)),function(x){
      (mean(names(ref.pop[x])==pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])*
      (mean(names(ref.pop[x])==pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)])-ref.pop[x])
                                                     }))
  
  r.return <- r.return+sum(ref.pop*(1-ref.pop)/(2*n-1))
  
  return(r.return)
  

}