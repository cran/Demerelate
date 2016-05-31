Li <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
  
{
  ## Calculations are not equivalent to Spagedi, since other bias corrections are applied. For single locus estimates 
  # It is based on Li et al 1993 and Lynch 1988 and is denoted as LL in Coancestry
  
  # Similarity is given by:
  # 1,1 - 1,2 : 0.75
  # 1,2 - 1,2 : 1
  # 1,2 - 1,3 : 0.5
  # 1,1 - 1,1 : 1
  # Rest      : 0
  
  N <- as.numeric(names(ref.pop[allele.column]))
  p <- ref.pop[[allele.column]]
  
  
  a2 <- (N*(sum(p^2))-1)/(N-1)
  
  # eq. 13  
  a3 <- (N^2*sum(p^3)-3*(N-1)*a2-1)/((N-1)*(N-2))
  
  u <- 2*a2-a3
    
  # Equation (9) in Li et al 1993
  r.return <-((sum(
                  c(pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)]%in%pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)]),
                  c(pop2[data[row,2],(allele.column*2+1):(allele.column*2+2)]%in%pop1[data[row,1],(allele.column*2+1):(allele.column*2+2)]))
                  /4)-u)/(1-u) 
 
 return(r.return)
  
}
