loiselle.w <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  
  r.return <- sum(ref.pop[[allele.column]]*(1-ref.pop[[allele.column]]))
  
    
  return(r.return)
  

}