Mxy <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
  
{
  
  # Estimate compare Blouin 1996
  
  # Similarity is given by:
  # 1,1 - 1,2 : 0.5
  # 1,2 - 1,2 : 1
  # 1,2 - 1,3 : 0.5
  # 1,1 - 1,1 : 1
  # Rest      : 0
  
  sum(c(
    (pop1[data[row,1],(allele.column*2+1)]),
    (pop1[data[row,1],(allele.column*2+2)])
  )== 
    c(
      (pop2[data[row,2],(allele.column*2+1)]),
      (pop2[data[row,2],(allele.column*2+2)])
    ))/2
  
}