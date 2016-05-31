# Updated 08.10.2015 Similarity is given is not correct
# Former version
# Estimate compare Li and Horvitz 1953
# Compare calculation in Wang 2002 equation (15). It is based on Li et al 1993 and Lynch 1988 and is denoted as LL in Coancestry

# Similarity is given by:
# 1,1 - 1,2 : 0.25
# 1,2 - 1,2 : 0.25
# 1,2 - 1,3 : 0.125
# 1,1 - 1,1 : 1
# Rest      : 0

Bxy <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
  
{  
  
  # Estimate compare Li and Horvitz 1953
  
  # Similarity is given by:
  # 1,1 - 1,2 : 0.5
  # 1,2 - 1,2 : 0.5
  # 1,2 - 1,3 : 0.25
  # 1,1 - 1,1 : 1
  # Rest      : 0
  
  sum(
    (pop1[data[row,1],(allele.column*2+1)]==pop2[data[row,2],(allele.column*2+1)]),
    (pop1[data[row,1],(allele.column*2+1)]==pop2[data[row,2],(allele.column*2+2)]),
    (pop1[data[row,1],(allele.column*2+2)]==pop2[data[row,2],(allele.column*2+1)]),
    (pop1[data[row,1],(allele.column*2+2)]==pop2[data[row,2],(allele.column*2+2)])
  )/4	
  
  
  
}
