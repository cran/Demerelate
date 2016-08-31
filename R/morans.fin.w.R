morans.fin.w <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  ## Test for p needed 05072016
  ## changed for speed issues 03082016
  ## in version 0.9 65s for demerelref
  ## now 12.181s

  pop <- rbind(pop1,pop2)
  ref.pop <- table(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))/length(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))

  var.p <- sapply(seq(1:length(ref.pop)),function(x){rbind(
    rowMeans(names(ref.pop[x])==pop1[,(allele.column*2+1):(allele.column*2+2)]),
    rowMeans(names(ref.pop[x])==pop2[,(allele.column*2+1):(allele.column*2+2)]))
  })
  
  var.p <-rowMeans((t(var.p)-colMeans(var.p))^2)
  
  return(sum(var.p))
  
}
