morans.w <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  
  #ref.pop <- ref.pop[[allele.column]]
  
  pop <- rbind(pop1,pop2)
  ref.pop <- table(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))/length(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))
  
  var.p <- sapply(seq(1:length(ref.pop)),function(x){rbind(
    rowMeans(names(ref.pop[x])==pop1[,(allele.column*2+1):(allele.column*2+2)]),
    rowMeans(names(ref.pop[x])==pop2[,(allele.column*2+1):(allele.column*2+2)]))
  })
  
  var.p <-rowMeans((t(var.p)-colMeans(var.p))^2)
  
  return(sum(var.p))
  
}