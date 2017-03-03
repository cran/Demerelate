morans.w <- function(pop1, pop2, allele.column, ref.pop=NA)
{

  pop <- rbind(pop1,pop2)
  
  ri <- .subset2(pop1,allele.column*2+1)
  rj <- .subset2(pop2,allele.column*2+2)

  var.p <- sapply(seq(1:length(ref.pop)),function(x){rowMeans(names(ref.pop[x])==data.frame(ri,rj))})
  
  var.p <-rowMeans((t(var.p)-colMeans(var.p))^2)
  
  return(sum(var.p))
  
}
