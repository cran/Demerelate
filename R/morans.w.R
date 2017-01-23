morans.w <- function(pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  ## Test for p needed 05072016
  ## changed for speed issues 03082016
  ## in version 0.9 65s for demerelref
  ## now 12.181s
  ## pop nicht gleich refpop wie in emp.calc
  pop <- rbind(pop1,pop2)
  
  ri <- .subset2(pop1,allele.column*2+1)
  rj <- .subset2(pop2,allele.column*2+2)

  var.p <- sapply(seq(1:length(ref.pop)),function(x){rowMeans(names(ref.pop[x])==data.frame(ri,rj))})
  
  var.p <-rowMeans((t(var.p)-colMeans(var.p))^2)
  
  return(sum(var.p))
  
}
