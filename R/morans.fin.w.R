morans.fin.w <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  # Morans I from Spagedi manual
  # Variance is calculated from pop1 over individuals, this should be changed if onlypairs==TRUE

  ref.pop <- ref.pop[[allele.column]]
  
  var.p <- sapply(seq(1:length(ref.pop)),function(x){rbind(
    apply(names(ref.pop[x])==pop1[,(allele.column*2+1):(allele.column*2+2)],1,mean),
    apply(names(ref.pop[x])==pop2[,(allele.column*2+1):(allele.column*2+2)],1,mean))
  })
  
  var.p <-apply((t(apply(var.p,1,'-',apply(var.p,2,mean))))^2,2,mean)
  
  return(sum(var.p))
  
}
