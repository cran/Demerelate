
### P. Kraemer 13.5.2015
## Basic functions from Ritland 1996, with Spagedi Manual (Hardy et al.)
# Last updated 15.12.2015 - Equation 5 from Ritland 1996
# Last updated 10.02.2016 - Equation (9) from Lynch and Ritland 1999 with locus weihgtsa after Sagei

ritland <- function(row, data, pop1, pop2, allele.column, ref.pop=NA)
{
  
  pm <- ref.pop[[allele.column]]
  
  # Relative frequency for specific alleles in ref.pop
  pxma <- pm[which(names(pm)==pop1[data[row,1],(allele.column*2+1)])]
  pxmb <- pm[which(names(pm)==pop1[data[row,1],(allele.column*2+2)])]
  pyma <- pm[which(names(pm)==pop2[data[row,2],(allele.column*2+1)])]
  pymb <- pm[which(names(pm)==pop2[data[row,2],(allele.column*2+2)])]
  
  p<-c(pxma,pxmb,pyma,pymb)[unique(names(c(pxma,pxmb,pyma,pymb)))]
  
  r.return <- sum(sapply(1:length(p),function(x){
    (sum(
      (as.numeric(names(p[x])==names(pxma)) && as.numeric(names(p[x])==names(pyma)))/p[x],
      (as.numeric(names(p[x])==names(pxma)) && as.numeric(names(p[x])==names(pymb)))/p[x],
      (as.numeric(names(p[x])==names(pxmb)) && as.numeric(names(p[x])==names(pyma)))/p[x],
      (as.numeric(names(p[x])==names(pxmb)) && as.numeric(names(p[x])==names(pymb)))/p[x]))})/4)-1
  
  return(2*r.return)

  ## 2* is calculated due to Lynch & Ritland 1999 in order to be compareable to loiselle
}
