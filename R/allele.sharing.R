# Written by P. Kraemer
# Last modified 06.10.2015

allele.sharing <- function(pop1, pop2, allele.column, data, value=NA, ref.pop)
	{ 

	# Function calculates sharing rate of alleles for pop1 and pop2 for each locus in column [allele.column] and [allele.column+1]

  simil.vector <- vector(mode="list",length=length(data[,1]))
  
  
    if (value=="wang" | value=="wang.fin")
       {
       simil.vector <- do.call("rbind",lapply(seq(1:length(data[,1])), FUN="wang", data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop))
       }
  
    if (value=="lxy" | value=="lxy.w")
	     {
       simil.vector <- sapply(seq(1:length(data[,1])), FUN=value, data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop)
       simil.vector <- list(simil.vector[1,],simil.vector[2,])
       }

	  if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li" | value=="ritland" | value=="ritland.w" | value=="loiselle")
       {
       simil.vector <- sapply(seq(1:length(data[,1])), FUN=value, data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop)
	  }
  
    if (value=="morans.fin" | value=="morans")
       {
       pop <- rbind(pop1,pop2)
       ref.pop <- table(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))/length(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))
       simil.vector <- sapply(seq(1:length(data[,1])), FUN="morans.fin", data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop)
    }
  
  if (value=="morans.fin.w" | value=="morans.w")
  {
    pop <- rbind(pop1,pop2)
    ref.pop <- table(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))/length(c(pop[,allele.column*2+1],pop[,allele.column*2+2]))
    simil.vector <- morans.w(pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop)
  }
  
  if (value=="morans")
  {
    if (identical(pop1,pop2)==TRUE){n<-nrow(pop1)}else{n<-nrow(pop)}
    simil.vector <- simil.vector+morans.w(pop1, pop2, allele.column=allele.column, ref.pop=ref.pop)/(n-1)
  }
	  return(simil.vector)
	}

