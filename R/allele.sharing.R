# Written by P. Kraemer
# Last modified 06.10.2015

allele.sharing <- function(pop1, pop2, allele.column, onlypairs=FALSE, value=NA, ref.pop)
	{ 

	# Function calculates sharing rate of alleles for pop1 and pop2 for each locus in column [allele.column] and [allele.column+1]
    

	if (onlypairs==FALSE)
    {
    # Preparing vectors/dfs for applying
	  data <- expand.grid(seq(1:length(pop1[,1])),seq(1:length(pop2[,1])))
	  data <- data[as.numeric(data[,1]) <= as.numeric(data[,2]),]
	  data <- data[data[,1]!=data[,2],]
	  }
    
	if (onlypairs==TRUE)
    {
	  data <- data.frame(seq(1:length(pop1[,1])),seq(1:length(pop2[,1])))
	  }
	  
    if (value=="wang" | value=="wang.fin")
       {
       simil.vector <- vector(mode="list",length=length(data[,1]))
       simil.vector <- do.call("rbind",lapply(seq(1:length(data[,1])), FUN=value, data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop))
       simil.vector <- data.frame(simil.vector, row.names=apply(expand.grid(pop1[,1],pop2[,1]),1,paste,collapse="_")[as.numeric(row.names(data))])
       }
  
    if (value=="lxy" | value=="lxy.w")
	     {
       simil.vector <- vector(mode="numeric",length=length(data[,1]))
       simil.vector <- sapply(seq(1:length(data[,1])), FUN=value, data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop)
       simil.vector.RE <- simil.vector[1,]
       simil.vector.RAT <- simil.vector[2,]
       simil.vector.RE <- data.frame(simil.vector.RE, row.names=apply(expand.grid(pop1[,1],pop2[,1]),1,paste,collapse="_")[as.numeric(row.names(data))])
       simil.vector.RAT <- data.frame(simil.vector.RAT, row.names=apply(expand.grid(pop1[,1],pop2[,1]),1,paste,collapse="_")[as.numeric(row.names(data))])
       simil.vector <- list(simil.vector.RE,simil.vector.RAT)
       }

	  if (value=="Mxy" | value=="Bxy" | value=="Sxy" | value=="rxy" | value=="Li" | value=="ritland" | value=="ritland.w" | value=="loiselle" | value=="loiselle.w" | value=="morans.fin" | value=="morans.fin.w" | value=="morans" | value=="morans.w")
       {
       simil.vector <- vector(mode="numeric",length=length(data[,1]))
       simil.vector <- sapply(seq(1:length(data[,1])), FUN=value, data=data, pop1=pop1, pop2=pop2, allele.column=allele.column, ref.pop=ref.pop)
       simil.vector <- data.frame(simil.vector, row.names=apply(expand.grid(pop1[,1],pop2[,1]),1,paste,collapse="_")[as.numeric(row.names(data))])        
       }
    
  
	  return(simil.vector)    
	}

