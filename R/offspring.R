offspring <- function(parent1,parent2,allele.column, pairs)
		{## remove pairs
		
			# Offspring is calculated from parent1 crossed with parent2 (class = input/data.frame)
			# Output == population according to format from Input.txt; offspr (class = data.frame)
			# allele.column... allele.column+1 as source for offspring
		
		locusA <- apply(parent1[,allele.column:(allele.column+1)],1,function(x){sample(x,prob=c(.5,.5),size=1)})
		locusB <- apply(parent2[,allele.column:(allele.column+1)],1,function(x){sample(x,prob=c(.5,.5),size=1)})
		
		offspring <- t(apply(data.frame(locusA,locusB),1,function(x){sort(x,decreasing=TRUE)}))
		
		offspring <- data.frame(offspring,offspring)

		offspring <- remove.na.rows(offspring)
		return(offspring)
			
		}
		
