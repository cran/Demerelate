# Changed 18.1.2016

Demerelate <- function(inputdata, tab.dist = "NA", ref.pop = "NA", object = FALSE, value = "Mxy", Fis = FALSE, file.output = FALSE, p.correct = FALSE, iteration = 1000, pairs = 1000, dis.data = "relative", NA.rm = TRUE)

{

    message("\n","---- Demerelate v0.9 ----","\n","\n")
    
    # Data are loaded for input
    if (object==FALSE) {
                          tab.pop <- input.txt(inputdata, "pop")
                          if (tab.dist!="NA") {tab.dist <- input.txt(tab.dist, "dist")}
                          if (ref.pop=="NA") {ref.pop <- tab.pop} else {ref.pop <- input.txt(ref.pop, "ref.pop")}
                        }
    
    if (object==TRUE) {
                          tab.pop <- inputdata
                          if (is(ref.pop)[1]!="data.frame") {ref.pop <- inputdata}
                       }
  
    
    if (NA.rm==TRUE)
    {
    if (is(tab.dist)[1]=="data.frame") {tab.dist <- tab.dist[complete.cases(tab.pop),]}
    ref.pop <- ref.pop[complete.cases(ref.pop),]
    tab.pop <- tab.pop[complete.cases(tab.pop),]
    }
    
    ## Error.checking
    if (is(tab.dist)[1]=="data.frame")
    {
      if (length(tab.pop[,1])!=length(tab.dist[,1]))
        {warning("The genetic inputdata and the geographic positioning data do not comply in length. Obviously, some data are missing here. You may want to check that. Alternatively you may run analyses without tab.dist data.")
         stop()}
      if (length(which((tab.pop[,1]%in%tab.dist[,1])==FALSE)>0)) 
        {warning("Not for every sample in your inpudata geographic data are available. You have to check completeness of tab.dist before continueing with analyses.")
          stop()} 
      if (length(which((tab.dist[,1]%in%tab.pop[,1])==FALSE)>0)) 
        {warning("Not for every sample in your geographic data inpudata are available. You have to check completeness of inpudata before continueing with analyses.")
         stop()} 
    
    }
    
    number.loci <- (length(tab.pop)-2)/2

    # Setting barcode and path for global output    
    barcode <- round(rnorm(1),4)
    
    
    if (file.output==TRUE) {
    out.name <- as.character(barcode)
    inputdata <- deparse(substitute(inputdata))

    if (object==FALSE) {directory.name <- paste("Demerelate",out.name,sep="")}
    if (object==TRUE) {directory.name <- paste("Demerelate",out.name,sep="")}
    
    dir.create(directory.name)}
    
    else {
      
      directory.name <- "NA"
      out.name <- "NA"
    }
    

    # Start Fis calculations on the whole dataset for each population defined in tab.pop[,1]

    if (Fis==TRUE)
    {
      fstatistic.return <- F.stat(tab.pop, TRUE, iteration, directory.name, out.name)
    }
    

    # Relatedness calculations for different populations in inputdata          

            tab.pop.pop <- split(tab.pop,tab.pop[,2])
            if (is(tab.dist)[1]=="data.frame") {tab.dist.dist <- split(tab.dist,tab.dist[,2])}
            empirical <- vector("list",length(length(tab.pop.pop)))
            chisquare <- vector("list",length(length(tab.pop.pop)))
            lin.out <- empirical

	    # Calculating threshold statistics by general linearized model and modelling sibling populations from reference populations
            message(paste("---","Relatedness calculations are performed for reference populations based on overall allelefrequencies","----", Sys.time()),"\n")
            relate.return <- relate.calc(tab.pop, pairs, file.output, value, directory.name, ref.pop)
            Thresholds <- relate.return[[5]]
    
            # Empirical allele sharing and statistical output for each population and distance calculations
            
            for (k in 1:length(tab.pop.pop))
              {
          
                      
               	stat.out <- stat.pops(Thresholds, tab.pop.pop[[k]], pairs, p.correct, directory.name, out.name, file.output, inputdata, object, value, iteration, ref.pop)
               	empirical[[k]] <- stat.out[[1]]        
                chisquare[[k]] <- stat.out[[2]] 
               
                        # Calculations on distances
                        dist.m <- "NA"
                        if (is(tab.dist)[1]=="data.frame") {dist.m <- geo.dist(tab.dist.dist[[k]], tab.dist.dist[[k]], FALSE, dis.data)}
                        names(empirical)[k] <- as.character(tab.pop.pop[[k]][1,2])
                        names(chisquare)[k] <- as.character(tab.pop.pop[[k]][1,2])
               	  
            
                        # Statistics and export of plots
                        
                        if (is(tab.dist)[1]=="data.frame") {
                          empirical.dist <- data.frame(do.call("rbind",strsplit(names(empirical[[k]]),"_")),empirical[[k]])
                          empirical.dist <- with(empirical.dist,
                                                 structure(
                                                   empirical.dist[,3],
                                                   Size = length(unique(c(as.character(empirical.dist[,1]),as.character(empirical.dist[,2])))),
                                                   Labels = unique(c(as.character(empirical.dist[,1]),as.character(empirical.dist[,2]))),
                                                   Diag = FALSE,
                                                   Upper = FALSE,
                                                   method = "user",
                                                   class = "dist"))
                                lin.out[[k]] <- Lin.reg.distance(dist.m, empirical.dist, pairs, tab.pop.pop[[k]], relate.return[[3]], relate.return[[2]], relate.return[[4]], file.output, directory.name, out.name, inputdata, object, value, iteration)
                                names(lin.out)[k]<-as.character(tab.pop.pop[[k]][1,2])
                                                              }
            
                        # Export geodistance matrices for each population

if (file.output==TRUE)
{
                if (is(tab.dist)[1]=="data.frame") {write.table(file=paste(".","/",directory.name,"/","Geo.distances.",tab.pop.pop[[k]][1,2],".txt",sep=""),x=round(dist.m,3), quote=FALSE, sep=" ")}
                
                
                write.table(file=paste(".","/",directory.name,"/","Empirical.relatedness.",tab.pop.pop[[k]][1,2],".txt",sep=""),x=empirical[[k]], quote=FALSE, sep=" ", col.names=value)
}

              }
    
    empirical_rel <- empirical
      

 	# If value continous boxplot is given as output to compare mean relatedness between populations

            	empirical[[length(empirical)+1]]<- relate.return[[4]]
                names(empirical)[length(empirical)] <- "Non-related"
            	empirical[[length(empirical)+1]]<- relate.return[[3]]
                names(empirical)[length(empirical)] <- "Half Siblings"
            	empirical[[length(empirical)+1]]<- relate.return[[2]]
                names(empirical)[length(empirical)] <- "Full Siblings"
                
                
                # Boxplot grafical output
            	if (file.output==TRUE)
            	{
                pdf(paste(".","/",directory.name,"/","Relatedness of populations.pdf",sep=""))
            
                par(las=1, cex.axis=0.5)
                boxplot(empirical, main="Mean relatedness of populations", names=c(names(tab.pop.pop),"Non-related","Half siblings","Full siblings"), ylab=paste(value,"estimate of relatedness"))
                
                dev.off()
            	}
                

                # Preparations on T test statistics for comparison between populations
                e.values <- unlist(lapply(empirical,as.vector))
                groups <- unlist(lapply((lapply(empirical,as.vector)),length))
                gr.vect <- as.character(vector()) 
            
            	for (i in 1:length(groups))
                  {
                    gr.vect <- c(gr.vect,(rep(names(empirical)[i],groups[i])))
                  }
      
if (file.output==TRUE)
{
      # Pairwise T test Output of statistics
      out.file <- file(paste(".","/",directory.name,"/","Pairwise.t.txt",sep=""),"w")
      
      writeLines(paste(
        "Demerelate - v.0.9", " ---","\n",
        "Relatedness outputfile on file: ", inputdata,"\n",
        "Analysis had been made based on ",pairs," pairs "," using the ",value," estimator.","\n",
        if (value=="Bxy"){paste("Calculations are based on Li and Horvitz 1953. The values represent an indication on relatedness based on allele sharing.","\n", sep="")},
        if (value=="Mxy"){paste("Calculations are based on Bluoin et al. 1996. The values represent relatedness assessment based on genotype sharing.","\n", sep="")},
        if (value=="rxy"){paste("Calculations are based on Queller and Goodnight 1989. The values represent relatedness value corrected for total allele diversity.","\n", sep="")},
        if (value=="Sxy"){paste("Calculations are based on Lynch 1988. The values represent an indication on relatedness based on proportion of shared alleles.","\n", sep="")},
        if (value=="lxy"){paste("Calculations are based on Lynch and Ritland 1999. Final values represent relatedness estimates and are averaged via RE-RAT procedure.","\n", sep="")},
        if (value=="ritland"){paste("Calculations are based on Ritland 1996 with slight modifications based on Lynch and Ritland 1999.","\n", sep="")},
        if (value=="wang.fin"){paste("Calculations are based on Wang 2002. The estimates should only be used for finite populations.","\n", sep="")},
        if (value=="wang"){paste("Calculations are based on Wang 2002. The estimates is corrected for sample size and for expected frequency of unrelated individuals according to Li et al 1993.","\n", sep="")},
        if (value=="morans"){paste("Calculations are based on Hardy and Vekemans 1999 which describes morans I as estimator for genetic relatedness.","\n", sep="")},
        if (value=="morans.fin"){paste("Calculations are based on Hardy and Vekemans 1999 but ommiting sample size corrections. This value is only applicable for finite populations.","\n", sep="")},
        if (value=="loiselle"){paste("Calculations are based on Loiselle et al. 1995 ans are implemented as described by Hardy and
Vekemans 2015.","\n", sep="")},
        if (value=="Li"){paste("The estimator is calculated according to Li et al. 1993","\n", sep="")},
        "For further information mind References at the end of this file.","\n","\n",
        "\n",
        "Summary on Pairwise T-Test","\n",
        "---",sep=""),con=out.file)
      
      writeLines(pairwise.t.test(e.values,gr.vect)[[1]],con=out.file)
      write.table(round(as.data.frame(pairwise.t.test(e.values,gr.vect)[[3]]),3),out.file, append=T, quote=F, sep="\t", col.names = NA)
      
      writeLines(paste(
        "---","\n","\n",
        "\n","\n","\n","References","\n",
        "Blouin, M., Parsons, M., Lacaille, V. and Lotz, S. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.","\n",
        "Hardy, O.J. and Vekemans, X. (1999) Isolation by distance in a contiuous population: reconciliation between spatial autocorrelation analysis and population genetics models. Heredity, 83, 145-154.","\n",
        "Li, C.C., Weeks, D.E. and Chakravarti, A. (1993) Similarity of DNA fingerprints due to chance and relatedness. Human Heredity, 43, 45-52.","\n",
        "Li, C.C. and Horvitz, D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics, 5, 107-17.","\n",
        "Loiselle, B.A., Sork, V.L., Nason, J. and Graham, C. (1995) Spatial genetic structure of a tropical understory shrub, Psychotria officinalis (Rubiaceae). American Journal of Botany, 82, 1420-1425.","\n",
        "Lynch, M. (1988) Estimation of relatedness by DNA fingerprinting. Molecular Biology and Evolution, 5(5), 584-599.","\n",
        "Lynch, M. and Ritland, K. (1999) Estimation of pairwise relatedness with molecular markers. Genetics, 152, 1753-1766.","\n",
        "Mantel, N. (1967) The detection of disease clustering and a generalized regression approach. Cancer Research, 27, 209-220.","\n",
        "Oksanen, J. et al. (2013) vegan: Community Ecology Package. R package version 2.0-8.","\n",
        "Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.","\n",
        "Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.","\n",
        "Ritland, K. (1999) Estimators for pairwise relatedness and individual inbreeding coefficients. Genetics Research, 67, 175-185.","\n",
        "Wang, J. (2002) An estimator for pairwise relatedness using molecular markers. Genetics, 160, 1203-1215.","\n",
        sep=""), con=out.file)
      
      close(out.file)
    
}    

 # Information on overall result location

 settings <- cbind(c("Barcode", "Date", "Object", "Relatedness-Value", "Randomized Pairs", "Fis", "Fis-Iterations"),c(barcode, date(), object, value, pairs, Fis, iteration))
 dem.results <- list(settings, empirical_rel, list(relate.return[[2]],relate.return[[3]],relate.return[[4]]), as.data.frame(pairwise.t.test(e.values,gr.vect)[[3]]), chisquare, Thresholds, if(Fis==TRUE){fstatistic.return}, lin.out)    
 names(dem.results) <- c("Settings", "Empirical_Relatedness", "Randomized_Populations_for_Relatedness_Statistics", "Relatedness_Statistics (T-test)", "Relatedness_Statistics (X^2-Test)", "Thresholds for relatedness", "Fis_Statistics", "Mantel_Correlation_of_Relatedness")     

if (file.output==FALSE) {return(dem.results)} 
if (file.output==TRUE) {return(list(dem.results,message(paste("\n ---","Please find your results in folder:", directory.name," ---- \n"))))}

    
}
