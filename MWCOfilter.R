
###	function to take peptide strings, compute mass, and filter by cutoff

MWCOfilter=function(temppeps,cutoff=10000){
	masses<-sapply(FUN=calculatePepM,temppeps)
	names(masses)<-temppeps
	flowthrough<-names(masses[masses<cutoff])
	overcutoff<-names(masses[masses>=cutoff])
	filtered<-list(flowthrough,overcutoff)
	return(filtered)
	}