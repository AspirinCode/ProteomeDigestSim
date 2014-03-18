
read.proteome=function(input=fastapicker(),dir=getwd(),cleave=TRUE,pie=FALSE){
    setwd(dir)
    ##### read .fasta file and fill up object
    #####    slots with useful tidbits
    require(seqinr)
    require(gplots)
    object<-new("proteome")
    object@fasta<-read.fasta(input,seqtype="AA",as.string=TRUE)
    object@seq<-getSequence(object@fasta)
    object@names<-names(object@fasta)
    x<-unlist(object@seq)
    aatotal<-length(x)
    ###empty summary matrix
    summary<-matrix(c(rep(0,40)),nrow=20,dimnames=list(
        c("A","S","T","G","V","C","N","L","I","M","P","Y","W","Q","F","D","E","H","K","R"),
        c("count","ratio")),
        ncol=2)
	summary[,1]<-c(length(x[x=="A"]),length(x[x=="S"]),length(x[x=="T"]),
      length(x[x=="G"]),length(x[x=="V"]),length(x[x=="C"]),
	length(x[x=="N"]),length(x[x=="L"]),length(x[x=="I"]),
	length(x[x=="M"]),length(x[x=="P"]),length(x[x=="Y"]),
	length(x[x=="W"]),length(x[x=="Q"]),length(x[x=="F"]),
	length(x[x=="D"]),length(x[x=="E"]),length(x[x=="H"]),
	length(x[x=="K"]),length(x[x=="R"]))
	for(i in 1:length(summary[,1])){summary[i,2]<-summary[i,1]/aatotal}
	par(mfcol=c(1,2))
	if(pie==TRUE){
		pie(summary[,2],main=paste(input,"stats"),radius=1)
		}
	textplot(summary)
	object@AAratios<-summary
	print(aatotal)
	if(cleave==TRUE){
		object<-cut.predict(object=object,residues=c("R","K"),protease="Trypsin")
		}
	return(object)
	} 
