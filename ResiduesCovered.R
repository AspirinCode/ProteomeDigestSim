#########################################################################
######	ResiduesCovered takes a proteome object with slot @cleaved and 
######	computes the number and ratio of ser, thr, and tyr covered
#########################################################################
######	Function to generate object of class pepsum from MSGFDB
######	output from tabdelim.txt files
######	index=TRUE indexes row numbers for <=0.01 Pepfdr 
#########################################################################
######	usage: read.MSGFDB(dir=getwd(),input="file.txt",index=TRUE)
######	unique=TRUE prints # of unique peptides and proteins, respectively
#########################################################################

ResiduesCovered=function(object=cervis,name="aspN"){
	peps<-object@cleaved
	peplengths<-nchar(peps[[1]])
	sized<-unlist(peps)[nchar(peps[[1]])>=7 & nchar(peps[[1]])<=35]
	unlisted<-unlist(strsplit(sized,""))
	
	residues<-rownames(object@AAratios)
	coveredlist<-list()
	for(x in residues){
		coveredlist[[x]]<-length(unlisted[unlisted==x])
		print(coveredlist)
		}

	pctcovered<-(unlist(coveredlist)/object@AAratios[1:20,1])*100
	##### output the numbers and ratios of each residue covered
	for(i in 1:20){
		print(residues[i])
		print(coveredlist[[i]]/object@AAratios[i,1])
		}
	unlist(coveredlist)
	barplot(pctcovered,ylim=c(0,100),ylab="Percent",xlab="residue",main=paste(name, ", theoretical residue coverage"))
	###line at the total proteome coverage
	abline(h=object@coverage[1]*100)
	}
