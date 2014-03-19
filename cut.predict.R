################################################################################
## cut predictor v2; for use with multiple cleavage points
###########################
####				9/1/2011 THIS WORKS with an arbitrary number of residues
##### version three, allow cleavage at n-terminal or c-terminal
####	use a vector of length == length(residues) where
####	0 equals n-terminal to the residue, and 1 equals c-terminal to the residue
#####	
####### 032712 this works to cleave at either n-term or c-term
##########		#############		###########
###	FIXED: the loop when term=0 needs to start on 2 instead of 1
###  	added a propensity for missed cleavage

cut.predict.test=function(object=spombe,residues=c("D"),term=c(0),protease="AspN",add=FALSE, missed.propensity=0.01){
	###	testing the cleavage n-terminal to residues[1]
	#########
	proteins.l<-length(object@fasta)
	peps<-c(rep(0,sum(nchar(unlist(object@seq)))))
	peps2<-c(rep(0,sum(nchar(unlist(object@seq)))))
	pep.lengths<-c(rep(0,sum(nchar(unlist(object@seq)))))
	peps.l<-length(peps)
	lengths<-c()
	print(residues[1])
	m<-1###the location in peps
	for(i in 1:proteins.l){
		pro<-unlist(object@seq[i])
		pro.l<-length(pro)
		###loop over AAs in pro[i]
		k<-1
		if(term[1]==0){
			for(j in 2:pro.l){
				if(pro[j]==residues[1] & runif(1)>=missed.propensity){
					end<-j-1
					peps[m]<-c(paste(pro[k:end],collapse=""))
					k<-j
					m<-m+1
					}
				if(j==pro.l){
					peps[m]<-c(paste(pro[k:j],collapse=""))
					m<-m+1
					}
				}
			}
		##### if cleavage is supposed to be at the n-terminal 
		if(term[1]==1){
			for(j in 1:pro.l){
				if(pro[j]==residues[1] & runif(1)>=missed.propensity){
					peps[m]<-c(paste(pro[k:j],collapse=""))
					k<-j+1
					m<-m+1
					}
				else if(j==pro.l){
					peps[m]<-c(paste(pro[k:j],collapse=""))
					m<-m+1
					}
				}
			}
		}
	if(length(residues)>1){
		for(n in 2:length(residues)){
			print(residues[1:n])
			peps.l<-length(peps[peps!=0])
			m<-1
			for(i in 1:peps.l){
        			pep<-unlist(strsplit(peps[i],split=""))
				pep.l<-length(pep)
				###loop over AAs in pro[i]
				k<-1
				for(j in 1:pep.l){
					if(pep[j]==residues[n] & runif(1)>=missed.propensity){
					##if the jth residue=residues
					###write pep 1:j
						peps2[m]<-c(paste(pep[k:j],collapse=""))
						k<-j+1
						m<-m+1
						}
					###if the jth residue == the last residue
					else if(j==pep.l){
						peps2[m]<-c(paste(pep[k:j],collapse=""))
						k<-j+1
						m<-m+1
						}
					}
				}
			peps<-peps2
			}
		}
	peps<-peps[peps!=0]
	### this works but there is a bug that allows the same set twice
	if(add){
		object@cleaved[[protease]]<-peps
		}
	if(add==FALSE){
		object@cleaved<-list(protease=peps)
		}
	return(object)
	}
