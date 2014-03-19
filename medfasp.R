#### use this function to simulate MED-FASP cleavage with iterative digestion
#### As used in the manuscript - 10,000 MWCO filtration during digestion with final filtration based
#### on peptide length
####  set integer seed for random number generation for use with the missed cleavage propensity
#set.seed(420, kind = NULL, normal.kind = NULL)
#runif(1)

medfasp=function(object=cere,proteases=c("lysC"),filterMW=TRUE,filterLEN=FALSE,LENcutoff=c(8,35),MWcutoff=5000, missed.propensity=c()){
	###	iteratively digest with n proteases, removing those below cutoff
	lowCO<-LENcutoff[1]
	highCO<-LENcutoff[2]

	### make sure that if missed.propensity!=0 the length matches # of proteases
	if(length(missed.propensity!=0)){
		if(length(proteases)==length(missed.propensity))print("using missed cleavages")
		if(length(proteases)!=length(missed.propensity)){
			print("# of proteases does not match length of missed.propensity")
			return()
			}
		}
	### first part call cut.predict ###
	if(proteases[1]=="trpC"){
		tempproteome<-cut.predict(object,residues=c("W"),term=c(1),protease="trpC")
		}
	if(proteases[1]=="CNBr"){
		tempproteome<-cut.predict(object,residues=c("M"),term=c(1),protease="CNBr")
		}
	if(proteases[1]=="trypsin"){
		tempproteome<-cut.predict(object,residues<-c("R","K"),term<-1)
		}
	if(proteases[1]=="aspN"){
		tempproteome<-cut.predict(object,residues<-c("D"),term<-0)
		}
	if(proteases[1]=="gluC"){
		tempproteome<-cut.predict(object,residues<-c("E"),term<-1)
		}		
	if(proteases[1]=="argC"){
		tempproteome<-cut.predict(object,residues<-c("R"),term<-1)
		}
	if(proteases[1]=="lysC"){
		tempproteome<-cut.predict(object,residues<-c("K"),term<-1)
		}
	if(proteases[1]=="NTCB"){
		tempproteome<-cut.predict(object,residues<-c("C"),term<-0)
		}
	### decide which peptides are smaller than cutoff
	temppeps<-unlist(tempproteome@cleaved)
	peptides.len<-length(temppeps)
	peplengths<-nchar(temppeps)
	names(temppeps)<-temppeps
	if(filterMW==TRUE){
		filtered<-MWCOfilter(temppeps,MWcutoff)
		flowthrough<-filtered[[1]]
		overcutoff<-filtered[[2]]
		}
	if(filterLEN==TRUE){
		flowthrough<-names(temppeps[peplengths<cutoff[2]])
		overcutoff<-names(temppeps[peplengths>cutoff[2]])
		}
	proteases.l<-length(proteases)
	if(proteases.l>=2){
		### loop over proteases from 2:length(proteases)
		for(q in 2:proteases.l){
			#print(q)
			proteins.l<-length(overcutoff)
			#print(proteins.l)
			### set residues and termini based on protease
			if(proteases[q]=="CNBr"){
				residues<-c("M")
				term<-1
			}
			if(proteases[q]=="trypsin"){
				residues<-c("R","K")
				term<-1
			}
			if(proteases[q]=="aspN"){
				residues<-c("D")
				term<-0
			}
			if(proteases[q]=="gluC"){
				residues<-c("E")
				term<-1
			}		
			if(proteases[q]=="argC"){
				residues<-c("R")
				term<-1
			}
			if(proteases[q]=="lysC"){
				residues<-c("K")
				term<-1
			}
			if(proteases[q]=="trpC"){
				residues<-c("W")
				term<-1
			}
			if(proteases[q]=="NTCB"){
				residues<-c("C")
				term<-0
			}
			## now take the peptides with mass over cutoff and cleave then similar to cut.predict()
			peps<-c(rep(0,sum(nchar(unlist(overcutoff)))))
			peps2<-c(rep(0,sum(nchar(unlist(overcutoff)))))
			m<-1
			residueslen<-length(residues)
				for(i in 1:proteins.l){
					pro<-unlist(strsplit(overcutoff[i],""))
					pro.l<-length(pro)
						###loop over AAs in pro[i]
						k<-1
						if(term[1]==0){
							for(j in 2:pro.l){
								if(pro[j]==residues[1]  & runif(1)>=missed.propensity){
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
						if(residueslen==1){
							for(j in 1:pro.l){
								if(pro[j]==residues  & runif(1)>=missed.propensity){
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
						if(residueslen==2){
							for(j in 1:pro.l){
								if(pro[j]==residues[1]  & runif(1)>=missed.propensity |pro[j]==residues[2] & runif(1)>=missed.propensity){
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
					# for debugging
					#print(i/proteins.l)
					}
			peps<-peps[peps!=0]
			peptides.len<-length(peps)
			peplengths<-nchar(peps)
			names(peps)<-peps
			#### 	filter new peptide populations
			if(filterLEN==TRUE){
				flowthrough<-names(peps[peplengths<cutoff[2]])
				overcutoff<-names(peps[peplengths>cutoff[2]])
				}
			if(filterMW==TRUE){
				filtered<-MWCOfilter(temppeps=peps,MWcutoff)
				flowthrough<-c(flowthrough,filtered[[1]])
				overcutoff<-filtered[[2]]
				}
			}
		}
	###### calculate the percent coverage based on peptides above lower LENcutoff
	notcovered<-sum(nchar(overcutoff)) + sum(nchar(flowthrough[nchar(flowthrough)<lowCO]))
	covered<-sum(nchar(flowthrough[nchar(flowthrough)>=lowCO]))
	print(notcovered+covered)
	print("ratio of amino acids within cutoff")
	print(covered/(notcovered+covered))
	return(flowthrough)
	}
