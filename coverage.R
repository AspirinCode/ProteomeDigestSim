######################################################################################
#################################################################
########		script to compute proteome coverage given cleavage 
######## 		output from cut.predict

coverage=function(object=Spombe){
	#### how many proteases were used?
	cleaved.l<-length(object@cleaved)
	###loop over the cleaved.l
	for(j in 1:cleaved.l){
		peps2<-unlist(object@cleaved[j])
		totalvec<-c(rep(0,sum(nchar(peps2))))
		peplengths<-nchar(peps2)
		pepvec.l<-length(peplengths)
		start<-1
		end<-peplengths[1]
		if(peplengths[1]>=7&peplengths[1]<=35){
			totalvec[start:end]<-rep(1,times=peplengths[1])
			}
		start<-end+1
		for(i in 2:pepvec.l){
			end<-peplengths[i]+end
			if(peplengths[i]>=7&&peplengths[i]<=35){
				totalvec[start:end]<-rep(1,times=peplengths[i])
				}
			start<-start+peplengths[i]
			print(i)
			}
		### this part puts the binary outcome into object@coverage
		if(j==1){object@coveragevec<-totalvec}
		if(j>=2){object@coveragevec<-object@coveragevec+totalvec}
		}
	notcovered<-length(totalvec[totalvec==0])
	total<-length(totalvec)
	pctnotcovered<-notcovered/total
	pctcovered<-1-pctnotcovered
	object@coverage<-c(covered=pctcovered,not=pctnotcovered)
	pie(object@coverage,main="predicted coverage")
	return(object)
	}
