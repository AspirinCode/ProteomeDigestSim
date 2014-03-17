#### function to take a character string of amino acids and compute peptide mass


calculatePepM=function(charSeq="TEDELQDKIHPF"){
	###table of exact amino acid masses
	AA.masses<-list(A=71.03711,
	R=156.10111,
	N=114.04293,
	D=115.02694,
	C=103.00919,
	Q=128.05858,
	E=129.04259,
	G=57.02146,
	H=137.05891,
	I=113.08406,
	L=113.08406,
	K=128.09496,
	M=131.04049,
	F=147.06841,
	P=97.05276,
	S=87.03203,
	T=101.04768,
	W=186.07931,
	Y=163.06333,
	V=99.06841)
	#######take the input character sequence	
	subsplt<-unlist(strsplit(charSeq,""))
	seq.len<-length(subsplt)
	sum=18.0105######water mass
	for(i in 1:seq.len){
		sum=sum+AA.masses[[subsplt[i]]]

		}
	print(sum)
	}