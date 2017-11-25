#################################################################### 
###   Taper window   (tapered by set percent at the both ends    ###
#################################################################### 

# When 100 is divided by the set percentage*2, it has to result in even integer

Taper<-function(X,percent)
{

	# determine the time series length
	n=length(X)

	# create an interger sequence from 1 to the lenght of the series
	t=seq_len(n)

	# Total percentage (both ends)
	Total_percent=percent*2

	Taper=rep(1,n)

	for(i in 1:round(percent/100*n)){
	Taper[i]=sin(100/Total_percent*pi*t[i]/n)^2
	}

	for(i in (n-round(percent/100*n)+1):n){
	Taper[i]=sin(100/Total_percent*pi*t[i]/n)^2
	}

	Tapered=as.numeric(X*Taper)

	out=list(Taper=Taper,Tapered=Tapered)
}