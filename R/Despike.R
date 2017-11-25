# Despike.r :     removes spikes from a 1-D data vector and fill the gaps (including skipped scans) by the linear interpolation
#            
#        X        = 1-D data vector  
#                    
#        usage: Y = Despike(data,X)

Despike<-function(data,X){

# Loop is ensuring that the test is repeated if needed (Foken 2008, p. 109) - it is important because first spikes can result in large standard deviation and as a consequence smaller spikes could remain
fulfilled='no'
while(fulfilled %in% 'no')
{
	Xm=mean(X,na.rm=TRUE)
	stda=sd(X,na.rm=TRUE)

# We assume 5 standard deviations as a threshold for data removal
	Fc1=-5*stda+Xm
	Fc2= 5*stda+Xm

	# Values which are outside given limits or are missing are assignes as NA
	X[X<Fc1|X>Fc2|X %in% NA]=NA
	
	M=mean(X,na.rm=TRUE)
	SD=sd(X,na.rm=TRUE)
	Low=-5*SD+M
	Up = 5*SD+M
		if(any(X<Low|X>Up,na.rm=TRUE)==TRUE)
		{
			fulfilled='no'
			}else{
			fulfilled='yes'
		}
}

	Y=Time_management(data,X)

	return(Y)}

