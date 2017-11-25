#################################################################### 
###   Linear trend removal                                       ###
#################################################################### 

Detrend<-function(X)
{
	# determine the time series length
	n1=length(X)

	# create an interger sequence from 1 to the length of the series
	X1=seq_len(n1)

	# calculate the linear trend and offset of the time series
	Intercept=as.numeric(lm(X~X1)$coefficients[1])
	Slope=as.numeric(lm(X~X1)$coefficients[2])

	# remove the trend and offset of time series
	Detrended=X-(X1*Slope+Intercept)

out=list(Intercept=Intercept,Slope=Slope,Detrended=Detrended)
}