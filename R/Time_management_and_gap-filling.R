# This script create continuos time series without any gaps and skipped scans
# It needs to be modifed for other time steps than 30-minutes

Time_management<-function(data,X){

# Load the original timestamp
	timestamp<-as.character(data$V1)

	options(digits.secs=3)
	date=as.POSIXlt(as.character(timestamp),origin = "1970-01-01 00:00:00.1",tz="UTC")

	Date<-format(date,"%Y-%m-%d")
	Hour<-as.numeric(format(date,"%H"))
	Minute<-as.numeric(format(date,"%M"))
	Seconds<-as.numeric(format(date,"%OS2"))
	Seconds<-round(Seconds/dSampling)*dSampling
	Date<-paste(Date,' ',Hour,':',Minute,':',Seconds,sep='')

	date<-as.POSIXlt(strftime(Date,'%Y-%m-%d %H:%M:%OS2'),tz="UTC")

# Create continuous timestamp
	length<-1/dSampling*60*Time_step

	hour=str_sub(paste(0,as.numeric(format(date[1],"%H")),sep=""),-2,-1)

	minute=as.numeric(format(date[1],"%M"))

	if(minute<Time_step){
		M=0
		}else{
		M=Time_step
	}

	Ini=paste(as.POSIXlt(as.character(timestamp[1]),origin='1970-01-01 00:00:00.0',tz='UTC',format='%Y-%m-%d'),
	' ',hour,':',M,':','00.0',sep='')

	Time<-rep(seq((as.POSIXlt(as.character(Ini),tz="UTC")),length.out=Time_step*60,by=1),Fs)

	Time<-Time[order(Time)]

	Time<-Time[2:length(Time)]
	Time<-as.POSIXlt(c(Time,Time[length(Time)]+1),tz="UTC")

	Date<-format(Time,"%Y-%m-%d")
	Hour<-as.numeric(format(Time,"%H"))
	Minute<-as.numeric(format(Time,"%M"))
	Seconds<-as.numeric(format(Time,"%S"))
	seconds<-rep(c(seq(dSampling,60-dSampling,by=dSampling),0),Time_step)
	Date<-paste(Date,' ',Hour,':',Minute,':',seconds,sep='')

	Time<-as.POSIXlt(strftime(Date,'%Y-%m-%d %H:%M:%OS2'),tz="UTC")

# Merge the original data with the new continuos NA series
	Original<-data.frame(as.character(date),X)
	names(Original)[1:2]<-c('Time','X')
	Data<-data.frame(as.character(Time),NA,check.names=FALSE)
	names(Data)[1:2]<-c('Time','X')
	merged<-rbind(Original,Data)

# Remove the duplicates
	filtered<-subset(merged,!duplicated(merged[,1]))

# Order the data frame accroding to continuous timestamp
	Final_order=filtered[order(as.POSIXlt(filtered[,1])),]

# First and last non-NA index
	First=min(which(!is.na(X)))
	Last=max(which(!is.na(X)))

# Use first and last non-NA values from the original time series as the first and last values of the final series
	Y=Final_order$X
	Y[1]=X[First]
	Y[length(Y)]=X[Last]

# Misssing data are linearly interpolated
	Y_interpol<-na.approx(Y)

return(Y_interpol)
}
