Ensemble<-function(N,f_Hz,Nyquist,SD)
{

X<-log10(f_Hz[2:(which(f_Hz==Nyquist))])	
Y<-log10(SD[2:(which(f_Hz==Nyquist))])	# Pre-multiplied form
intervals=seq(from=min(X),to=max(X),length.out=N+1);intervals<-intervals[2:length(intervals)]
F=numeric()
S=numeric()
for(i in 1:length(intervals))
{
if(i %in% 1)
{
F[i]=mean(subset(X,X<=intervals[i]))
S[i]=mean(subset(Y,X<=intervals[i]))
}else{
F[i]=mean(subset(X,X>intervals[i-1]&X<=intervals[i]))
S[i]=mean(subset(Y,X>intervals[i-1]&X<=intervals[i]))
}
}

out<-list(F=F,S=S)
}
