##############################################################################################
#                                                                                            #
#      ---------------------------------------------                                         #
#       Plotting spectral densities of temeperature                                          #
#      ---------------------------------------------                                         #
#                                                                                            #
#       Author: Milan Fischer                                                                #
#       email: fischer.milan@gmail.com                                                       #
#                                                                                            #
##############################################################################################

	rm(list=ls())
	library(plotrix)
	library(stringr)
	library(signal)
	library(zoo)
	library(zoom)

	#Functions to be used
	source('Despike.R')
	source('Ensemble_averages.R')
	source('FFT.R')
	source('Linear_detrend.R')
	source('Spectral_density.R')
	source('Taper.R')
	source('Time_management_and_gap-filling.R')

	#Sampling frequency
	Fs=10

	#Number of bins for ensemble spectra averaging
	N=20

	#Time step (minutes)
	Time_step=30

	#the sampling period (s)
	dSampling=1/Fs

	#Move one directory up
	setwd('..'); MainWD_link=getwd()
	
	#Move to directory with the data
	setwd(paste(getwd(),'/Data',sep=''))

	#Create a list of raw data files
	file<-list.files(path=getwd(),pattern="\\.dat$")  

PlotSpectra<-function(file){

	data<-read.table(file,sep=",",skip=4,na.strings="NAN")

	#Timestamp (determined as the the end of the interval)
	timestamp=str_sub(as.character(data$V1[length(data$V1)]),1,16)

	#Air temperature (°C)

	ts_FREE=data$V2
	ts_IN=data$V3
	ts_WALL=data$V4

	# Despike the temperature time series
	Tde_FREE=Despike(data,ts_FREE)
	Tde_IN=Despike(data,ts_IN)
	Tde_WALL=Despike(data,ts_WALL)

#################################################################### 
### Trend removal                                                  #
#################################################################### 

	T_dlm_FREE=Detrend(Tde_FREE)$Detrended
	T_dlm_IN=Detrend(Tde_IN)$Detrended
	T_dlm_WALL=Detrend(Tde_WALL)$Detrended

#################################################################### 
### Windowing                                                      #
#################################################################### 

	T_dlmW_FREE=Taper(T_dlm_FREE,100/12/2)$Tapered
	T_dlmW_IN=Taper(T_dlm_IN,100/12/2)$Tapered
	T_dlmW_WALL=Taper(T_dlm_WALL,100/12/2)$Tapered

####################################################################
### Fourier transformation                                         #
#################################################################### 

	out=FFT(T_dlmW_FREE)
	even_or_odd=out$even_or_odd
	f_Hz=out$f_Hz
	Nyquist=out$Nyquist
	Fts_FREE=fft(T_dlmW_FREE)/length(T_dlmW_FREE)
	rm(out)


	out=FFT(T_dlmW_IN)
	even_or_odd=out$even_or_odd
	f_Hz=out$f_Hz
	Nyquist=out$Nyquist
	Fts_IN=fft(T_dlmW_IN)/length(T_dlmW_IN)
	rm(out)

	out=FFT(T_dlmW_WALL)
	even_or_odd=out$even_or_odd
	f_Hz=out$f_Hz
	Nyquist=out$Nyquist
	Fts_WALL=fft(T_dlmW_WALL)/length(T_dlmW_WALL)
	rm(out)
	
####################################################################
### Spectral density (Stull 1988, p. 312-313)                      #
#################################################################### 

	SD_FREE=Spectral_density(even_or_odd,f_Hz,Nyquist,Fts_FREE)
	SD_IN=Spectral_density(even_or_odd,f_Hz,Nyquist,Fts_IN)
	SD_WALL=Spectral_density(even_or_odd,f_Hz,Nyquist,Fts_WALL)

# Parseval's theorem (energy conservation-the spectrum intergral is equal to scalar variance)
	
	df_Hz=f_Hz[2]-f_Hz[1]	

	sum((T_dlmW_FREE-mean(T_dlmW_FREE))^2)/length(T_dlmW_FREE)	# Biased variance
	sum((abs(Fts_FREE))^2)
	sum(SD_FREE*df_Hz)

	sum((T_dlmW_IN-mean(T_dlmW_IN))^2)/length(T_dlmW_IN)	# Biased variance
	sum((abs(Fts_IN))^2)
	sum(SD_IN*df_Hz)

	sum((T_dlmW_WALL-mean(T_dlmW_WALL))^2)/length(T_dlmW_WALL)	# Biased variance
	sum((abs(Fts_WALL))^2)
	sum(SD_WALL*df_Hz)

####################################################################
### Ensemble averages                                              #
#################################################################### 

	out=Ensemble(N,f_Hz,Nyquist,SD_FREE)
	F_FREE=out$F
	S_FREE=out$S
	rm(out)

	out=Ensemble(N,f_Hz,Nyquist,SD_IN)
	F_IN=out$F
	S_IN=out$S
	rm(out)

	out=Ensemble(N,f_Hz,Nyquist,SD_WALL)
	F_WALL=out$F
	S_WALL=out$S
	rm(out)

################################
### Plotting                 ###
################################

	Intertial_range=(10^F_FREE)^(-5/3)

	#Log-log presentation of power spectra

	par(mfrow=c(1,1),oma=c(3.5,3.5,0,0)+1,mar=c(0,0,0,0),xpd=NA)
	plot(log10(f_Hz[2:(which(f_Hz==Nyquist))]),log10(SD_FREE[2:(which(f_Hz==Nyquist))]),type='l',xaxt = 'n',xlab='Frequency (Hz)',
	yaxt = 'n',ylab='Spectral density (K*K)',xlim=c(log10(f_Hz[2]),log10(Nyquist)),ylim=c(-12,2),cex=1.2)
	axis(1, at=c(log10(0.001),log10(0.01),log10(0.1),log10(1),log10(5)), labels=c(0.001,0.01,0.1,1,5))
	axis(2, at=c(-12,-10,-8,-6,-4,-2,0,2), labels=c(10^-12,10^-10,10^-8,10^-6,10^-4,10^-2,10^0,10^2))
	clip(x1=log10(f_Hz[2])-(log10(Nyquist)+log10(f_Hz[2]))*0.04, x2=log10(Nyquist)+(log10(Nyquist)-log10(f_Hz[2]))*0.04, y1=-14-14*0.04,y2=2+14*0.04)
	lines(F_FREE,S_FREE,col='blue',lwd=3)
	lines(F_IN,S_IN,col='grey',lwd=2)
	lines(F_WALL,S_WALL,col='red',lwd=1)

	text(x=-1.051317,y=1.75,labels='-5/3',col='red')
	text(x=0.1,y=1.9,labels=timestamp,cex = 1.2)

	Fig_offset=-0.5
	fit=lm((log10(Intertial_range[2:20])-Fig_offset)~F_FREE[2:20])
	abline(fit,col='red',lwd=1,lty=2)

	box()
	T_FREE<-bquote(T['mean'] == .(sprintf("%.3f",mean(ts_FREE)))*';')
	T_IN<-bquote(T['mean'] == .(sprintf("%.3f",mean(ts_IN)))*';')
	T_WALL<-bquote(T['mean'] == .(sprintf("%.3f",mean(ts_WALL)))*';')
	text(x=-2.55,y=-9.8,labels=T_FREE,cex = 1.2,col='blue')
	text(x=-2.55,y=-10.8,labels=T_IN,cex = 1.2,col='grey')
	text(x=-2.55,y=-11.8,labels=T_WALL,cex = 1.2,col='red')

	Var_FREE<-bquote(T['var'] == .(sprintf("%.3f",var(ts_FREE))))
	Var_IN<-bquote(T['var'] == .(sprintf("%.3f",var(ts_IN))))
	Var_WALL<-bquote(T['var'] == .(sprintf("%.3f",var(ts_WALL))))
	text(x=-1.45,y=-9.8,labels=Var_FREE,cex = 1.2,col='blue')
	text(x=-1.45,y=-10.8,labels=Var_IN,cex = 1.2,col='grey')
	text(x=-1.45,y=-11.8,labels=Var_WALL,cex = 1.2,col='red')

}

#################################
### Aplying the main function ###
#################################

pdf(paste(MainWD_link,'/Output/Power_spectrum.pdf',sep=''))

for(i in 1:length(file)){
try(PlotSpectra(file[i]),silent=TRUE)
}
dev.off()
