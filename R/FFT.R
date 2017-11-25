FFT<-function(T_dlmW)
{

	# Nyquist frequency (note the difference between even and odd length of time series)
	even_or_odd=ifelse(length(T_dlmW) %% 2==0,'even','odd')

	#Frequency (Hz)
	if(even_or_odd=='even'){
	f_Hz=Fs*seq(from=0,to=1,length.out=length(T_dlmW)+1)
	f_Hz[(length(T_dlmW)/2+2):length(f_Hz)]=f_Hz[(length(T_dlmW)/2+2):length(f_Hz)]-Fs
	f_Hz=f_Hz[1:length(T_dlmW)]
	}else{
	f_Hz=Fs*seq(from=0,to=1,length.out=length(T_dlmW)+2)
	f_Hz[(length(T_dlmW)/2+1.5):length(f_Hz)]=f_Hz[(length(T_dlmW)/2+1.5):length(f_Hz)]-Fs
	f_Hz=f_Hz[-(length(T_dlmW)/2+1.5)][1:length(T_dlmW1)] # Exclude the pseudo-nyquist frequecny and shorten the vector to the scalar length
	}

	Nyquist=max(f_Hz)

# Fast Fourier transformation:
	Fts<-fft(T_dlmW)/length(T_dlmW)

out<-list(even_or_odd=even_or_odd,f_Hz=f_Hz,Nyquist=Nyquist,Fts=Fts)
}