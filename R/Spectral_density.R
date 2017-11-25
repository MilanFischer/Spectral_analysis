Spectral_density<-function(even_or_odd,f_Hz,Nyquist,Fts)
{	
	df_Hz<-f_Hz[2]-f_Hz[1]
	if(even_or_odd=='even'){
	E=c(2*(abs(Fts[1:(which(f_Hz==Nyquist)-1)]))^2,(abs(Fts[which(f_Hz==Nyquist)]))^2)
	}else{
	E=2*(abs(Fts[1:(which(f_Hz==Nyquist))]))^2
	}

	# Spectral density  (see Stull 1988, p. 313)
	SD=E/df_Hz

	return(SD)
}