# Functions to estimate the derivative(s) of the (mean) spectrum

calc_deriv_simple <- function(x){
	#add statement to check that length of x >= 3
	dx = rep(NA,length(x))
	dx[1] = x[2] - x[1]
	dx[length(x)] = x[length(x)] - x[length(x)-1]
	for(i in 2:(length(x)-1)){
		dx[i] = (x[i+1] - x[i-1])/2
	}
	return(dx)
}

calc_doppler_component_simple <- function(lambda,flux){
	# add statement to check that length of lambda == length of flux
	dlambdadpix = calc_deriv_simple(lambda)
	dfluxdpix = calc_deriv_simple(flux)
	doppler_basis = dfluxdpix*(lambda/dlambdadpix)
	return(doppler_basis)
}

subtract_first_pc_simple <- function(X, fixed_comp){
		
	Xtmp = matrix(NA,dim(X)[1],dim(X)[2])
	fixed_comp_norm = 1/sum(fixed_comp^2)
	Xtmp = X - X%*%fixed_comp%*%t(fixed_comp)*fixed_comp_norm
	
	return(Xtmp)
}
	
make_noisy_spectrum <- function(spec, wavelengths, snr=200, sampling=3){
	
	snr_per_pixel = snr/sqrt(sampling)
	integrated_flux = trapz(wavelengths,spec)
	scale_fac = integrated_flux/sum(diff(wavelengths))
	scaled_spec = (spec*snr*snr)/(scale_fac)
	noisy_spec = scaled_spec + rnorm(length(spec))*sqrt(scaled_spec)
	want_out = (scale_fac*noisy_spec)/(snr*snr)
	
	return(want_out)
	
}

getLower <- function(samps) quantile(samps,0.025)
getUpper <- function(samps) quantile(samps,0.975)