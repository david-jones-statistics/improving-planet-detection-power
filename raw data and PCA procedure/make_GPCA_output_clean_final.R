library(caTools)
source("GPCA_functions.R")

# rewrite this to allow user to enter the lat
spot_files <- list.files("reduced_spec_incl90_lat40_size1pc_high_density")
n_spectra <- length(spot_files)

# get the phases from the file names
spot_phase <- NULL
for(i in 1:n_spectra){
	str_tmp <- strsplit(spot_files[i],"phase_")[[1]][2]
	spot_phase[i] <- as.numeric(strsplit(str_tmp,".c")[[1]][1])
}

# look at the first spectra to get the number of wavelengths
tmp <- read.csv(paste("reduced_spec_incl90_lat40_size1pc_high_density/",spot_files[1],sep=""),header=TRUE)
n_wavelengths = length(tmp$Wavelength)
wavelengths <- tmp$Wavelength

# build matrix with rows corresponding to observations and columns to wavelength
spot_data <- matrix(NA,n_spectra,n_wavelengths)
spot_data[1,] <- tmp$Intensity
for(i in 2:n_spectra){
	spot_data[i,] <- read.csv(paste("reduced_spec_incl90_lat40_size1pc_high_density/",spot_files[i],sep=""),header=TRUE)$Intensity
}

# reorder the data by phase
ind <- c(63:1,64:125)
#ind <- c(13:1,14:25)
spot_data <- spot_data[ind,]
spot_phase <- spot_phase[ind]

#
# Not doing any interpolation
#

# Get re-normalized data
spot_data_renorm <- matrix(NA,n_spectra,n_wavelengths)
integrated_flux_quiet = trapz(wavelengths,spot_data[1,])
for(i in 1:n_spectra){
	integrated_flux = trapz(wavelengths,spot_data[i,])
	spot_data_renorm[i,] <- spot_data[i,]*(integrated_flux_quiet/integrated_flux) # this is correct?
}

doppler_comp_simple = calc_doppler_component_simple(wavelengths,apply(spot_data_renorm,2,mean))
doppler_phase_comp_fix = spot_data_renorm%*%doppler_comp_simple
X_fix = subtract_first_pc_simple(spot_data_renorm,doppler_comp_simple)

# doppler_comp_simple = calc_doppler_component_simple(wavelengths,apply(spot_data,2,mean))
# doppler_phase_comp_fix = spot_data%*%doppler_comp_simple
# X_fix = subtract_first_pc_simple(spot_data,doppler_comp_simple)
PCA_fix = prcomp(X_fix,center=TRUE,scale=TRUE)
PC2_fix <- PCA_fix$x[,1]
PC3_fix <- PCA_fix$x[,2]
PC4_fix <- PCA_fix$x[,3]

#rm(spot_data_renorm)

n_boot = 200
num_eigen <- NULL
doppler_phase_comp <- matrix(NA,length(spot_phase), n_boot)
PCA2 <- matrix(NA,length(spot_phase),n_boot)
PCA3 <- matrix(NA,length(spot_phase),n_boot)
PCA4 <- matrix(NA,length(spot_phase),n_boot)

for(j in 1:n_boot){
	
	print(j)
	
	noisy_obs = t(apply(spot_data,1,make_noisy_spectrum,wavelengths=wavelengths,snr=200,sampling=3))
	
	noisy_obs_renorm <- matrix(NA,n_spectra,n_wavelengths)
	integrated_flux_quiet = trapz(wavelengths,noisy_obs[1,])
	for(i in 1:n_spectra){
		integrated_flux = trapz(wavelengths,noisy_obs[i,])
		noisy_obs_renorm[i,] <- noisy_obs[i,]*(integrated_flux_quiet/integrated_flux) # this is correct?
	}
	noisy_obs <- noisy_obs_renorm
	rm(noisy_obs_renorm)
	
	doppler_phase_comp[,j] = noisy_obs%*%doppler_comp_simple
	X = subtract_first_pc_simple(noisy_obs,doppler_comp_simple)
	PCA_proj = predict(PCA_fix,newdata=X)
	PCA2[,j] = PCA_proj[,1]
	PCA3[,j] = PCA_proj[,2]
	PCA4[,j] = PCA_proj[,3]
		
}


PC1_med <- doppler_phase_comp_fix
PC2_med <- PC2_fix
PC3_med <- PC3_fix
PC4_med <- PC4_fix

PC1_sd <- apply(doppler_phase_comp,1,sd)
PC2_sd <- apply(PCA2,1,sd)
PC3_sd <- apply(PCA3,1,sd)
PC4_sd <- apply(PCA4,1,sd)


PC_coords <- cbind(PC1_med,PC2_med,PC3_med,PC4_med)
PC_SDs <- cbind(PC1_sd,PC2_sd,PC3_sd,PC4_sd)
phases <- spot_phase

colnames(PC_coords) <- c("PC1","PC2","PC3","PC4")
colnames(PC_SDs) <- c("SD1","SD2","SD3","SD4")


  #quartz()
 pdf("old_clean_basis.pdf",width=8,height=8)
  par(mfrow=c(2,2))
  plot(spot_phase,doppler_phase_comp[,1],pch=19,col="blue", main="", ylab="Apparent RV", cex=0.3, xlab="phase", type="l", ylim=range(doppler_phase_comp))
  for(i in 2:n_boot){
 	lines(spot_phase,doppler_phase_comp[,i],pch=19,col="blue")
 }
 points(spot_phase,PC1_med,pch=19,col="red", cex=0.75)
 plot(spot_phase,PCA2[,1],pch=19,col="blue", main="", ylab="PC 1", cex=0.3, xlab="phase", type="l", ylim=range(PCA2))
 for(i in 2:n_boot){
 	lines(spot_phase,PCA2[,i],pch=19,col="blue")
 }
 points(spot_phase,PC2_med,pch=19,col="red", cex=0.75)
  plot(spot_phase,PCA3[,1],pch=19,col="blue", main="", ylab="PC 2", cex=0.3, xlab="phase", type="l", ylim=range(PCA3))
 for(i in 2:n_boot){
 	lines(spot_phase,PCA3[,i],pch=19,col="blue")
 }
 points(spot_phase,PC3_med,pch=19,col="red", cex=0.75)
  plot(spot_phase,PCA4[,1],pch=19,col="blue", main="", ylab="PC 3", cex=0.3, xlab="phase", type="l", ylim=range(PCA4))
 for(i in 2:n_boot){
 	lines(spot_phase,PCA4[,i],pch=19,col="blue")
 }
 points(spot_phase,PC4_med,pch=19,col="red", cex=0.75)
 dev.off()


write.csv(PC_coords,file="PCcoords_incl90_lat40_brightness_correct_HIGH_DENSITY_CLEAN.csv",row.names=FALSE)
write.csv(PC_SDs,file="PCstdevs_incl90_lat40_brightness_correct_HIGH_DENSITY_CLEAN.csv",row.names=FALSE)
write.csv(phases,file="phases_incl90_lat40_brightness_correct_HIGH_DENSITY_CLEAN.csv",row.names=FALSE)





# ####
# pdf("PCA_spec_4512_to_4518_RENORMALIZED_v2.pdf",width=8,height=6.18)
# wv_ind <- wavelengths > 4512 & wavelengths < 4518
# #pdf("vectors1.pdf",width=8,height=6)
# par(mfrow=c(5,1))
# par(cex.axis=1,cex.lab=1.2)
# par(mar=c(0,6.1,0,0),oma=c(4,0.5,0.5,0.5), las=1)
# par(mgp=c(4,0.75,0))
# plot(wavelengths[wv_ind],spot_data_renorm[1,wv_ind],main="",xlab="",ylab="spectrum", col="black",xaxt='n', type="l")
# plot(wavelengths[wv_ind],doppler_comp_simple[wv_ind],main="",xlab="",ylab="doppler\n component", col="black",xaxt='n',type="l")
# #par(mar=c(0,4.1,0,2.1))
# plot(wavelengths[wv_ind],PCA_fix$rotation[wv_ind,1],xlab="",ylab="PC 1",col="black",xaxt='n', type="l")
# plot(wavelengths[wv_ind],PCA_fix$rotation[wv_ind,2],xlab="",ylab="PC 2",col="black",xaxt='n', type="l")
# plot(wavelengths[wv_ind],PCA_fix$rotation[wv_ind,3],xlab="wavelength",ylab="PC 3",col="black", type="l")
# mtext("wavelength (Angstroms)", side=1,outer=FALSE,cex=0.8,line=2.5)
# dev.off()
