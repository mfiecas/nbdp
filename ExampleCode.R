library(RcppArmadillo)
library(Rcpp)
library(astsa)
library(beyondWhittle)
sourceCpp('C:/Users/Brian/Google Drive/PhD Work/LongitudinalEEG/NDPrcpp.cpp')

################   Fit BDP Model to Sunspot Data   ####################
#Read in sunspot data and standardize to mean 0 unit variance
sunspot.norm=(sqrt(sunspot.year)-mean(sqrt(sunspot.year)))/sd(sqrt(sunspot.year))
#Fit BDP model
Nsamp=1000
burnin=100
#see the SpectralBDP function in NDPrcpp.cpp for complete list of settings
Fit1 = SpectralBDP(sunspot.norm,Nsamp)
#calculate curves from posterior samples
Fit1_curves=GetCurvesBDP(Fit1$omegas, Fit1$Lambda, Fit1$P, Fit1$Z, 0.5)
#plot pointwise confidence bands and estimate of log spectral density
plot(Fit1$omegas*2*pi,log(apply(Fit1_curves[burnin:Nsamp,],2,median)), col=1, type="l",
     xlab="Frequency",ylab="Power", main="Estimated Spectral Density")
lines(Fit1$omegas*2*pi,log(apply(Fit1_curves[burnin:Nsamp,],2,quantile,probs=0.025)), lty=2)
lines(Fit1$omegas*2*pi,log(apply(Fit1_curves[burnin:Nsamp,],2,quantile,probs=0.975)), lty=2)




#function to calculate true spectral density function for AR process 
ar_spec = function(omegas, phi, sigma){
  n_omega = length(omegas)
  out = numeric(n_omega)
  p = length(phi)
  for (i in 1:n_omega) {
    temp = 0.0
    for (j in 1:p) {
      temp = temp + phi[j]*(exp(complex(real = 0, imaginary = -2*pi*omegas[i]))^j)
    }
    out[i] = (sigma^2)/(abs((1 - temp))^2)
  }
  return(out)
}

#############Simulate data and test NBDP model
#Time series length
T=100        
#specify subjects per group and group AR parameters for simulated data
N1=5;arG1=c(0.3);
N2=5;arG2=c(0.3,-0.3);
N3=5;arG3=c(-0.3);
N=N1+N2+N3
#simulate time series
ARs=list()
TimeSeries=matrix(NA,T,N)
for(i in 1:N1){ARs[[i]]=arG1}
for(i in (N1+1):(N1+N2)){ARs[[i]]=arG2}
for(i in (N1+N2+1):N){ARs[[i]]=arG3}
for(i in 1:N){
  Temp=arima.sim(list(ar=ARs[[i]]), n=T)
  TimeSeries[,i]=(Temp-mean(Temp))/sd(Temp)
}

#Fit NBDP model
Nsamp = 1000
burnin = 500
#see the SpectralNBDP function in NDPrcpp.cpp for complete list of settings
NBDPFit = SpectralNBDP(X=TimeSeries, Nsamp=Nsamp, K=10, Sup=10)
#Second function takes last iteration from model output and runs Nsamp more iterations
NBDPFit = SpectralNBDP2(NBDPFit, Nsamp, T, Sup=10)
par(mfrow=c(3,5))
for(n in 1:N){
  Curves=GetCurvesNBDP(NBDPFit$omegas, NBDPFit$zeta[,n], NBDPFit$Lambda, NBDPFit$P, 
                       NBDPFit$Z, MaxFreq=0.5)
  plot(NBDPFit$omegas,apply(Curves[burnin:Nsamp,],2,median),type="l", ylim=c(0,2))
  true_curves = ar_spec(NBDPFit$omegas,ARs[[n]],1)
  lines(NBDPFit$omegas, true_curves/mean(true_curves), col=2)
}

##calculate heritability of data comes from twin pairs
#MZ_t1_ind and MZ_t2_ind contain the indices (from 0 to N-1) for the first and second twin in each pair
  #i.e. MZ_t1_ind[i] and MZ_t2_ind[i] are the i-th MZ twin pair
#DZ_t1_ind and DZ_t2_ind are equivalent to the MZ version except they contain the DZ twin pair indices
#indices should be between 0 and N-1 where there are N subjects
#Note that you need a sample much bigger than 15 for this to produce output
MZ_t1_ind = c()
MZ_t2_ind = c()
DZ_t1_ind = c()
DZ_t2_ind = c()
Herit = GetHeritNBDP(NBDPFit, MaxFreq=0.5, burnin=burnin, MZ_t1_ind, MZ_t2_ind, 
                     DZ_t1_ind, DZ_t2_ind, log_spec=FALSE)
par(mfrow=c(1,1))
plot(NBDPFit$omegas,ifelse(apply(Herit$HeritCurve,2,median)<1,apply(Herit$HeritCurve,2,median),1),type="l", ylab="Heritability",
     xlab="Frequency (Hz)",main="Power Spectrum Heritability",ylim=c(0,1))
lines(omegas,ifelse(apply(Herit$HeritCurve,2,quantile,probs=c(0.025))<1,apply(Herit$HeritCurve,2,quantile,probs=c(0.025)),1),lty=2)
lines(omegas,ifelse(apply(Herit$HeritCurve,2,quantile,probs=c(0.975))<1,apply(Herit$HeritCurve,2,quantile,probs=c(0.975)),1),lty=2)




