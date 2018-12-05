#include <iostream>
#include <string>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


//Function to calculate the estimated spectral density curve for a given set of parameters
// [[Rcpp::export]]
vec Get_f(const vec& omegas, const vec& P, const vec& Z, int Lambda, 
          int L, double MaxFreq, int N_freq){
  //omegas is a vector of fourier frequencies
  //P is an L length vector of BDP weights
  //Z is an L length vector of BDP atoms
  //Lambda is the number of Bernstein polynomial components to use
  //L is the truncation level of the BDP
  //MaxFreq is the maximum frequency to consider for the spectral density
  //N_freq is the length of omegas (i.e. number of Fourrier frequencies)
  
  vec f(N_freq, fill::zeros);
  vec weights(Lambda, fill::zeros);
  
  //calculate Bernstein polynomial components weights from BDP weights
  for(int l=0; l<L; l++){
    weights(floor(Z(l)*Lambda)) += P(l);
  }
  
  //calculate estimated spectral density given parameters and add tiny number to avoid errors from rounding to 0
  for(int lambda=0; lambda<Lambda; lambda++){
    if(weights(lambda) > 0){
      for(int n_omega=0; n_omega<N_freq; n_omega++){
        f(n_omega) += weights(lambda)*R::dbeta(omegas(n_omega)/MaxFreq,lambda+1,Lambda-lambda,false) + pow(10,-323);
      }
    }
  }      
  //return the estimated spectral density curve
  return f;
}


//Function to calculate the Whittle likelihood for a given curve and periodogram
// [[Rcpp::export]]
double Whittle(const vec& f, const vec& Pgram){
  //f is a vector containing the estimated spectral density curve
  //Pgram is a vector of corresponding periodogram values
  return sum(-log(f) - Pgram/f);
}


//Function to calculate BDP weights from stick break proportions V for truncation level L
// [[Rcpp::export]]
vec GetWeights(const vec& V, int L){
  //V is vector of stick breaking proportions
  //L is the length of V
  vec CumProds=cumprod(1-V);
  vec out(L);
  out(0)=V(0);
  out.tail(L-1)=CumProds.head(L-1)%V.tail(L-1);
  return out;
}


//Function to calculate stick breaking proportions from BDP weights for truncation level L
// [[Rcpp::export]]
vec GetWeightsInv(const vec& P, int L){
  //P is vector of BDP weights
  //L is the length of P
  vec out(L);
  out(0)=P(0);
  out(L-1)=1;
  vec CumSum = cumsum(1-P);
  for(int l=1; l<L-1; l++){
    out(l) = P(l)/CumSum(l-1);
  }
  return out;
}


//Function to propose new element of V (weights) or Z (atoms) using uniform proposal
// [[Rcpp::export]]
vec NewVZ(const vec& VZ, const vec& epsilon, double dim){
  //VZ is vector of current weights or atoms
  //epsilon is a vector of uniform halfwidths for proposals
  //dim is the element of the weight or atom vector which is being proposed
  
  vec VZnew=VZ;
  //propose new value and take modulus to ensure it falls between 0 and 1
  double prop=fmod(VZ(dim)+epsilon(dim)*(2*randu()-1),1);
  if(prop < 0){prop = prop+1;}
  VZnew(dim) = prop;
  return VZnew;
}


//Function to propose MH move distance for V or Z using mixture of uniforms proposal
// [[Rcpp::export]]
double VZmove(const vec& epsilon,double dim, double MS_mixP){
  //epsilon is a vector of uniform halfwidths for proposals
  //dim is the element of the weight or atom vector which is being proposed
  //MS_mixP is the mixture weight for the unif(0,1) component
  
  if(randu() < MS_mixP){
    return randu()-0.5;
  }else{
    return 2*epsilon(dim)*(randu()-0.5);
  }
}


//function to generate a multivariate normal random sample
// [[Rcpp::export]]
mat mvrnormArma(int n, vec mu, mat sigma) {
  //n is the number of samples to draw
  //mu is the mean vector
  //sigma is the covariance matrix
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}


//function to calculate likelihood of atom coefficients for DBDP
// [[Rcpp::export]]
double Zcoef_lik(int L, const mat& Zcoef_cur, const vec& Zcoef_mu,
                 const mat& Zcoef_SInv){
  //L is the DBDP truncation level
  //Zcoef_cur is the atom coefficients for which we want the likelihood
  //Zcoef_mu is the mean of the posterior normal for the atom coefficients
  //Zcoef_SInv is the inverse covariance of the posterior normal for the atom coefficients
  double likeli=0;
  mat temp(1,1,fill::zeros);
  for(int l=0; l<L; l++){
    temp = (Zcoef_cur.row(l)-Zcoef_mu.t())*Zcoef_SInv*(Zcoef_cur.row(l).t()-Zcoef_mu);
    likeli = likeli + temp(0,0);
  }
  likeli=-0.5*likeli;
  return likeli;
}


//Turn DBDP atom regression coefficients into atoms (Z)
// [[Rcpp::export]]
mat ZcoefToZ(const mat& Y, const mat& Zcoef_cur, int L, int N){
  //Y is a matrix of covariates with a row for each time series
  //Zcoef_cur is a matrix of coefficients with predictors across the rows and atoms down the rows
  //L is the DBDP truncation level
  //N is the number of time series
  mat Z(L, N);
  vec temp(L);
  for(int n=0; n < N; n++){
    temp = exp(Zcoef_cur*Y.row(n).t());
    Z.col(n) = temp/(ones(L)+temp);
  }
  return Z;
}


//functions to sample from wishart and inverse wishart distributions
// [[Rcpp::export]]
mat rwishart(int df, const mat& S){
  // Dimension of returned wishart
  int m = S.n_rows;
  // Z composition: sqrt chisqs on diagonal, random normals below diagonal, misc above diagonal
  mat Z(m,m);
  // Fill the diagonal
  for(int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i));
  }
  // Fill the lower matrix with random guesses
  for(int j = 0; j < m; j++){
    for(int i = j+1; i < m; i++){
      Z(i,j) = R::rnorm(0,1);
    }
  }
  // Lower triangle * chol decomp
  mat C = trimatl(Z).t() * chol(S);
  // Return random wishart
  return C.t()*C;
}
// [[Rcpp::export]]
mat riwishart(int df, const mat& S){
  return rwishart(df,S.i()).i();
}


//function to calculate multivariate normal densities
// [[Rcpp::export]]
double dmvnrm_arma(rowvec x, rowvec mean, mat sigma){
  int xdim = x.n_cols;
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log(2.0 * M_PI);
  vec z = rooti*trans(x - mean);
  return(constants - 0.5 * sum(z%z) + rootisum);
}




//function to sample number of Bernstein polynomials, Lambda, using MH step
// [[Rcpp::export]]
List sample_Lambda(int N, int L, const vec& omegas, const mat& Pgrams, const vec& P_cur, 
                   const mat& Z_cur, int Lambda_cur, string Lambda_prior, int Lambda_init, 
                   int Lambda_max, double MaxFreq, string Lambda_prop_dist, 
                   const mat& f_cur, int N_freq, bool Dep){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a matrix containing the latest iteration BDP atoms with a column for each subject in DDP
  //Lambda_cur is the current number of Bernstein polynomial components
  //Lambda_prior is the prior for the number of Bernstein polynomial components
  //Lambda_init is the initial number of Bernstein polynomial components in the MCM chain
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //MaxFreq is the maximum frequency to consider for the spectral density
  //Lambda_prop_dist is the proposal distribution to be used for Lambda, up1down1 or poisson. If incorrectly/not specified Lambda will be constant
  //f_cur is a matrix containing the current estimated spectral density curve with a column/subject in DDP
  //N_freq is the number of Fourier frequencies
  //Dep is a boolean to determine if DDP or DP is being used
  
  int Lambda_new;
  int Lambda_prop = Lambda_cur;
  double accept_Lambda = 0;
  //propose new value for Lambda using poisson proposal and account for proposal density in MH prob
  if(Lambda_prop_dist=="up1down1"){
    if(Lambda_cur==1){Lambda_prop=2;}else{
      if(Lambda_cur==Lambda_max){Lambda_prop=Lambda_max-1;}else{
        if(randu()<0.5){Lambda_prop=Lambda_cur-1;}else{Lambda_prop=Lambda_cur+1;}
      }
    }
  }
  if(Lambda_prop_dist=="poisson"){
    Lambda_prop=R::rpois(Lambda_cur);
    if(Lambda_prop==0){Lambda_prop=1;}else{if(Lambda_prop>Lambda_max){Lambda_prop=Lambda_max;}}
    accept_Lambda = R::dpois(Lambda_cur,Lambda_prop,true) - 
                    R::dpois(Lambda_prop,Lambda_cur,true);
  }
  
  //account for prior and likelihood contributions to the MH prob
  mat f_new(f_cur.n_rows, f_cur.n_cols);
  mat f_prop(f_cur.n_rows, f_cur.n_cols);
  if(Lambda_prop==Lambda_cur){
    Lambda_new = Lambda_cur;
    f_new = f_cur;
  }else{
    if(Dep == false){
      f_prop.col(0) =  Get_f(omegas, P_cur, Z_cur.col(0), Lambda_prop, L, MaxFreq, N_freq);
    }else{
      for(int i=0; i<N; i++){
        f_prop.col(i) =  Get_f(omegas, P_cur, Z_cur.col(i), Lambda_prop, L, MaxFreq, N_freq);
      }
    }
    if(Lambda_prior=="flat"){accept_Lambda += 0;}
    if(Lambda_prior=="poisson"){accept_Lambda += R::dpois(Lambda_prop,Lambda_init,true) - 
                                                 R::dpois(Lambda_cur,Lambda_init,true);}
    if(Lambda_prior=="expon"){accept_Lambda += -0.05*pow(Lambda_prop,2) + 
                                                0.05*pow(Lambda_cur,2);}

    if(N != 0){
      for(int i=0; i<N; i++){
        if(Dep == false){
          accept_Lambda += Whittle(f_prop.col(0), Pgrams.col(i)) - 
                           Whittle(f_cur.col(0), Pgrams.col(i));
        }else{
          accept_Lambda += Whittle(f_prop.col(i), Pgrams.col(i)) - 
                           Whittle(f_cur.col(i), Pgrams.col(i));
        }
      }
    }
    
    //accept or reject proposed Lambda
    if(accept_Lambda > log(randu())){
      Lambda_new = Lambda_prop;
      f_new = f_prop;
    }else{
      Lambda_new = Lambda_cur;
      f_new = f_cur;
    }
  }  
  
  return List::create(
      _["Lambda_new"] = Lambda_new,
      _["f_new"] = f_new
  );
}      


//function to sample stick breaking weights, V, using MH step and uniform proposal
// [[Rcpp::export]]
List sample_V(int N, int L, const vec& omegas, const mat& Pgrams, const vec& V_cur, 
              const vec& P_cur, const mat& Z_cur, int Lambda_cur, const vec& epsilon, 
              double alpha_L_cur, double MaxFreq, const mat& f_cur, int N_freq, bool Dep){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //V_cur is a vector containing the latest iteration BDP stick breaking weights
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a matrix containing the latest iteration BDP atoms with a column/subject in DDP
  //Lambda_cur is the current number of Bernstein polynomial components
  //epsilon is a vector of uniform halfwidths for proposals
  //alpha_L_cur is the latest value of the alpha_L concentration parameter
  //MaxFreq is the maximum frequency to consider for the spectral density
  //f_cur is a matrix containing the current spectral density estimate w/ column/subject in DDP
  //N_freq is the # of Fourier frequencies
  //Dep is a boolean designating if DDP or DP is being used

  vec V_prop(L);
  vec P_prop(L);
  mat f_prop(N_freq, f_cur.n_cols);
  vec V_new = V_cur;
  vec P_new = P_cur;
  mat f_new = f_cur;
  //calculate likelihood of current weights and curves
  double V_new_lik = (alpha_L_cur-1)*sum(log(1-V_new.head(L-1)));
  if(N != 0){
    for(int i=0; i<N; i++){     
      if(Dep == false){
        V_new_lik += Whittle(f_new.col(0), Pgrams.col(i));
      }else{
        V_new_lik += Whittle(f_new.col(i), Pgrams.col(i));
      }
    }
  }
  double V_prop_lik;
  double accept_V;  
  //loop through updating each of the L-1 random values of stick breaking weights V
  for(int l=0; l<(L-1); l++){
    //propose new V for l-th dimension
    V_prop = NewVZ(V_new, epsilon, l);
    P_prop = GetWeights(V_prop, L);  
    //get proposed curves given proposed Ps
    if(Dep == false){
      f_prop.col(0) = Get_f(omegas, P_prop, Z_cur.col(0), Lambda_cur, L, MaxFreq, N_freq);
    }else{
      for(int i=0; i<N; i++){
        f_prop.col(i) = Get_f(omegas, P_prop, Z_cur.col(i), Lambda_cur, L, MaxFreq, N_freq);
      }
    }
    
    //calculate acceptance probability
    V_prop_lik = (alpha_L_cur-1)*sum(log(1-V_prop.head(L-1)));
    if(N != 0){
      for(int i=0; i<N; i++){
        if(Dep == false){
          V_prop_lik += Whittle(f_prop.col(0), Pgrams.col(i));
        }else{
          V_prop_lik += Whittle(f_prop.col(i), Pgrams.col(i));
        }
      }
    }
    accept_V = V_prop_lik - V_new_lik;
    //accept or reject new V value using MH step
    if(accept_V > log(randu())){
      V_new = V_prop;
      P_new = P_prop;
      V_new_lik = V_prop_lik;
      f_new = f_prop;
    }
  }
  return List::create(
    _["V_new"] = V_new,
    _["P_new"] = P_new,
    _["f_new"] = f_new
  ); 
}


//function to sample atoms, Z, using MH step and uniform proposal
// [[Rcpp::export]]
List sample_Z(int N, int L, const vec& omegas, const mat& Pgrams, const vec& P_cur, 
             const vec& Z_cur, int Lambda_cur, const vec& epsilon, double MaxFreq, 
             const vec& f_cur, int N_freq){
  //N is the number of subjects
  //L is BDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //P_cur is a vector containing the latest iteration BDP weights
  //Z_cur is a vector containing the latest iteration BDP atoms
  //Lambda_cur is the current number of Bernstein polynomial components
  //epsilon is a vector of uniform halfwidths for proposals
  //MaxFreq is the maximum frequency to consider for the spectral density
  //f_cur is a vector containing the current spectral density estimate
  //N_freq is the # of Fourier frequencies

  vec Z_prop(L);
  vec Z_new = Z_cur;
  vec f_prop(N_freq);
  vec f_new = f_cur;
  double Z_new_lik = 0;
  if(N != 0){
    for(int i=0; i<N; i++){
      Z_new_lik += Whittle(f_cur, Pgrams.col(i));
    }
  }

  double Z_prop_lik;
  double accept_Z;  
  //loop through updating each of the L BDP atoms, Z
  for(int l=0; l<L; l++){
    Z_prop_lik = 0;
    //propose new l-th atom
    Z_prop = NewVZ(Z_new, epsilon, l);
    f_prop = Get_f(omegas, P_cur, Z_prop, Lambda_cur, L, MaxFreq, N_freq);
    //calculate acceptance probability
    if(N != 0){
      for(int i=0; i<N; i++){
        Z_prop_lik += Whittle(f_prop, Pgrams.col(i));
      }
    }
    accept_Z = Z_prop_lik - Z_new_lik;
    //acceptance or reject using MH step
    if(accept_Z > log(randu())){
      Z_new = Z_prop; 
      Z_new_lik = Z_prop_lik;
      f_new = f_prop; 
    }
  }  
  return List::create(
      _["Z_new"] = Z_new,
      _["f_new"] = f_new
  );
}


//function to sample coefficients for atoms, Z, from the DBDP
// [[Rcpp::export]]
List sample_Zcoef(int N, int L, const vec& omegas, const mat& Pgrams, const mat& Y,
                  const vec& P_cur, const mat& Zcoef_cur, const mat& Z_cur,
                  int Lambda_cur, const vec& Zcoef_mu_cur, const mat& Zcoef_SInv_cur,
                  const mat& Zcoef_prop_sig, double MaxFreq, const mat& f_cur, int N_freq){
  //N is the number of subjects
  //L is DBDP truncation level
  //omegas is a vector of fourier frequencies
  //Pgrams is a matrix of corresponding periodogram values with a column for each subject
  //Y is a matrix of covariates with a row for each time series
  //P_cur is a vector containing the latest iteration BDP weights
  //Zcoef_cur is a matrix containing the latest iteration DBDP atom coefficients
  //Z_cur is a matrix containing the latest iteration DBDP atoms with a column for each subject
  //Lambda_cur is the current number of Bernstein polynomial components
  //Zcoef_mu_cur is the posterior mean of the DBDP atom coefficients
  //Zcoef_SInv_cur is the posterior covariance of the DBDP atom coefficients
  //Zcoef_prop_sigma is the covariance matrix of the DBDP atom coefficient proposal distribution
  //MaxFreq is the maximum frequency to consider for the spectral density
  //f_cur is a matrix containing the current estimated spectral densities where each column is a subject
  //N_freq is the number of Fourier frequencies

  mat Zcoef_new = Zcoef_cur;
  mat Zcoef_prop = Zcoef_cur;
  mat Z_new = Z_cur;
  mat Z_prop = Z_cur;
  mat f_new = f_cur;
  mat f_prop(N_freq,N);
  double Z_new_lik = Zcoef_lik(L, Zcoef_new, Zcoef_mu_cur, Zcoef_SInv_cur);
  if(N != 0){
    for(int n=0; n<N; n++){
      Z_new_lik += Whittle(f_new.col(n), Pgrams.col(n));
    }
  }
  double Z_prop_lik;
  double accept_Z;
  for(int l=0; l<L; l++){
    Zcoef_prop.row(l) = mvrnormArma(1, Zcoef_new.row(l).t(), Zcoef_prop_sig);
    Z_prop = ZcoefToZ(Y, Zcoef_prop, L, N);
    Z_prop_lik = Zcoef_lik(L, Zcoef_prop, Zcoef_mu_cur, Zcoef_SInv_cur);
    if(N != 0){
      for(int n=0; n<N; n++){
        f_prop.col(n) = Get_f(omegas, P_cur, Z_prop.col(n), Lambda_cur, L, MaxFreq, N_freq);
        Z_prop_lik += Whittle(f_prop.col(n), Pgrams.col(n));
      }
    }
    accept_Z = Z_prop_lik - Z_new_lik;
    if(accept_Z > log(randu())){
      Zcoef_new = Zcoef_prop;
      Z_new = Z_prop;
      Z_new_lik = Z_prop_lik;
      f_new = f_prop;
    }
    Zcoef_prop = Zcoef_new;
    Z_prop = Z_new;
  }
  return List::create(
    _["Zcoef_new"] = Zcoef_new,
    _["Z_new"] = Z_new,
    _["f_new"] = f_new
  );
}


//function to sample cluster indicators, zeta
// [[Rcpp::export]]
int sample_zeta(int K, const vec& pi_cur, const vec& Pgram, const mat& f_cur, int N_freq){
  //K is the truncation level for clustering DP (top-level)
  //pi_cur is vector of clustering DP weights
  //Pgram is a vector of corresponding periodogram values for a subject
  //f_cur contains the estimated spectral density curves w/ a column for each cluster
  //N_freq is the number of Fourier frequencies
  
  int zetanew;
  vec zetaPs(K);
  double temprand;
  vec zetaPmax(K);
  vec logPi = log(pi_cur);
  
  //calculate multinomial weights for each subject for each cluster
  zetaPs.fill(0);
  for(int k=0; k < K; k++){
    zetaPs(k) = logPi(k) + Whittle(f_cur.col(k), Pgram);
  }
  zetaPmax.fill(max(zetaPs));
  zetaPs = exp(zetaPs - zetaPmax);
  
  //sample from multinomial dist to get new cluster assignments for each subject
  temprand=randu();
  vec cumProbs=cumsum(zetaPs/sum(zetaPs));
  int check=0;
  int counter=0;
  while(check == 0){
    if(temprand < cumProbs(counter)){zetanew=counter;check=1;}
    counter++;
  }
  return zetanew;
}


//MCMC Sampling Algorithm for single subject spectral BDP
// [[Rcpp::export]]
List SpectralBDP(const vec& X, int Nsamp, const vec& eps=0, int L=20, 
                 double SamplingRate=1, double MaxFreq=0.5, int Sup=500, 
                 int Lambda_init=30, string Lambda_prior="flat", int Lambda_max=300, 
                 string Lambda_prop_dist="poisson", bool alpha_const=true, 
                 double alpha_init=1, double a_alpha_L=1, double b_alpha_L=1){
  //X is a vector containing the time series to be analyzed
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
    //either "expon", "poisson", or "flat"
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prop_dist is the proposal distribution to be used for Lambda MH steps
    //either "poisson" or "up1down1". If mispecified or "const" then lambda will be constant 
  //alpha_const is a boolean to determine if the concetration parameter is set to constant
  //alpha_init is the initial value for concentration parameter and the value for all samples if constant
  //a_alpha_L and b_alpha_L are the gamma hyperparameters to be used if alpha is to be sampled

  //determine Fourier frequencies to be considered
  double T = X.n_elem;
  int N_freq_temp = floor((T-1)/2);
  int N_freq = 0;
  vec omegas(N_freq_temp);
  for(int n_freq=0; n_freq < N_freq_temp; n_freq++){
    omegas(n_freq) = (n_freq+1)/T*SamplingRate;
    if((n_freq+1)/T*SamplingRate < MaxFreq){N_freq += 1;}
  }
  omegas=omegas.head(N_freq);
  
  //calculate periodograms
  vec Pgram=square(abs(fft(X)))/T;
  Pgram=Pgram.head(N_freq);

  vec epsilon(L);
  if(eps(0)==0){
    for(int l=0; l<L; l++){
      epsilon(l) = (l+1)/(l+1+2*sqrt(T));
    }
  }else{
    epsilon = eps;
  }
  
  //initialize concentration parameter alpha
  vec alpha_L(Nsamp);
  alpha_L(0)=alpha_init;

  //initialize BDP (from Choudhuri)
  vec Lambda(Nsamp);
  Lambda(0)=Lambda_init;
  List LambdaList(2);

  mat V(Nsamp, L);
  V.submat(0,0,0,L-2).fill(0.5);
  V(0,L-1)=1;

  mat P(Nsamp, L);
  P.row(0)=GetWeights(V.row(0).t(), L).t();
  List VPList(3);

  mat Z(Nsamp, L);
  Z.row(0).randu();
  List ZList(2);
  
  vec f = Get_f(omegas, P.row(0).t(), Z.row(0).t(), Lambda(0), L, MaxFreq, N_freq);

  for(int iter=1; iter<Nsamp; iter++){
    //Update number of Bernstein polynomials, Lambda
    LambdaList = sample_Lambda(1, L, omegas, Pgram, P.row(iter-1).t(), Z.row(iter-1).t(), 
                               Lambda(iter-1), Lambda_prior, Lambda_init, Lambda_max, 
                               MaxFreq, Lambda_prop_dist, f, N_freq, false);
    Lambda(iter) = as<int>(LambdaList("Lambda_new"));
    f = as<vec>(LambdaList("f_new"));

    //Update weights V and P
    VPList = sample_V(1, L, omegas, Pgram, V.row(iter-1).t(), P.row(iter-1).t(), 
                      Z.row(iter-1).t(), Lambda(iter), epsilon, alpha_L(iter-1), 
                      MaxFreq, f, N_freq, false);
    V.row(iter) = as<vec>(VPList("V_new")).t();
    P.row(iter) = as<vec>(VPList("P_new")).t();
    f = as<vec>(VPList("f_new"));

    //Update Z
    ZList = sample_Z(1, L, omegas, Pgram, P.row(iter).t(), Z.row(iter-1).t(), 
                     Lambda(iter), epsilon, MaxFreq, f, N_freq);
    Z.row(iter) = as<vec>(ZList("Z_new")).t();
    f = as<vec>(ZList("f_new"));
  
    //Update alpha_L via Gibbs step
    if(alpha_const==false){
      alpha_L(iter)=R::rgamma(a_alpha_L+L-1, b_alpha_L - sum(log(V.row(iter).head(L-1))));
    }
    if(alpha_const==true){
      alpha_L(iter)=alpha_L(iter-1);
    }
 
    if((iter-1)/Sup!=iter/Sup){cout << iter << " posterior samples complete" << endl;}
  }

  return List::create(
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["V"] = V,
    _["Z"] = Z,
    _["omegas"] = omegas,
    _["Pgram"] = Pgram
  ); 
}  


//MCMC Sampling Algorithm for DDP, all subjects in 1 cluster
// [[Rcpp::export]]
List SpectralDBDP(const mat& X, const mat& Y, int Nsamp, double Zcoef_prop_sig=0.2, 
                  const vec& eps=0, int L=20, double SamplingRate=1, 
                  double MaxFreq=0.5, int Sup=500, int Lambda_init=30, int Lambda_max=300,
                  string Lambda_prior="flat", string Lambda_prop_dist="poisson",
                  bool alpha_L_const=true, double a_alpha_L=1, double b_alpha_L=1,
                  double alpha_L_init=1){
  //X is a matrix where each column is a time series to be analyzed
  //Y is a matrix of covariates with a row for each time series
  //Nsamp is the number of posterior samples
  //L is the truncation level for the DBDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //Zcoef_prop_sig is the variance to be used for the coefficient proposal draws
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
    //either "expon", "poisson", or "flat"
  //Lambda_prop_dist is the proposal distribution to be used for Lambda MH steps
    //either "poisson" or "up1down1". If mispecifiec or "const" then lambda will be constant
  //alpha_L_const is a boolean to determine if the DBDP concetration parameter is set to constant
  //alpha_L_init is the initial value for DBDP concentration parameter and the value for all samples if constant
  //a_alpha_L and b_alpha_L are the gamma hyperparameters to be used if alpha_L is to be sampled
  
  //calculate Fourier frequnecies to be used
  int N = X.n_cols;
  double T = X.n_rows;
  int N_freq_temp = floor((T-1)/2);
  int N_freq = 0;
  int N_pred = Y.n_cols;
  vec omegas(N_freq_temp);
  for(int n_freq=0; n_freq < N_freq_temp; n_freq++){
    omegas(n_freq) = (n_freq+1)/T*SamplingRate;
    if(omegas(n_freq) < MaxFreq){N_freq += 1;}
  }
  omegas=omegas.head(N_freq);
  
  //calculate periodograms
  mat Pgrams(N_freq, N);
  vec PgramTemp(T);
  for(int n=0; n < N; n++){
    PgramTemp = square(abs(fft(X.col(n))))/T;
    Pgrams.col(n) = PgramTemp.head(N_freq);
  }
  
  vec epsilon(L);
  if(eps(0)==0){
    for(int l=0; l<L; l++){
      epsilon(l) = (l+1)/(l+1+2*sqrt(T));
    }
  }else{
    epsilon = eps;
  }
  
  vec prop_sig_vec(N_pred);
  prop_sig_vec.fill(Zcoef_prop_sig);
  mat Zcoef_prop_sig_mat = diagmat(prop_sig_vec);
  
  //initiate structures to hold posterior samples for parameters
  vec alpha_L(Nsamp);
  alpha_L(0)=alpha_L_init;
  vec Lambda(Nsamp);
  mat V(Nsamp, L);
  mat P(Nsamp, L);
  cube Zcoef(L, N_pred, Nsamp);
  mat Z(L, N);
  mat f(N_freq, N);
  
  //initialize parameter estimates randomly for each cluster if no warm start
  Lambda(0)=Lambda_init;
  List LambdaList(2);
  V.submat(0,0,0,L-2).fill(0.5);
  V.col(L-1).fill(1);
  P.row(0)=GetWeights(V.row(0).t(), L).t();
  List VPList(3);
  
  Zcoef.slice(0).fill(0);
  Z=ZcoefToZ(Y, Zcoef.slice(0), L, N);
  List ZList(2);
  
  mat Zcoef_mu(Nsamp, N_pred, fill::zeros);
  vec Zcoef_mu_pmu(N_pred);
  mat Zcoef_mu_pS(N_pred, N_pred);
  cube Zcoef_S(N_pred, N_pred, Nsamp);
  Zcoef_S.slice(0).eye();
  mat SInv = inv(Zcoef_S.slice(0));
  mat S0Inv(N_pred, N_pred);
  S0Inv.eye();
  mat IWvar(N_pred, N_pred);
  
  for(int i=0; i<N; i++){
    f.col(i) = Get_f(omegas, P.row(0).t(), Z.col(i), Lambda(0), L, MaxFreq, N_freq);
  }
  
  for(int iter=1; iter<Nsamp; iter++){
    //Update number of Bernstein polynomials, Lambda
    LambdaList = sample_Lambda(N, L, omegas, Pgrams, P.row(iter-1).t(), Z, Lambda(iter-1),
                               Lambda_prior, Lambda_init, Lambda_max, MaxFreq,
                               Lambda_prop_dist, f, N_freq, true);
    Lambda(iter) = as<int>(LambdaList("Lambda_new"));
    f = as<mat>(LambdaList("f_new"));
    
    //Update DBDP stick-breaking weights, V
    VPList = sample_V(N, L, omegas, Pgrams, V.row(iter-1).t(), P.row(iter-1).t(), Z,
                      Lambda(iter), epsilon, alpha_L(iter-1), MaxFreq, f, N_freq, true);
    V.row(iter) = as<vec>(VPList("V_new")).t();
    P.row(iter) = as<vec>(VPList("P_new")).t();
    f = as<mat>(VPList("f_new"));
    
    //Update atom coefficients and their priors
    ZList = sample_Zcoef(N, L, omegas, Pgrams, Y, P.row(iter).t(), Zcoef.slice(iter-1),
                         Z, Lambda(iter), Zcoef_mu.row(iter-1).t(), SInv,
                         Zcoef_prop_sig_mat, MaxFreq, f, N_freq);
    Zcoef.slice(iter) = as<mat>(ZList("Zcoef_new"));
    Z = as<mat>(ZList("Z_new"));
    f = as<mat>(ZList("f_new"));
    
    Zcoef_mu_pS = inv(S0Inv + L*SInv);
    Zcoef_mu_pmu = Zcoef_mu_pS * (SInv*sum(Zcoef.slice(iter),0).t());
    Zcoef_mu.row(iter) = mvrnormArma(1, Zcoef_mu_pmu, Zcoef_mu_pS);
    
    IWvar.eye();
    for(int l=0; l < L; l++){
      IWvar = IWvar + ( (Zcoef_mu.row(iter).t() - Zcoef.slice(iter).row(l).t()) *
        (Zcoef_mu.row(iter) - Zcoef.slice(iter).row(l)) );
    }
    Zcoef_S.slice(iter) = riwishart(L+N_pred, IWvar);
    SInv = inv(Zcoef_S.slice(iter));
    
    //Update alpha_L via Gibbs step
    if(alpha_L_const==false){
      alpha_L(iter)=R::rgamma(a_alpha_L+L-1, b_alpha_L - sum(log(V.row(iter).head(L-1))));
    }
    if(alpha_L_const==true){
      alpha_L(iter)=alpha_L(iter-1);
    }
    
    if((iter-1)/Sup!=iter/Sup){cout << iter << " posterior samples complete" << endl;}
  }
  
  return List::create(
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["V"] = V,
    _["Zcoef"] = Zcoef,
    _["omegas"] = omegas,
    _["Pgrams"] = Pgrams
  );
}


//MCMC Sampling Algorithm for nested bernstein polynomials dirichlet process spectral estimator using split and merge sampler
// [[Rcpp::export]]
List SpectralNBDP(const mat& X, int Nsamp, const vec& eps=0, int L=20, int K=35, double SamplingRate=1, 
                  double MaxFreq=0.5, int Sup=500, int Lambda_init=30, int Lambda_max=300, 
                  string Lambda_prior="flat", double alpha_K=1, double alpha_L=1, 
                  bool warmStart=false, int Nwarm=1000){
  //X is a matrix where each column is a time series to be analyzed
  //Nsamp is the number of posterior samples
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //L is the truncation level for the BDP
  //K is the truncation level for the upper level DP
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
    //either "expon", "poisson", or "flat"
  //alpha_K and alpha_L are the DP concentration parameters
  //warmStart is a boolean to determine if each cluster should be initialed at a fit from the BDP model
  //Nwarm is the number of BDP iterations to run for each cluster when warming up the clusters
  
  //determine Fourier frequencies to consider
  int N = X.n_cols;
  double T = X.n_rows;
  int N_freq_temp = floor((T-1)/2);
  int N_freq = 0;
  vec omegas(N_freq_temp);
  for(int n_freq=0; n_freq < N_freq_temp; n_freq++){
    omegas(n_freq) = (n_freq+1)/T*SamplingRate;
    if(omegas(n_freq) < MaxFreq){N_freq += 1;}
  }
  omegas=omegas.head(N_freq);
  
  //calculate periodograms
  mat Pgrams(N_freq, N);
  vec PgramTemp(T);
  for(int n=0; n < N; n++){
    PgramTemp = square(abs(fft(X.col(n))))/T;
    Pgrams.col(n) = PgramTemp.head(N_freq);
  }
  
  vec epsilon(L);
  if(eps(0)==0){
    for(int l=0; l<L; l++){
      epsilon(l) = (l+1)/(l+1+2*sqrt(T));
    }
  }else{
    epsilon = eps;
  }
  
  //initiate structures to hold posterior samples for parameters
  mat Lambda(Nsamp,K);
  cube P(L,K,Nsamp);
  cube Z(L,K,Nsamp);
  mat zeta(Nsamp, N);
  
  //initialize parameter estimates randomly for each cluster
  mat V(L,K);
  
  if(warmStart==true){   
    uvec perm_id(N);
    for(int i=0; i<N; i++){perm_id(i)=i;}
    perm_id = shuffle(perm_id);
    //initialize parameter estimates randomly for each cluster if warm start
    if(N >= K){
      List SubjFit;
      for(int k=0; k<K; k++){
        cout<<"Warming up cluster "<<k+1<<endl;
        SubjFit = SpectralBDP(X.col(perm_id(k)), Nwarm, epsilon, L, SamplingRate, MaxFreq, 
                              100000, Lambda_init, Lambda_prior, Lambda_max);
        Lambda(0,k) = as<vec>(SubjFit("Lambda"))(Nwarm-1);
        V.col(k) = as<mat>(SubjFit("V")).row(Nwarm-1).t();
        Z.slice(0).col(k) = as<mat>(SubjFit("Z")).row(Nwarm-1).t();
      }
    }else{
      List SubjFit;
      for(int k=0; k<N; k++){
        cout<<"Warming up cluster "<<k+1<<endl;
        SubjFit = SpectralBDP(X.col(perm_id(k)), Nwarm, epsilon, L, SamplingRate, MaxFreq, 
                              100000, Lambda_init, Lambda_prior, Lambda_max);
        Lambda(0,k) = as<vec>(SubjFit("Lambda"))(Nwarm-1);
        V.col(k) = as<mat>(SubjFit("V")).row(Nwarm-1).t();
        Z.slice(0).col(k) = as<mat>(SubjFit("Z")).row(Nwarm-1).t();
      }
      Lambda.row(0).tail(K-N).fill(Lambda_init);
      V.submat(0,N,L-2,K-1).fill(0.5);
      V.submat(L-1,N,L-1,K-1).fill(1);
      Z.slice(0).submat(0,N,L-1,K-1).randu();
    }
  }else{
    //initialize parameter estimates randomly for each cluster if no warm start
    Lambda.row(0).fill(Lambda_init);
    V.submat(0,0,L-2,K-1).fill(0.5);
    V.row(L-1).fill(1);
    Z.slice(0).randu();
  }
  
  for(int k=0; k<K; k++){
    P.slice(0).col(k)=GetWeights(V.col(k), L);
  }
  
  vec u(K);
  vec pis(K);
  u.head(K-1).fill(0.5);
  u(K-1) = 1;
  pis = GetWeights(u, K);
  
  mat f(N_freq, K);
  for(int k=0; k<K; k++){
    f.col(k) = Get_f(omegas, P.slice(0).col(k), Z.slice(0).col(k), Lambda(0,k), L, MaxFreq, N_freq);
  }
  
  vec n(K, fill::zeros); 
  for(int i=0; i<N; i++){
    zeta(0,i) = sample_zeta(K, pis, Pgrams.col(i), f, N_freq);
    n(zeta(0,i)) += 1;
  }   
  field<uvec> inK(K);    
  for(int k=0; k < K; k++){   
    inK(k) = find(zeta.row(0).t() == k);
  }
  
  List LambdaList(2);
  List VPList(3);
  List ZList(2);
  
  for(int iter=1; iter<Nsamp; iter++){
    for(int k=0; k<(K-1); k++){
      u(k) = R::rbeta(1+n(k), alpha_K + sum(n.tail(K-k-1)));
    }
    pis=GetWeights(u, K);  
      
    n.fill(0);
    for(int i=0; i<N; i++){
      zeta(iter,i) = sample_zeta(K, pis, Pgrams.col(i), f, N_freq);
      n(zeta(iter,i)) += 1;
    } 
    
    field<uvec> inK(K);
    for(int k=0; k<K; k++){
      uvec inK2 = find(zeta.row(iter).t() == k);
      inK(k) = inK2;
    }
      
    //Update number of Bernstein polynomials (Lambda) in each cluster
    for(int k=0; k<K; k++){
      LambdaList = sample_Lambda(n(k), L, omegas, Pgrams.cols(inK(k)), P.slice(iter-1).col(k),
                                 Z.slice(iter-1).col(k), Lambda(iter-1,k), Lambda_prior, 
                                 Lambda_init, Lambda_max, MaxFreq, "poisson", f.col(k), 
                                 N_freq, false);
      Lambda(iter,k) = as<int>(LambdaList("Lambda_new"));
      f.col(k) = as<vec>(LambdaList("f_new"));
    }
      
    //Update BDP stick breaking weights, V, for each cluster
    for(int k=0; k<K; k++){
      VPList = sample_V(n(k), L, omegas, Pgrams.cols(inK(k)), V.col(k), 
                        P.slice(iter-1).col(k), Z.slice(iter-1).col(k), Lambda(iter,k),
                        epsilon, alpha_L, MaxFreq, f.col(k), N_freq, false);
      V.col(k) = as<vec>(VPList("V_new"));
      P.slice(iter).col(k) = as<vec>(VPList("P_new"));
      f.col(k) = as<vec>(VPList("f_new"));
    }
      
    //Update BDP atoms, Z, for each cluster
    for(int k=0; k<K; k++){
      ZList = sample_Z(n(k), L, omegas, Pgrams.cols(inK(k)), P.slice(iter).col(k),
                       Z.slice(iter-1).col(k), Lambda(iter,k), epsilon, MaxFreq, f.col(k), N_freq);
      Z.slice(iter).col(k) = as<vec>(ZList("Z_new"));
      f.col(k) = as<vec>(ZList("f_new"));
    }
    
    if((iter+1)%Sup == 0){cout << iter+1 << " posterior samples complete" << endl;}
  }
  
  return List::create(
    _["zeta"] = zeta,
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["Z"] = Z,
    _["omegas"] = omegas,
    _["Pgrams"] = Pgrams
  ); 
}    


//MCMC Sampling Algorithm to continue from previous NBDP fit
// [[Rcpp::export]]
List SpectralNBDP2(const List& Fit_old, int Nsamp, int T, const vec& eps=0, double SamplingRate=1, 
                   double MaxFreq=0.5, int Sup=500, int Lambda_init=30, int Lambda_max=300, 
                   string Lambda_prior="flat", double alpha_K=1, double alpha_L=1){
  //Fit_old is output from SpectralNBDP() or SpectralNBDP2() function
  //Nsamp is the number of posterior samples
  //eps is a vector of uniform halfwidths for V and Z proposals
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
    //either "expon", "poisson", or "flat"
  //alpha_K and alpha_L are the DP concentration parameters

  vec omegas = as<vec>(Fit_old("omegas"));
  mat Pgrams = as<mat>(Fit_old("Pgrams"));
  mat zeta_old = as<mat>(Fit_old("zeta"));
  mat Lambda_old = as<mat>(Fit_old("Lambda"));
  cube P_old = as<cube>(Fit_old("P"));
  cube Z_old = as<cube>(Fit_old("Z"));
  int Nsamp_old = zeta_old.n_rows;
  int N_freq = Pgrams.n_rows;
  int N = Pgrams.n_cols;
  int L = P_old.n_rows;
  int K = P_old.n_cols;
  vec epsilon(L);
  if(eps(0)==0){
    for(int l=0; l<L; l++){
      epsilon(l) = (l+1)/(l+1+2*sqrt(T));
    }
  }else{
    epsilon = eps;
  }

  mat Lambda(Nsamp+Nsamp_old,K);
  cube P(L,K,Nsamp+Nsamp_old);
  cube Z(L,K,Nsamp+Nsamp_old);
  mat zeta(Nsamp+Nsamp_old, N);
  
  Lambda.rows(0,Nsamp_old-1) = Lambda_old;
  P.slices(0,Nsamp_old-1) = P_old;
  Z.slices(0,Nsamp_old-1) = Z_old;
  zeta.rows(0,Nsamp_old-1) = zeta_old;
  
  mat V(L, K);
  for(int k=0; k<K; k++){
    V.col(k) = GetWeightsInv(P.slice(Nsamp_old-1).col(k), L);
  }
  
  mat f(N_freq, K);
  for(int k=0; k<K; k++){
    f.col(k) = Get_f(omegas, P.slice(Nsamp_old-1).col(k), Z.slice(Nsamp_old-1).col(k), 
                     Lambda(Nsamp_old-1,k), L, MaxFreq, N_freq);
  }
  
  vec n(K, fill::zeros); 
  for(int i=0; i<N; i++){
    n(zeta(Nsamp_old-1,i)) += 1;
  }   
  field<uvec> inK(K);    
  for(int k=0; k < K; k++){   
    inK(k) = find(zeta.row(Nsamp_old-1).t() == k);
  }
  
  List LambdaList(2);
  List VPList(3);
  List ZList(2);
  vec u(K);
  u(K-1)=1;
  vec pis(K);
  
  for(int iter=Nsamp_old; iter<(Nsamp_old+Nsamp); iter++){
    
    for(int k=0; k<(K-1); k++){
      u(k) = R::rbeta(1+n(k), alpha_K + sum(n.tail(K-k-1)));
    }
    pis=GetWeights(u, K);  
        
    n.fill(0);
    for(int i=0; i<N; i++){
      zeta(iter,i) = sample_zeta(K, pis, Pgrams.col(i), f, N_freq);
      n(zeta(iter,i)) += 1;
    } 
      
    field<uvec> inK(K);
    for(int k=0; k<K; k++){
      uvec inK2 = find(zeta.row(iter).t() == k);
      inK(k) = inK2;
    }
      
    //Update number of Bernstein polynomials (Lambda) in each cluster
    for(int k=0; k<K; k++){
      LambdaList = sample_Lambda(n(k), L, omegas, Pgrams.cols(inK(k)), P.slice(iter-1).col(k),
                                 Z.slice(iter-1).col(k), Lambda(iter-1,k), Lambda_prior, 
                                 Lambda_init, Lambda_max, MaxFreq, "poisson", f.col(k), 
                                 N_freq, false);
      Lambda(iter,k) = as<int>(LambdaList("Lambda_new"));
      f.col(k) = as<vec>(LambdaList("f_new"));
    }
      
    //Update BDP stick breaking weights, V, for each cluster
    for(int k=0; k<K; k++){
      VPList = sample_V(n(k), L, omegas, Pgrams.cols(inK(k)), V.col(k), 
                        P.slice(iter-1).col(k), Z.slice(iter-1).col(k), Lambda(iter,k),
                        epsilon, alpha_L, MaxFreq, f.col(k), N_freq, false);
      V.col(k) = as<vec>(VPList("V_new"));
      P.slice(iter).col(k) = as<vec>(VPList("P_new"));
      f.col(k) = as<vec>(VPList("f_new"));
    }
      
    //Update BDP atoms, Z, for each cluster
    for(int k=0; k<K; k++){
      ZList = sample_Z(n(k), L, omegas, Pgrams.cols(inK(k)), P.slice(iter).col(k),
                       Z.slice(iter-1).col(k), Lambda(iter,k), epsilon, MaxFreq, f.col(k), N_freq);
      Z.slice(iter).col(k) = as<vec>(ZList("Z_new"));
      f.col(k) = as<vec>(ZList("f_new"));
    }
    
    if((iter+1)%Sup == 0){cout << iter+1 << " posterior samples complete" << endl;}
  }
  
  return List::create(
    _["zeta"] = zeta,
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["Z"] = Z,
    _["omegas"] = omegas,
    _["Pgrams"] = Pgrams
  ); 
}    


//MCMC Sampling Algorithm for nested bernstein polynomials dirichlet process spectral estimator using split and merge sampler
// [[Rcpp::export]]
List SpectralNDBDP(const mat& X, const mat& Y, int Nsamp, double Zcoef_prop_sig=0.2, 
                   const vec& eps=0, int L=20, int K=35, double SamplingRate=1, 
                   double MaxFreq=0.5, int Sup=500, int Lambda_init=30, int Lambda_max=300, 
                   string Lambda_prior="flat", double alpha_K=1, double alpha_L=1,
                   bool warmStart=false, int Nwarm=1000){
  //X is a matrix where each column is a time series to be analyzed
  //Y is a matrix of covariates with a row for each subject
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //Zcoef_prop_sig is the Zcoefficient proposal density variance
  //K is the truncation level for the upper level DP
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
    //either "expon", "poisson", or "flat"
  //alpha_K and alpha_L are the concentration parameters for the DPs
  //warmStart is a boolean to determine if each cluster should be initialed at a fit from the BDP model
  //Nwarm is the number of BDP iterations to run for each cluster when warming up the clusters
  

  //determine Fourier frequencies to consider
  int N = X.n_cols;
  double T = X.n_rows;
  int N_freq_temp = floor((T-1)/2);
  int N_freq = 0;
  int N_pred = Y.n_cols;
  vec omegas(N_freq_temp);
  for(int n_freq=0; n_freq < N_freq_temp; n_freq++){
    omegas(n_freq) = (n_freq+1)/T*SamplingRate;
    if(omegas(n_freq) < MaxFreq){N_freq += 1;}
  }
  omegas=omegas.head(N_freq);

  //calculate periodograms
  mat Pgrams(N_freq, N);
  vec PgramTemp(T);
  for(int n=0; n < N; n++){
    PgramTemp = square(abs(fft(X.col(n))))/T;
    Pgrams.col(n) = PgramTemp.head(N_freq);
  }
  
  vec epsilon(L);
  if(eps(0)==0){
    for(int l=0; l<L; l++){
      epsilon(l) = (l+1)/(l+1+2*sqrt(T));
    }
  }else{
    epsilon = eps;
  }
  
  vec prop_sig_vec(N_pred);
  prop_sig_vec.fill(Zcoef_prop_sig);
  mat Zcoef_prop_sig_mat = diagmat(prop_sig_vec);

  //initiate structures to hold posterior samples for parameters
  mat Lambda(Nsamp,K);
  cube P(L,K,Nsamp);
  field<cube> Zcoef(Nsamp);
  mat zeta(Nsamp, N);
  
  //initialize parameter estimates randomly for each cluster
  mat V(L,K);
  cube Zcoeftemp(L, N_pred, K);
  
  if(warmStart==true){
    //initialize parameter estimates randomly for each cluster if warm start
    int NperG = floor(N/K);
    uvec perm_id(N);
    for(int i=0; i<N; i++){perm_id(i)=i;}
    perm_id = shuffle(perm_id);
    if(NperG >= 10){
      List GroupFit;
      for(int k=0; k<K; k++){
        cout<<"Warming up cluster "<<k+1<<endl;
        GroupFit = SpectralDBDP(X.cols(perm_id.subvec(k*NperG,k*NperG+NperG-1)), 
                                Y.rows(perm_id.subvec(k*NperG,k*NperG+NperG-1)),
                                Nwarm, Zcoef_prop_sig, epsilon, L, SamplingRate, MaxFreq, 
                                100000, Lambda_init, Lambda_max, Lambda_prior);
        Lambda(0,k) = as<vec>(GroupFit("Lambda"))(Nwarm-1);
        V.col(k) = as<mat>(GroupFit("V")).row(Nwarm-1).t();
        Zcoeftemp.slice(k) = as<cube>(GroupFit("Zcoef")).slice(Nwarm-1);
      }
    }else{
      List GroupFit;
      int InitGs = floor(N/10);
      for(int k=0; k<InitGs; k++){
        cout<<"Warming up cluster "<<k+1<<endl;
        GroupFit = SpectralDBDP(X.cols(perm_id.subvec(k*10,k*10+9)), 
                                Y.rows(perm_id.subvec(k*10,k*10+9)), 
                                Nwarm, Zcoef_prop_sig, epsilon, L, SamplingRate, MaxFreq, 
                                100000, Lambda_init, Lambda_max, Lambda_prior);
        Lambda(0,k) = as<vec>(GroupFit("Lambda"))(Nwarm-1);
        V.col(k) = as<mat>(GroupFit("V")).row(Nwarm-1).t();
        Zcoeftemp.slice(k) = as<cube>(GroupFit("Zcoef")).slice(Nwarm-1);
      }
      Lambda.row(0).tail(K-InitGs).fill(Lambda_init);
      V.submat(0,InitGs,L-2,K-1).fill(0.5);
      V.submat(L-1,InitGs,L-1,K-1).fill(1);
      for(int k=InitGs; k<K; k++){
        Zcoeftemp.slice(k).fill(0);
      }
    }
  }else{
    //initialize parameter estimates randomly for each cluster if no warm start
    Lambda.row(0).fill(Lambda_init);
    V.submat(0,0,L-2,K-1).fill(0.5);
    V.row(L-1).fill(1);
    for(int k=0; k<K; k++){
      Zcoeftemp.slice(k).fill(0);
    }
  }
  
  for(int k=0; k<K; k++){
    P.slice(0).col(k)=GetWeights(V.col(k), L);
  }
  Zcoef(0) = Zcoeftemp;
  cube Z(L, N, K);
  for(int k=0; k<K; k++){
    Z.slice(k) = ZcoefToZ(Y, Zcoef(0).slice(k), L, N);
  }
  
  mat Zcoef_mu(N_pred, K, fill::zeros);
  vec Zcoef_mu_pmu(N_pred);
  mat Zcoef_mu_pS(N_pred, N_pred);
  cube Zcoef_S(N_pred, N_pred, K);
  for(int k=0; k<K; k++){
    Zcoef_S.slice(k).eye();
  }
  mat SInv(N_pred, N_pred);
  mat S0Inv(N_pred, N_pred);
  S0Inv.eye();
  mat IWvar(N_pred, N_pred);
  
  vec u(K);
  vec pis(K);
  u.head(K-1).fill(0.5);
  u(K-1) = 1;
  pis = GetWeights(u, K);
  
  vec n(K, fill::zeros);
  mat f_subj(N_freq, K);
  mat f_cur(N_freq, N);
  for(int i=0; i < N; i++){
    for(int k=0; k < K; k++){
      f_subj.col(k) = Get_f(omegas, P.slice(0).col(k), Z.slice(k).col(i), Lambda(0,k), 
                            L, MaxFreq, N_freq);
    }
    zeta(0,i) = sample_zeta(K, pis, Pgrams.col(i), f_subj, N_freq);
    n(zeta(0,i)) += 1;
    f_cur.col(i) = f_subj.col(zeta(0,i));
  }
  field<uvec> inK(K);
  for(int k=0; k < K; k++){
    inK(k) = find(zeta.row(0).t() == k);
  }
  
  List VPList(3);
  List LambdaList(2);
  List ZList(3);
  
  for(int iter=1; iter<Nsamp; iter++){
    
    for(int k=0; k<(K-1); k++){
      u(k) = R::rbeta(1+n(k), alpha_K + sum(n.tail(K-k-1)));
    }
    pis = GetWeights(u, K);
    n.fill(0);
    for(int i=0; i < N; i++){
      for(int k=0; k < K; k++){
        f_subj.col(k) = Get_f(omegas, P.slice(iter-1).col(k), Z.slice(k).col(i), 
                              Lambda(iter-1,k), L, MaxFreq, N_freq);
      }
      zeta(iter,i) = sample_zeta(K, pis, Pgrams.col(i), f_subj, N_freq);
      n(zeta(iter,i)) += 1;
      f_cur.col(i) = f_subj.col(zeta(iter,i));
    }
      
    field<uvec> inK(K);
    for(int k=0; k<K; k++){
      uvec inK2 = find(zeta.row(iter).t() == k);
      inK(k) = inK2;
    }
    
    //Update number of Bernstein polynomials (Lambda) in each cluster
    for(int k=0; k<K; k++){
      LambdaList = sample_Lambda(n(k), L, omegas, Pgrams.cols(inK(k)), P.slice(iter-1).col(k),
                 Z.slice(k).cols(inK(k)), Lambda(iter-1,k), Lambda_prior, Lambda_init,
                 Lambda_max, MaxFreq, "poisson", f_cur.cols(inK(k)), N_freq, true);
      Lambda(iter,k) = as<int>(LambdaList("Lambda_new"));
      f_cur.cols(inK(k)) = as<mat>(LambdaList("f_new"));
    }
    
    //Update BDP stick breaking weights, V, for each cluster
    for(int k=0; k<K; k++){
      VPList = sample_V(n(k), L, omegas, Pgrams.cols(inK(k)), V.col(k),
                       P.slice(iter-1).col(k), Z.slice(k).cols(inK(k)), Lambda(iter,k),
                       epsilon, alpha_L, MaxFreq, f_cur.cols(inK(k)), N_freq, true);
      V.col(k) = as<vec>(VPList("V_new"));
      P.slice(iter).col(k) = as<vec>(VPList("P_new"));
      f_cur.cols(inK(k)) = as<mat>(VPList("f_new"));
    }
    
    //Update BDP atom coefficients, Zcoef, for each cluster
    for(int k=0; k<K; k++){
      SInv = inv(Zcoef_S.slice(k));
      ZList = sample_Zcoef(n(k), L, omegas, Pgrams.cols(inK(k)), Y.rows(inK(k)), 
                           P.slice(iter).col(k), Zcoef(iter-1).slice(k),Z.slice(k).cols(inK(k)), 
                           Lambda(iter,k), Zcoef_mu.col(k), SInv, Zcoef_prop_sig_mat, 
                           MaxFreq, f_cur.cols(inK(k)), N_freq);
      Zcoeftemp.slice(k) = as<mat>(ZList("Zcoef_new"));
      f_cur.cols(inK(k)) = as<mat>(ZList("f_new"));
      Z.slice(k) = ZcoefToZ(Y, Zcoeftemp.slice(k), L, N);

      Zcoef_mu_pS = inv(S0Inv + L*SInv);
      Zcoef_mu_pmu = Zcoef_mu_pS * (SInv*sum(Zcoeftemp.slice(k),0).t());
      Zcoef_mu.col(k) = mvrnormArma(1, Zcoef_mu_pmu, Zcoef_mu_pS).t();
      IWvar.eye();
      for(int l=0; l < L; l++){
        IWvar = IWvar + ( (Zcoef_mu.col(k) - Zcoeftemp.slice(k).row(l).t()) *
          (Zcoef_mu.col(k).t() - Zcoeftemp.slice(k).row(l)) );
      }
      Zcoef_S.slice(k) = riwishart(L+N_pred, IWvar);
    }
    Zcoef(iter) = Zcoeftemp;
    
    if((iter+1)%Sup == 0){cout << iter+1 << " posterior samples complete" << endl;}
  }

  return List::create(
    _["zeta"] = zeta,
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["Zcoef"] = Zcoef,
    _["omegas"] = omegas,
    _["Pgrams"] = Pgrams
  );
}


//MCMC Sampling Algorithm for nested bernstein polynomials dirichlet process spectral estimator using split and merge sampler
// [[Rcpp::export]]
List SpectralNDBDP2(const List& Fit_old, const mat& Y, int Nsamp, int T, double Zcoef_prop_sig=0.2, 
                    const vec& eps=0, double SamplingRate=1, double MaxFreq=0.5, int Sup=500, 
                    int Lambda_init=30, int Lambda_max=300, string Lambda_prior="flat", 
                    double alpha_K=1, double alpha_L=1){
  //Y is a matrix of covariates with a row for each subject
  //Nsamp is the number of posterior samples
  //L is the truncation level for the BDP
  //epsilon is a vector of uniform halfwidths for V and Z proposals
  //Zcoef_prop_sig is the Zcoefficient proposal density variance
  //K is the truncation level for the upper level DP
  //SamplingRate is the sampling rate of the time series
  //MaxFreq is the maximum frequnecy to consider for the analysis, at most 1/2 the sampling rate
  //Sup is how often status updates are printed to the console, every Sup iterations
  //Lambda_init is the initial value for Lambda and the value for all samples if Lambda is set to constant
  //Lambda_max is the maximum number of Bernstein polynomial components to consider
  //Lambda_prior is the prior on the number of Bernstein polynomial components, Lambda.
    //either "expon", "poisson", or "flat"
  //alpha_K and alpha_L are the concentration parameters for the DPs
  //warmStart is a boolean to determine if each cluster should be initialed at a fit from the BDP model
  //Nwarm is the number of BDP iterations to run for each cluster when warming up the clusters

  vec omegas = as<vec>(Fit_old("omegas"));
  mat Pgrams = as<mat>(Fit_old("Pgrams"));
  mat zeta_old = as<mat>(Fit_old("zeta"));
  mat Lambda_old = as<mat>(Fit_old("Lambda"));
  cube P_old = as<cube>(Fit_old("P"));
  field<cube> Zcoef_old = as<field<cube>>(Fit_old("Zcoef"));
  int Nsamp_old = zeta_old.n_rows;
  int N_freq = Pgrams.n_rows;
  int N = Pgrams.n_cols;
  int N_pred = Y.n_cols;
  int K = P_old.n_cols;
  int L = P_old.n_rows;
  vec epsilon(L);
  if(eps(0)==0){
    for(int l=0; l<L; l++){
      epsilon(l) = (l+1)/(l+1+2*sqrt(T));
    }
  }else{
    epsilon = eps;
  }
  
  mat Lambda(Nsamp+Nsamp_old,K);
  cube P(L,K,Nsamp+Nsamp_old);
  field<cube> Zcoef(Nsamp+Nsamp_old);
  mat zeta(Nsamp+Nsamp_old, N);
  
  zeta.rows(0,Nsamp_old-1) = zeta_old;
  Lambda.rows(0,Nsamp_old-1) = Lambda_old;
  P.slices(0,Nsamp_old-1) = P_old;
  for(int i=0; i<Nsamp_old; i++){
    Zcoef(i) = Zcoef_old(i);
  }
  
  cube Zcoeftemp = Zcoef(Nsamp_old-1);  
  mat V(L, K);
  for(int k=0; k<K; k++){
    V.col(k) = GetWeightsInv(P.slice(Nsamp_old-1).col(k), L);
  }
  
  vec prop_sig_vec(N_pred);
  prop_sig_vec.fill(Zcoef_prop_sig);
  mat Zcoef_prop_sig_mat = diagmat(prop_sig_vec);
  
  cube Z(L, N, K);
  for(int k=0; k<K; k++){
    Z.slice(k) = ZcoefToZ(Y, Zcoeftemp.slice(k), L, N);
  }
  
  mat Zcoef_mu(N_pred, K, fill::zeros);
  vec Zcoef_mu_pmu(N_pred);
  mat Zcoef_mu_pS(N_pred, N_pred);
  cube Zcoef_S(N_pred, N_pred, K);
  mat SInv(N_pred, N_pred);
  mat S0Inv(N_pred, N_pred);
  S0Inv.eye();
  mat IWvar(N_pred, N_pred);
  for(int k=0; k<K; k++){
    Zcoef_S.slice(k).eye();
    Zcoef_mu_pS = inv(S0Inv + L*SInv);
    Zcoef_mu_pmu = Zcoef_mu_pS * (SInv*sum(Zcoeftemp.slice(k),0).t());
    Zcoef_mu.col(k) = mvrnormArma(1, Zcoef_mu_pmu, Zcoef_mu_pS).t();
    IWvar.eye();
    for(int l=0; l < L; l++){
      IWvar = IWvar + ( (Zcoef_mu.col(k) - Zcoeftemp.slice(k).row(l).t()) *
        (Zcoef_mu.col(k).t() - Zcoeftemp.slice(k).row(l)) );
    }
    Zcoef_S.slice(k) = riwishart(L+N_pred, IWvar);
  }
  
  vec n(K, fill::zeros);
  mat f_subj(N_freq, K);
  mat f_cur(N_freq, N);
  for(int i=0; i < N; i++){
    for(int k=0; k < K; k++){
      f_subj.col(k) = Get_f(omegas, P.slice(Nsamp_old-1).col(k), Z.slice(k).col(i), 
                            Lambda(Nsamp_old-1,k), L, MaxFreq, N_freq);
    }
    n(zeta(Nsamp_old-1,i)) += 1;
    f_cur.col(i) = f_subj.col(zeta(Nsamp_old-1,i));
  }
  field<uvec> inK(K);
  for(int k=0; k < K; k++){
    inK(k) = find(zeta.row(Nsamp_old-1).t() == k);
  }
  
  List VPList(3);
  List LambdaList(2);
  List ZList(3);
  vec u(K);
  u(K-1) = 1;
  vec pis(K);
  
  for(int iter=Nsamp_old; iter<(Nsamp_old+Nsamp); iter++){
    
    for(int k=0; k<(K-1); k++){
      u(k) = R::rbeta(1+n(k), alpha_K + sum(n.tail(K-k-1)));
    }
    pis = GetWeights(u, K);
    n.fill(0);
    for(int i=0; i < N; i++){
      for(int k=0; k < K; k++){
        f_subj.col(k) = Get_f(omegas, P.slice(iter-1).col(k), Z.slice(k).col(i), 
                              Lambda(iter-1,k), L, MaxFreq, N_freq);
      }
      zeta(iter,i) = sample_zeta(K, pis, Pgrams.col(i), f_subj, N_freq);
      n(zeta(iter,i)) += 1;
      f_cur.col(i) = f_subj.col(zeta(iter,i));
    }
      
    field<uvec> inK(K);
    for(int k=0; k<K; k++){
      uvec inK2 = find(zeta.row(iter).t() == k);
      inK(k) = inK2;
    }
      
    //Update number of Bernstein polynomials (Lambda) in each cluster
    for(int k=0; k<K; k++){
      LambdaList = sample_Lambda(n(k), L, omegas, Pgrams.cols(inK(k)), P.slice(iter-1).col(k),
                                 Z.slice(k).cols(inK(k)), Lambda(iter-1,k), Lambda_prior, Lambda_init,
                                 Lambda_max, MaxFreq, "poisson", f_cur.cols(inK(k)), N_freq, true);
      Lambda(iter,k) = as<int>(LambdaList("Lambda_new"));
      f_cur.cols(inK(k)) = as<mat>(LambdaList("f_new"));
    }
      
    //Update BDP stick breaking weights, V, for each cluster
    for(int k=0; k<K; k++){
      VPList = sample_V(n(k), L, omegas, Pgrams.cols(inK(k)), V.col(k),
                        P.slice(iter-1).col(k), Z.slice(k).cols(inK(k)), Lambda(iter,k),
                        epsilon, alpha_L, MaxFreq, f_cur.cols(inK(k)), N_freq, true);
      V.col(k) = as<vec>(VPList("V_new"));
      P.slice(iter).col(k) = as<vec>(VPList("P_new"));
      f_cur.cols(inK(k)) = as<mat>(VPList("f_new"));
    }
    
    //Update BDP atom coefficients, Zcoef, for each cluster
    for(int k=0; k<K; k++){
      SInv = inv(Zcoef_S.slice(k));
      ZList = sample_Zcoef(n(k), L, omegas, Pgrams.cols(inK(k)), Y.rows(inK(k)), 
                           P.slice(iter).col(k), Zcoef(iter-1).slice(k),Z.slice(k).cols(inK(k)), 
                           Lambda(iter,k), Zcoef_mu.col(k), SInv, Zcoef_prop_sig_mat, 
                           MaxFreq, f_cur.cols(inK(k)), N_freq);
      Zcoeftemp.slice(k) = as<mat>(ZList("Zcoef_new"));
      f_cur.cols(inK(k)) = as<mat>(ZList("f_new"));
      Z.slice(k) = ZcoefToZ(Y, Zcoeftemp.slice(k), L, N);
      
      Zcoef_mu_pS = inv(S0Inv + L*SInv);
      Zcoef_mu_pmu = Zcoef_mu_pS * (SInv*sum(Zcoeftemp.slice(k),0).t());
      Zcoef_mu.col(k) = mvrnormArma(1, Zcoef_mu_pmu, Zcoef_mu_pS).t();
      IWvar.eye();
      for(int l=0; l < L; l++){
        IWvar = IWvar + ( (Zcoef_mu.col(k) - Zcoeftemp.slice(k).row(l).t()) *
          (Zcoef_mu.col(k).t() - Zcoeftemp.slice(k).row(l)) );
      }
      Zcoef_S.slice(k) = riwishart(L+N_pred, IWvar);
    }
    Zcoef(iter) = Zcoeftemp;
   
    if((iter+1)%Sup == 0){cout << iter+1 << " posterior samples complete" << endl;}
  }
  
  return List::create(
    _["zeta"] = zeta,
    _["Lambda"] = Lambda,
    _["P"] = P,
    _["Zcoef"] = Zcoef,
    _["omegas"] = omegas,
    _["Pgrams"] = Pgrams
  );
}




//function to calculate posteior curve samples from posteior samples for BDP model
// [[Rcpp::export]]
mat GetCurvesBDP(const vec& omegas, const vec& Lambda_samps, const mat& P_samps, 
                 const mat& Z_samps, double MaxFreq){
  //omegas is a vector of fourier frequencies
  //Lambda_samps is a vector of samples of the number of Bernstein polynomial components
  //P_samps is a matrix where each row is an L length sample of BDP weights
  //Z_samps is a matrix where each row is an L length sample of BDP atoms
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  int N_freq = omegas.n_elem;
  int L = P_samps.n_cols;
  int Nsamp = P_samps.n_rows;
  mat CurveEsts(Nsamp,N_freq,fill::zeros);
  for(int s=0; s<Nsamp; s++){
    CurveEsts.row(s) = Get_f(omegas, P_samps.row(s).t(), Z_samps.row(s).t(),
                  Lambda_samps(s), L, MaxFreq, N_freq).t();
  }
  return CurveEsts;
}


//function to calculate posteior curve samples from posteior variables for a single subject
// [[Rcpp::export]]
mat GetCurvesDBDP(const vec& omegas, const vec& Lambda_samps, const mat& P_samps,
                  const cube& Zcoef_samps, double MaxFreq, const vec& Y){
  //omegas is a vector of fourier frequencies
  //Lambda_samps is a vector of the # of bernstein polynomial components
  //P_samps is a matrix where each row contains an iteration of DBDP weights
  //Zcoef_samps is a cube, each slice containing an iteration of DBDP atom coefficients, coefficients in the columns
  //MaxFreq is the maximum frequency to consider for the spectral density
  //Y is a vector of covariates for the given time series
  int Nsamp=Lambda_samps.n_elem;
  int L=P_samps.n_cols;
  vec Z_now(L);
  int N_freq=omegas.n_elem;
  mat CurveEsts(Nsamp,N_freq,fill::zeros);
  for(int s=0; s<Nsamp; s++){
    Z_now = ZcoefToZ(Y.t(), Zcoef_samps.slice(s), L, 1);
    CurveEsts.row(s) = Get_f(omegas, P_samps.row(s).t(), Z_now, Lambda_samps(s), L, MaxFreq, N_freq).t();
  }
  return CurveEsts;
}


//function to calculate posteior curve samples for a single subject for NBDP
// [[Rcpp::export]]
mat GetCurvesNBDP(const vec& omegas, const vec& zeta_samps, const mat& Lambda_samps, 
                  const cube& P_samps, const cube& Z_samps, double MaxFreq){
  //omegas is a vector of fourier frequencies
  //zeta_samps is a vector of cluster assignments for a given subject
  //Lambda_samps is a matrix of samples of the number of Bernstein polynomial components with a column for each cluster
  //P_samps is a cube of BDP weights with a column per cluster and a slice per iteration
  //Z_samps is a cube of BDP weights with a column per cluster and a slice per iteration
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  int Nsamp=zeta_samps.n_elem;
  int L=P_samps.n_rows;
  int k_now;
  int Lambda_now;
  vec P_now(L);
  vec Z_now(L);
  int N_freq = omegas.n_elem;
  mat CurveEsts(Nsamp,N_freq,fill::zeros);
  for(int s=0; s<Nsamp; s++){
    k_now=zeta_samps(s);
    Lambda_now=Lambda_samps(s,k_now);
    P_now=P_samps.slice(s).col(k_now);
    Z_now=Z_samps.slice(s).col(k_now);
    CurveEsts.row(s) = Get_f(omegas, P_now, Z_now, Lambda_now, L, MaxFreq, N_freq).t();
  }
  return CurveEsts;
}


//function to calculate posteior curve samples for a single subject for NDBDP
// [[Rcpp::export]]
mat GetCurvesNDBDP(const vec& omegas, const vec& zeta_samps, const mat& Lambda_samps, 
                   const cube& P_samps, const field<cube>& Zcoef_samps, double MaxFreq, 
                   const vec& Y){
  //omegas is a vector of fourier frequencies
  //zeta_samps is a vector of cluster assignments for a given subject
  //Lambda_samps is a matrix of samples of the number of Bernstein polynomial components with a column for each cluster
  //P_samps is a cube of BDP weights with a column per cluster and a slice per iteration
  //Zcoef_samps is a field of cubes of BDP atom coefficients
  //MaxFreq is the maximum frequency to consider for the spectral density
  
  int Nsamp=zeta_samps.n_elem;
  int L=P_samps.n_rows;
  int k_now;
  int Lambda_now;
  vec P_now(L);
  vec Z_now(L);
  int N_freq = omegas.n_elem;
  mat CurveEsts(Nsamp,N_freq,fill::zeros);
  for(int s=0; s<Nsamp; s++){
    k_now=zeta_samps(s);
    Lambda_now=Lambda_samps(s,k_now);
    P_now=P_samps.slice(s).col(k_now);
    Z_now = ZcoefToZ(Y.t(), Zcoef_samps(s).slice(k_now), L, 1);
    CurveEsts.row(s) = Get_f(omegas, P_now, Z_now, Lambda_now, L, MaxFreq, N_freq).t();
  }
  return CurveEsts;
}


//function to calculate heritability of spectrum for NBDP
// [[Rcpp::export]]
List GetHeritNBDP(const List& Fit, double MaxFreq, int burnin, const uvec& MZ_t1_ind,
                  const uvec& MZ_t2_ind, const uvec& DZ_t1_ind, const uvec& DZ_t2_ind,
                  bool log_spec){
  //Fit is a List object containing output from the SpectralNBDP() function
  //MaxFreq is the maximum frequency to consider for the spectral density
  //burnin is the number of samples to be discarded as burnin
  //MZ_t1_ind and MZ_t2_ind contain indices for MZ twin pairs twins 1 and 2, respectively
  //DZ_t1_ind and DZ_t2_ind contain indices for DZ twin pairs twins 1 and 2, respectively
  //log_spec is a boolean to determine if heritability should be calculated for f or log(f)

  vec omegas = as<vec>(Fit("omegas"));
  mat zeta_samps = as<mat>(Fit("zeta"));
  mat Lambda_samps = as<mat>(Fit("Lambda"));
  cube P_samps = as<cube>(Fit("P"));
  cube Z_samps = as<cube>(Fit("Z"));

  int N = zeta_samps.n_cols;
  int Nsamp = zeta_samps.n_rows;
  int L = P_samps.n_rows;
  int K = P_samps.n_cols;
  int N_freq = omegas.n_elem;
  int N_MZ = MZ_t1_ind.n_elem;
  int N_DZ = DZ_t1_ind.n_elem;

  vec P_now(L);
  vec Z_now(L);
  int Lambda_now;

  mat MZ_t1(N_MZ, N_freq);
  mat MZ_t2(N_MZ, N_freq);
  mat DZ_t1(N_DZ, N_freq);
  mat DZ_t2(N_DZ, N_freq);
  mat Vars(Nsamp-burnin, N_freq);
  mat HeritCurve(Nsamp-burnin, N_freq);

  mat SubjCurves(N, N_freq);
  mat Curves(K,N_freq);
  for(int s=0; s<(Nsamp-burnin); s++){
    for(int k=0; k<K; k++){
      Lambda_now = Lambda_samps(s+burnin,k);
      P_now = P_samps.slice(s+burnin).col(k);
      Z_now = Z_samps.slice(s+burnin).col(k);
      if(log_spec == true){
        Curves.row(k) = log(Get_f(omegas, P_now, Z_now, Lambda_now, L, MaxFreq, N_freq)).t();
      }else{
        Curves.row(k) = Get_f(omegas, P_now, Z_now, Lambda_now, L, MaxFreq, N_freq).t();  
      }
    }
    for(int i=0; i<N; i++){
      SubjCurves.row(i) = Curves.row(zeta_samps(s+burnin, i));
    }
    MZ_t1 = SubjCurves.rows(MZ_t1_ind);
    MZ_t2 = SubjCurves.rows(MZ_t2_ind);
    DZ_t1 = SubjCurves.rows(DZ_t1_ind);
    DZ_t2 = SubjCurves.rows(DZ_t2_ind);
    for(int n = 0; n < N_freq; n++){
      HeritCurve(s,n) = 2*(as_scalar(cor(MZ_t1.col(n), MZ_t2.col(n))) -
                           as_scalar(cor(DZ_t1.col(n), DZ_t2.col(n))));
      Vars(s,n) = var(SubjCurves.col(n));
    }
  }
  return List::create(
    _["HeritCurve"] = HeritCurve,
    _["CurveVariance"] = Vars
  );
}


//function to calculate heritability of spectrum for NBDP
// [[Rcpp::export]]
List GetHeritNDBDP(const List& Fit, const mat& Y, double MaxFreq, int burnin,
                   const uvec& MZ_t1_ind, const uvec& MZ_t2_ind,
                   const uvec& DZ_t1_ind, const uvec& DZ_t2_ind,
                   bool log_spec){
  //Fit is a List object containing output from the SpectralNDBDP() function
  //MaxFreq is the maximum frequency to consider for the spectral density
  //burnin is the number of samples to be discarded as burnin
  //MZ_t1_ind and MZ_t2_ind contain indices for MZ twin pairs twins 1 and 2, respectively
  //DZ_t1_ind and DZ_t2_ind contain indices for DZ twin pairs twins 1 and 2, respectively
  //log_spec is a boolean to determine if heritability should be calculated for f or log(f)

  vec omegas = as<vec>(Fit("omegas"));
  mat zeta_samps = as<mat>(Fit("zeta"));
  mat Lambda_samps = as<mat>(Fit("Lambda"));
  cube P_samps = as<cube>(Fit("P"));
  field<cube> Zcoef_samps = as<field<cube>>(Fit("Zcoef"));
  
  int Nsamp = zeta_samps.n_rows;
  int N = zeta_samps.n_cols;
  int L = P_samps.n_rows;
  int N_freq = omegas.n_elem;
  int N_MZ = MZ_t1_ind.n_elem;
  int N_DZ = DZ_t1_ind.n_elem;

  int Lambda_now;
  vec P_now(L);
  vec Z_now(L);
  int K_now;

  mat MZ_t1(N_MZ, N_freq);
  mat MZ_t2(N_MZ, N_freq);
  mat DZ_t1(N_DZ, N_freq);
  mat DZ_t2(N_DZ, N_freq);
  mat Vars(Nsamp-burnin, N_freq);
  mat HeritCurve(Nsamp-burnin, N_freq);

  mat Curves(N,N_freq);
  for(int s=0; s<(Nsamp-burnin); s++){
    for(int i=0; i<N; i++){
      K_now = zeta_samps(s+burnin,i);
      Lambda_now = Lambda_samps(s+burnin,K_now);
      P_now = P_samps.slice(s+burnin).col(K_now);
      Z_now = ZcoefToZ(Y.row(i), Zcoef_samps(s+burnin).slice(K_now), L, 1);
      if(log_spec == true){
        Curves.row(i) = log(Get_f(omegas, P_now, Z_now, Lambda_now, L, MaxFreq, N_freq)).t();
      }else{
        Curves.row(i) = Get_f(omegas, P_now, Z_now, Lambda_now, L, MaxFreq, N_freq).t();  
      }
    }
    
    MZ_t1 = Curves.rows(MZ_t1_ind);
    MZ_t2 = Curves.rows(MZ_t2_ind);
    DZ_t1 = Curves.rows(DZ_t1_ind);
    DZ_t2 = Curves.rows(DZ_t2_ind);
    for(int n = 0; n < N_freq; n++){
      HeritCurve(s,n) = 2*(as_scalar(cor(MZ_t1.col(n), MZ_t2.col(n))) -
                           as_scalar(cor(DZ_t1.col(n), DZ_t2.col(n))));
      Vars(s,n) = var(Curves.col(n));
    }
  }
  return List::create(
    _["HeritCurve"] = HeritCurve,
    _["CurveVariance"] = Vars
  );
}




