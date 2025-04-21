#include <TMB.hpp>
#include <Eigen/Dense>  // Include Eigen for matrix operations

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  //DATA_INTEGER(d);
  DATA_VECTOR(y);
  DATA_MATRIX(x);
  DATA_SCALAR(scale_icept);
  DATA_SCALAR(scale_global);
  DATA_SCALAR(nu_global);
  DATA_SCALAR(nu_local);
  DATA_SCALAR(slab_scale);
  DATA_SCALAR(slab_df);
  //DATA_INTEGER(p);
  
 // Transformed data
 //  Type d0=5;
  Type slab_scale2=slab_scale*slab_scale;
  Type half_slab_df=0.5*slab_df;
  vector<Type> mu(n);
  mu.setZero();
  matrix<Type> x2=x.array()*x.array();

  PARAMETER_VECTOR(z);
  PARAMETER(logtau);
  PARAMETER_VECTOR(loglambda);
  PARAMETER(logcaux);
  PARAMETER(logxi);
  int d=loglambda.size();
  
  // Transformed parameters
  Type tau = exp(logtau);
  vector<Type> lambda = exp(loglambda);
  Type caux = exp(logcaux);
  Type xi = exp(logxi);
  Type eta_one = scale_global*tau;
  Type m_squared = slab_scale2*caux;
  vector<Type> kappa_squared = m_squared*lambda*lambda/(m_squared+eta_one*eta_one*lambda*lambda);
  Type eta_two = eta_one*eta_one/m_squared*xi;
  matrix<Type> dkappa2(d,d);
  dkappa2.setZero();
  for(int i=0; i<d; i++) dkappa2(i,i)=kappa_squared(i);
  
  matrix<Type> K1 = (x*dkappa2) * x.transpose();
  matrix<Type> K2 = (x2*dkappa2) * x2.transpose();
 matrix<Type> K1sq(n,n); // element wise (K1+1)(K1+1)
 K1sq.setOnes(); 
 K1sq+=K1;
 K1sq=K1sq.array()*K1sq.array();
 matrix<Type> I(n,n);
 I.setOnes(); // need this for scalar +/- matrix calcs below?
  matrix<Type> K(n,n); 
  K.setOnes();
  K = 0.5*eta_two*eta_two*K1sq - 
  0.5*eta_two*eta_two*K2 +
  (eta_one*eta_one-eta_two*eta_two)*K1 +
  scale_icept*scale_icept*I -
  (0.5*eta_two*eta_two)*I;
  for(int i=0; i<n; i++) K(i,i)+=1e-5;
   // Compute the Cholesky decomposition using Eigen's LLT decomposition
  Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> llt(K);
  matrix<Type> L_K = llt.matrixL();
  vector<Type> func_vals=mu+ L_K*z;
  REPORT(L_K);
   vector<Type> pred=1/(1+exp(-func_vals));
  Type lp = 
  // Jacobians
  logtau +sum(loglambda)+logcaux + logxi +
    // priors
	dt(lambda, nu_local, true).sum()+
	dt(tau, nu_global, true) +
	 (dgamma(1/caux, half_slab_df, 1/half_slab_df, true)-2*log(caux))+
	 (dgamma(1/xi, half_slab_df, 1/half_slab_df, true)-2*log(xi))+
	 dnorm(z,0,1,true).sum(); 
	 
	 for(int i=0; i<n; i++) lp+=dbinom(y(i), Type(1), pred(i),true);

  REPORT(K1);
  REPORT(K2);
  REPORT(K);
  REPORT(func_vals);
  REPORT(lp);
  return(-lp);

}
