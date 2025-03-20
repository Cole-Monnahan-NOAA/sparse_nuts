#include <TMB.hpp>
#include <Eigen/Dense>  // Include Eigen for matrix operations

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(x);
  DATA_VECTOR(k);

  // fixed effects -- bounds added in R
  PARAMETER(logrho);
  PARAMETER(logalpha);
  PARAMETER_VECTOR(f_tilde);

  Type alpha = exp(logalpha);
  Type rho = exp(logrho);
  matrix<Type> cov(11,11);

  for(int i=0; i<11; i++){
    for(int j=0; j<11; j++){
      cov(i,j)= alpha*alpha*exp(-pow(x(i)-x(j),2)/(2*rho*rho));
    }
  }
  for(int i=0; i<11; i++) cov(i,i)+=1e-10;

  // Compute the Cholesky decomposition using Eigen's LLT decomposition
  Eigen::LLT<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> llt(cov);
  matrix<Type> L_chol = llt.matrixL();
  vector<Type> f=L_chol * f_tilde;
  REPORT(L_chol);

  Type lp =
    // priors
    dgamma(rho, Type(25), Type(0.25), true)+
    dnorm(alpha, Type(0), Type(2),true)+
    sum(dnorm(f_tilde,Type(0),Type(1),true)) + // hyperdistribution
    logrho + logalpha + // Jacobians
    sum(dpois(k, exp(f), true)); // likelihood
  REPORT(f);
  REPORT(lp);
  return(-lp);

}
