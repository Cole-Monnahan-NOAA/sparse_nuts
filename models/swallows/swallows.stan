// This model was adapted from section 14.5 of Korner-Nievergelt et al
// 2015. It is a Cormack-Jolly-Seber model.

// Last updated on 3/29/2016
data {
  int<lower=0> I;                       // number of individuals
  int<lower=2> K;                       // capture events
  int<lower=0> nfam;                    // number of families
  int<lower=0,upper=1> CH[I,K];         // CH[i,k]: individual i captured at k
  vector[I] carez;                     // covariable, duration of parental care, z-trans.
  int<lower=1,upper=4> year[I];         // index of year
  vector[K] agec;                      // age of fledling, centered
  int<lower=0, upper=nfam> family[I];   // index of group variable
  int<lower=0> last[I];                 // vector of last time the bird was seen
  int ones[I,K];			// matrix of 1's, for use in "ones trick"
  int ones2[I];				// vector of 1's, for use in "ones trick"
}


parameters {
  real<lower=0> sigmayearphi;  // between-year standard deviation in logit(phi)
  real<lower=0> sigmaphi;      // between family standard deviation in logit(phi)
  real<lower=0> sigmap;        // between family standard deviation in logit(p)
  real a[K-1];                   // intercept of phi
  real a1;                     // coef of phi
  real b0[4];                  // intercepts per year for p
  real b1[4];                  // slope for age per year for p
  real fameffphi_raw[nfam];        // family effects for phi
  real fameffp_raw[nfam];          // family effects for p
  real yeareffphi_raw[4];          // year effect on phi
}

transformed parameters {
  real<lower=0,upper=1>p[I,K];       // capture probability
  real<lower=0,upper=1>phi[I,K-1];   // survival probability
  real<lower=0,upper=1>chi[I,K+1];   // probability never seen again
  {
  int k;
  for(i in 1:I){ // loop over each individual
    // calculate phi as a function of fixed and random effects
    for(t in 1:(K-1)) {
      phi[i,t] = inv_logit(a[t]+ a1*carez[i]+ sigmayearphi*yeareffphi_raw[year[i]]+
			    sigmaphi*fameffphi_raw[family[i]]);
    }
    // calculate p as a function of fixed and random effectsa
    p[i,1] = 1;  // first occasion is marking occasion
    for(t in 2:K){
      p[i,t] = inv_logit(b0[year[i]] + b1[year[i]]*agec[t]+
			  sigmap*fameffp_raw[family[i]]);
    }
    // probabilitiy of never being seen after last observation. ind here is
    // a reverse index so this loop goes from K:2, recursively calculating
    // backward.
    chi[i,K+1] = 1.0;
    k = K;
    while (k > 1) {
      chi[i,k] = (1 - phi[i,k-1]) + phi[i,k-1] * (1 - p[i,k]) * chi[i,k+1];
      k = k - 1;
    }
    chi[i,1] = (1 - p[i,1]) * chi[i,2];
  }
  }
}


model {
  // priors, vectorized
  b0~normal(0,5);
  b1~normal(0,5);
  a~normal(0,1.5);
  a1~normal(0,5);
  sigmaphi~cauchy(0,1);
  sigmayearphi~normal(0,3);
  sigmap~cauchy(0,1);
  // random effects, vectorized
  fameffphi_raw~normal(0, 1);
  fameffp_raw~normal(0,1);
  yeareffphi_raw~normal(0, 1);
  // likelihood
  for (i in 1:I) {
    // probability of survival, known alive since k<last
    for (t in 2:last[i]) {
      ones[i,t]~bernoulli(phi[i, t-1]);
    }
    // probability of observation given known alive
    for(t in 1:last[i]){
      CH[i,t]~bernoulli(p[i,t]);
    }
    // probability of no observations after time period last
    ones2[i]~bernoulli(chi[i,last[i]+1]);
  }
}


