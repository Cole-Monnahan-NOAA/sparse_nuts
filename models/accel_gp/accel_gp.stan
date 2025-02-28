// generated with brms 2.22.0
functions {
  /* Spectral density function of a Gaussian process
   * with squared exponential covariance kernel
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns:
   *   numeric vector of length NB of the SPD evaluated at 'x'
   */
  vector spd_gp_exp_quad(data array[] vector x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    real constant = square(sdgp) * sqrt(2 * pi())^D;
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real neg_half_lscale2 = -0.5 * square(lscale[1]);
      constant = constant * lscale[1]^D;
      for (m in 1:NB) {
        out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
      }
    } else {
      // multi-dimensional non-isotropic GP
      vector[Dls] neg_half_lscale2 = -0.5 * square(lscale);
      constant = constant * prod(lscale);
      for (m in 1:NB) {
        out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m])));
      }
    }
    return out;
  }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  // data related to GPs
  int<lower=1> Kgp_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_1;  // GP dimension
  // number of basis functions of an approximate GP
  int<lower=1> NBgp_1;
  // number of latent GP groups
  int<lower=1> Nsubgp_1;
  // indices of latent GP groups per observation
  array[N] int<lower=1> Jgp_1;
  // approximate GP basis matrices
  matrix[Nsubgp_1, NBgp_1] Xgp_1;
  // approximate GP eigenvalues
  array[NBgp_1] vector[Dgp_1] slambda_1;
  // data related to GPs
  int<lower=1> Kgp_sigma_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_sigma_1;  // GP dimension
  // number of basis functions of an approximate GP
  int<lower=1> NBgp_sigma_1;
  // number of latent GP groups
  int<lower=1> Nsubgp_sigma_1;
  // indices of latent GP groups per observation
  array[N] int<lower=1> Jgp_sigma_1;
  // approximate GP basis matrices
  matrix[Nsubgp_sigma_1, NBgp_sigma_1] Xgp_sigma_1;
  // approximate GP eigenvalues
  array[NBgp_sigma_1] vector[Dgp_sigma_1] slambda_sigma_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[Kgp_1] sdgp_1;  // GP standard deviation parameters
  array[Kgp_1] vector<lower=0>[1] lscale_1;  // GP length-scale parameters
  vector[NBgp_1] zgp_1;  // latent variables of the GP
  real Intercept_sigma;  // temporary intercept for centered predictors
  vector<lower=0>[Kgp_sigma_1] sdgp_sigma_1;  // GP standard deviation parameters
  array[Kgp_sigma_1] vector<lower=0>[1] lscale_sigma_1;  // GP length-scale parameters
  vector[NBgp_sigma_1] zgp_sigma_1;  // latent variables of the GP
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, -13.3, 35.6);
  lprior += student_t_lpdf(sdgp_1 | 3, 0, 35.6)
    - 1 * student_t_lccdf(0 | 3, 0, 35.6);
  lprior += inv_gamma_lpdf(lscale_1[1][1] | 1.494197, 0.056607);
  lprior += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
  lprior += student_t_lpdf(sdgp_sigma_1 | 3, 0, 35.6)
    - 1 * student_t_lccdf(0 | 3, 0, 35.6);
  lprior += inv_gamma_lpdf(lscale_sigma_1[1][1] | 1.494197, 0.056607);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // scale latent variables of the GP
    vector[NBgp_1] rgp_1 = sqrt(spd_gp_exp_quad(slambda_1, sdgp_1[1], lscale_1[1])) .* zgp_1;
    vector[Nsubgp_1] gp_pred_1 = Xgp_1 * rgp_1;
    // scale latent variables of the GP
    vector[NBgp_sigma_1] rgp_sigma_1 = sqrt(spd_gp_exp_quad(slambda_sigma_1, sdgp_sigma_1[1], lscale_sigma_1[1])) .* zgp_sigma_1;
    vector[Nsubgp_sigma_1] gp_pred_sigma_1 = Xgp_sigma_1 * rgp_sigma_1;
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] sigma = rep_vector(0.0, N);
    mu += Intercept + gp_pred_1[Jgp_1];
    sigma += Intercept_sigma + gp_pred_sigma_1[Jgp_sigma_1];
    sigma = exp(sigma);
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zgp_1);
  target += std_normal_lpdf(zgp_sigma_1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
  // actual population-level intercept
  real b_sigma_Intercept = Intercept_sigma;
}

