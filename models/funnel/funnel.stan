// taken from https://mc-stan.org/docs/cmdstan-guide/diagnose_utility.html#running-the-diagnose-command
parameters {
  real y;
  vector[9] x;
}
model {
  y ~ normal(0, 3);
  x ~ normal(0, exp(y / 2));
}
