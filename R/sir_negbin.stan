functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real E = y[2];
      real I = y[3];
      real R = y[4];
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real epsilon = x_r[1];


      real dS = -beta * I * S / N;
      real dE =  beta * I * S / N - epsilon * E;
      real dI =  epsilon * E      - gamma * I;
      real dR =  gamma * I;
      
      return {dS, dE, dI, dR};
  }
}
data {
  int<lower=1> n_days;
  real y0[4];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
  real epsilon;
}
transformed data {
  real x_r[1];
  int x_i[1] = { N };
  x_r[1] = epsilon;
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}
transformed parameters{
  real y[n_days, 4];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  cases ~ neg_binomial_2(col(to_matrix(y), 2), phi);
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2), phi);
}

