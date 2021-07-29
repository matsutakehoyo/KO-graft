data {
  int<lower=1> n;
  int<lower=1> k;
  vector[k] x[n];
}

parameters {
  vector[k] mu;
  vector<lower=0>[k] sigma;
  real<lower=1> nu; 
  cholesky_factor_corr[k] Lcorr;// cholesky factor (L_u matrix for R)
}

transformed parameters {
  corr_matrix[k] R; // correlation matrix
  cov_matrix[k] cov; // VCV matrix
  R = multiply_lower_tri_self_transpose(Lcorr);
  cov = quad_form_diag(R, sigma); // quad_form_diag: diag_matrix(sigma) * R * diag_matrix(sigma)
}

model {
  x ~ multi_student_t(nu, mu, cov);
  // x ~ multi_normal(mu, cov);
  sigma ~ normal(0,1);
  Lcorr ~ lkj_corr_cholesky(2);
  mu ~ normal(0,1);
  nu ~ gamma(2,.1);
}

