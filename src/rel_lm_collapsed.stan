functions {
  // bijective softmax transformation
  vector softmax_id(vector alpha) {
    vector[num_elements(alpha) + 1] alphac1;
    for (k in 1:num_elements(alpha))
      alphac1[k] = alpha[k];
    alphac1[num_elements(alphac1)] = 0;
    return softmax(alphac1);
  }
}
data {
  int<lower=1> N; // number of samples in dataset
  int<lower=2> D; // number of species (e.g., OTUs in dataset)
  int<lower=1> Q; // number of covariates
  int<lower=0> Y[D, N]; // counts
  matrix[Q, N] X;

  // Priors
  matrix[D-1, Q] Theta;
  cov_matrix[Q] Gamma;
  cov_matrix[D-1] Xi;
  real<lower=D-2> upsilon;
}
transformed data{
  int<lower=0> YT[N, D];
  matrix[N, Q] XT = X';
  cov_matrix[D-1] K = inverse_spd(Xi);
  matrix[Q,D-1] ThetaT = Theta';
  matrix[N, N] I = diag_matrix(rep_vector(1, N));
  matrix[N, N] A = inverse_spd(I + XT*Gamma*XT');
  matrix[D-1, N] XTB0T = (XT*ThetaT)';
  real upsilonN = (upsilon+N+D-2)/2; // updated from (upsilon+N)/2
  // For Generated Quantities Block
  real nuN = upsilon+N;
  matrix[Q, Q] Lambda0 = inverse_spd(Gamma);
  matrix[Q, Q] LambdaNInv = inverse_spd(XT'*XT + Lambda0);
  matrix[Q, Q] LLambdaNInv = cholesky_decompose(LambdaNInv);
  for (i in 1:D){
    for (j in 1:N)
      YT[j,i] = Y[i,j];
  }
}
parameters {
  //matrix<lower=-30, upper=30>[D-1, N] eta;
  matrix[D-1, N] eta;
}
transformed parameters {
  matrix[D-1, N] E_T;
  simplex[D] pi[N];
  E_T = eta - XTB0T;
  for (j in 1:N)
    pi[j] = softmax_id(eta[,j]);
}
model {
  //target += -upsilonN*log_determinant(I + A*quad_form(Xi, E_T));
  target += -upsilonN*log_determinant(I + A*quad_form(K, E_T));
  for(i in 1:N)
    YT[i] ~ multinomial(pi[i]);
}
generated quantities {
  matrix[N, D-1] eta_T = eta';
  matrix[Q, D-1] BN = LambdaNInv*(XT'*eta_T + Lambda0*ThetaT);
  matrix[D-1, D-1] VN = Xi + crossprod(eta_T-XT*BN) + quad_form(Lambda0, BN-ThetaT);
  cov_matrix[D-1] Sigma = inv_wishart_rng(nuN, VN);
  matrix[Q, D-1] eB;
  matrix[Q, D-1] BT;
  matrix[D-1, Q] B;
  for (i in 1:Q){
    for (j in 1:(D-1))
      eB[i,j] = normal_rng(0,1);
  }
  BT = BN + LLambdaNInv*eB*cholesky_decompose(Sigma)';
  B = BT';
}
