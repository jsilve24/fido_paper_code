functions {
  // bijective softmax transformation
  vector softmax_id(vector alpha) {
    vector[num_elements(alpha) + 1] alphac1;
    for (k in 1:num_elements(alpha))
      alphac1[k] = alpha[k];
    alphac1[num_elements(alphac1)] = 0;
    return softmax(alphac1);
  }
  
  real matrix_normal_precomputed_lpdf(matrix X, matrix M, matrix U_inv, matrix V_inv, 
                                      real U_logdet, real V_logdet){
    int r = rows(X);
    int q = cols(X);
    matrix[r,q] Xresid;
    real tr;
    real log2pi = 1.8378770664093453;
    real numerator;
    real denominator;
    
    // precompute inverse and determinant
    Xresid = X-M;
    numerator = -0.5*trace_gen_quad_form(V_inv, U_inv, Xresid);
    denominator = -0.5*r*q*log2pi - 0.5*r*V_logdet - 0.5*q*U_logdet;
    
    return(numerator + denominator);
  }
  
  matrix matrix_normal_cholesky_rng(matrix M, matrix L_U, matrix L_V){
    int rp = rows(M);
    int q = cols(M);
    matrix[rp, q] X;
    
    for (i in 1:rp){
      for (j in 1:q){
        X[i,j] = normal_rng(0,1);
      }
    }
    // print("M:", M, "L_U", L_U, "L_V: ", L_V);
    // print("M+L_U*X*L_V'", M+L_U*X*L_V');
    return(M+L_U*X*L_V');
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
transformed data {
  int noffdiag = (D*D - 3*D + 2)/2;
  matrix[N, N] IN = diag_matrix(rep_vector(1, N));
  cov_matrix[D-1] XiInv = inverse_spd(Xi);
  cov_matrix[Q] GammaInv = inverse_spd(Gamma);
  cholesky_factor_cov[Q] LGammaInv = cholesky_decompose(GammaInv);
  real Gamma_logdet = log_determinant(Gamma);
}
parameters {
  matrix[D-1, Q] B;
  matrix[D-1, N] eta;
  // for bartlet decomposition of wishart (following wikipedia notation)
  vector<lower=0>[D-1] c;
  vector[noffdiag] n;
}
transformed parameters {
  simplex[D] pi[N];
  cholesky_factor_cov[D-1] LSigmaInv;
  
  // Bartlet Decomposition
  { matrix[D-1, D-1] W;
    W = diag_matrix(sqrt(c));
    if (D-1 > 1){
      int pos=1;
      for (j in 1:(D-1)){
        for (i in (j+1):(D-1)){
          W[i,j] = n[pos];
          pos=pos+1;
        }
      }  
      LSigmaInv = W;
    }
    LSigmaInv = LGammaInv * LSigmaInv;  
  }
  
  // Multinomial
  for (j in 1:N)
    pi[j] = softmax_id(eta[,j]);
}
model {
  real Sigma_logdet = 1/(prod(c)^2);
  // Priors 
  n ~ normal(0,1);
  for (i in 1:(D-1)){
    c[i] ~ chi_square(upsilon - i + 1);
  }
  {matrix[D-1, D-1] SigmaInv = tcrossprod(LSigmaInv);
   target += matrix_normal_precomputed_lpdf(B | Theta, SigmaInv, GammaInv, 
                                                Sigma_logdet, Gamma_logdet);
   target += matrix_normal_precomputed_lpdf(eta | B*X, SigmaInv, IN, 
                                                 Sigma_logdet, 1.0);
  }
 
  for(i in 1:N)
    Y[,i] ~ multinomial(pi[i]);
}
generated quantities {
  cov_matrix[D-1] Sigma;
  Sigma = inverse_spd(tcrossprod(LSigmaInv));
}
