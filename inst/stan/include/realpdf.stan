  real realpoisson_lpdf(row_vector x, vector lambda){
    // pdf for poisson extended to positive real numbers
    vector[num_elements(x)] lpdf = (to_vector(x) .* log(lambda)) - lambda - to_vector(lgamma(x + 1));
    return sum(lpdf);
  }
    real realnegbinom2_lpdf(vector x, real mu, real psi){
      vector[num_elements(x)] lpdf;
      for (n in 1:num_elements(x)){
        lpdf[n] = lchoose(x[n] + psi - 1, x[n]) + x[n]*log(mu) - x[n]*log(mu + psi) + psi*log(psi) - psi*log(mu + psi);

      }
    return sum(lpdf);
    }

    // real realnegbinom2_lpdf(row_vector x, vector mu, vector psi){
    //   vector[num_elements(x)] lpdf;
    //   for (n in 1:num_elements(x)){
    //     lpdf[n] = lchoose(x[n] + psi[n] - 1, x[n]) + x[n]*log(mu[n]) - x[n]*log(mu[n] + psi[n]) + psi[n]*log(psi[n]) - psi[n]*log(mu[n] + psi[n]);
    //
    //   }
//
//       return sum(lpdf);
//   }

  real realnegbinom3_lpdf(row_vector x, vector alpha, vector theta){
      vector[num_elements(x)] lpdf;
      vector[num_elements(x)] lpdf2;
      for (n in 1:num_elements(x)){
        lpdf[n] = lgamma(x[n] + alpha[n]) - lgamma(alpha[n]) - lgamma(x[n] + 1) + x[n]*log(theta[n]) + alpha[n]*log(1 - theta[n]);

      }

      return sum(lpdf);
  }


  real realmultinom_lpdf(matrix x, matrix theta){
    // pdf for multinomial extended to positive real numbers
    real  N[rows(x)];
    vector[rows(x)] lpdf;
    vector[rows(x)] lmcoef;
    int emts = rows(x);

    for (k in 1:emts){
      N[k] = sum(x[k]);
      lmcoef[k] = lgamma(N[k] +1 ) - sum(lgamma(x[k] + 1));
      lpdf[k] = lmcoef[k] + sum(x[k].*(log(theta[k])));
    }
    return sum(lpdf);
  }
  real realdirmultinom_lpdf(matrix x, matrix theta){
    // pdf for dirichlet-multinomial extended to positive real numbers
    real  N[rows(x)];
    real theta_0[rows(x)];
    vector[rows(x)] lpdf;
    vector[rows(x)] lmcoef;
    int emts = rows(x);

    for (k in 1:emts){
      N[k] = sum(x[k]);
      theta_0[k] = sum(theta[k]);
      if(theta_0[k] == 0) { // in case there are no counts leave loglikelihood unaltered by observations
        lmcoef[k] = 0;
        lpdf[k] = 0;
      } else {
        lmcoef[k] = lgamma(theta_0[k]) + lgamma(N[k] + 1) - lgamma(N[k] + theta_0[k]);
        lpdf[k] = lmcoef[k] + sum(lgamma(theta[k] + x[k]) - lgamma(theta[k]) - lgamma(x[k] + 1));
      }
    }
    return sum(lpdf);
  }
