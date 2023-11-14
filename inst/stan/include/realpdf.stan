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
