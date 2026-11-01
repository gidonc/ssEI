// Rosen, Jiang, King & Tanner (2001) multinomial-Dirichlet model, as in eiPack::ei.MD.bayes (no covariate), for HMC.
// Unit j, row r: beta[j, r] ~ Dirichlet(alpha[r]); alpha[r, c] ~ Gamma(lambda1, lambda2) (shape, rate; eiPack default 4, 2).
// Column counts: Y[j] ~ Multinomial(sum_r X[j, r] beta[j, r]), X = row proportions (eiPack instead augments the cell counts; same mean).
data {
  int<lower=1> J; int<lower=2> R; int<lower=2> C;
  matrix[J, R] X;
  array[J, C] int<lower=0> Y;
  real<lower=0> lambda1; real<lower=0> lambda2;
}
parameters {
  array[R] vector<lower=0>[C] alpha;
  array[J, R] simplex[C] beta;
}
model {
  for (r in 1:R) alpha[r] ~ gamma(lambda1, lambda2);
  for (j in 1:J) {
    vector[C] th = rep_vector(0, C);
    for (r in 1:R) { beta[j, r] ~ dirichlet(alpha[r]); th += X[j, r] * beta[j, r]; }
    Y[j] ~ multinomial(th / sum(th));
  }
}
