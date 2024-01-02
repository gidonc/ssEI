// implementing the onion method to avoid problems with larger correlation matricies
// onion method implementation described here: https://discourse.mc-stan.org/t/lkj-corr-lpdf-correlation-matrix-is-not-positive-definite/31995/5

  array[] vector create_shapes(int K, real eta) {
    real alpha = eta + (K - 2) / 2.0;
    array[2] vector[K - 1] shape;

    shape[1, 1] = alpha;
    shape[2, 1] = alpha;
    for (k in 2 : (K - 1)) {
      alpha -= 0.5;
      shape[1, k] = k / 2.0;
      shape[2, k] = alpha;
    }

    return shape;
  }

  matrix lkj_onion(int K, row_vector l, vector R2, data array[] vector shape) {
    matrix[K, K] L = rep_matrix(0, K, K); // cholesky_factor corr matrix
    {
      int start = 1;
      int end = 2;

      L[1, 1] = 1.0;
      L[2, 1] = 2.0 * R2[1] - 1.0;
      L[2, 2] = sqrt(1.0 - square(L[2, 1]));
      for (k in 2 : (K - 1)) {
        int kp1 = k + 1;
        row_vector[k] l_row = segment(l, start, k);
        real scale = sqrt(R2[k] / dot_self(l_row));
        L[kp1, 1 : k] = l_row[1 : k] * scale;
        L[kp1, kp1] = sqrt(1.0 - R2[k]);
        start = end + 1;
        end = start + k - 1;
      }
    }
    return L;
  }
