
      real robust_hinge_min(vector x, real delta){
        return(-1*delta*log_sum_exp(-1 * x/delta));
      }

      real robust_hinge_floor_zero(real x, real delta) {
        return(delta*log1p_exp(x/delta));
      }

  matrix make_helmert_basis(int N) {
    matrix[N, N - 1] V = rep_matrix(0.0, N, N - 1);
    for (k in 1:(N - 1)) {
      real r_k = k;
      real norm_const = sqrt(r_k * (r_k + 1.0));
      for (i in 1:k) {
        V[i, k] = 1.0 / norm_const;
      }
      V[k + 1, k] = -r_k / norm_const;
    }
    return V;
  }
  matrix kronecker_prod(matrix A, matrix B) {
    int rA = rows(A);
    int cA = cols(A);
    int rB = rows(B);
    int cB = cols(B);
    matrix[rA * rB, cA * cB] C;

    for (i in 1:rA) {
      for (j in 1:cA) {
        for (k in 1:rB) {
          for (l in 1:cB) {
            C[(i - 1) * rB + k, (j - 1) * cB + l] = A[i, j] * B[k, l];
          }
        }
      }
    }
    return C;
  }



  // ---------------------------------------------------------------------------------------------
  // Bounds of ONE cell in a sequential allocation -- used by every allocation in this model
  // (realised table, expected table, and the inverse used for anchors).
  //   sr   = what is left of the row,   sc = what is left of this column,
  //   rest = what is left of the later columns (they must be able to absorb the rest of the row).
  // True bounds:  lower = max(0, sr - rest),  upper = min(sc, sr),
  //               so the true range is  min(sc, sr, rest, sc + rest - sr).
  // Returns [lower, width, width] where the cell is lower + p * width, 0 < p < 1.
  //
  // delta > 0 : KINK SMOOTHING. Each max / min is replaced by a smooth version, with the corner
  //   rounded over a distance of about  delta x (the cell's own range).  Properties, for any inputs:
  //     * inward only: lower >= true lower, lower + width <= true upper  -> every table is feasible,
  //       margins are met exactly, and no slack tolerance is needed;
  //     * width > 0 whenever the true range is > 0;
  //     * at most 1.4 x delta of the cell's true range is given up (and only near a tie);
  //     * scale free: the same delta means the same thing in shares or in counts, and for small rows.
  // delta = 0 : the previous behaviour (hinges with legacy_floor / legacy_min), unchanged.
  vector seq_cell_bounds(real sr, real sc, real rest, real delta, real legacy_floor, real legacy_min) {
    vector[3] out;
    if (delta > 0) {
      real k = 1 / delta;
      real l_sr = log(fmax(sr, 1e-300));
      real l_sc = log(fmax(sc, 1e-300));
      real l_rest = log(fmax(rest, 1e-300));
      real l_oth = log(fmax(sc + rest - sr, 1e-300));
      // smooth version (from below, always > 0) of the true range min(sc, sr, rest, sc + rest - sr)
      real l_w = -log_sum_exp([-k * l_sc, -k * l_sr, -k * l_rest, -k * l_oth]') / k;
      // sharpness of the rounding: corner rounded over (range x delta), so narrow ranges get sharper corners
      real ke = k * exp(fmin(l_sr - l_w, 30));
      // upper = smooth min(sc, sr) from below;  lower = sr - smooth min(sr, rest) = smooth max(0, sr - rest) from above
      real up = sc * exp(-log1p_exp(ke * (l_sc - l_sr)) / ke);
      real lo = sr * (-expm1(-log1p_exp(ke * (l_sr - l_rest)) / ke));
      out[1] = lo;
      out[2] = up - lo;
      out[3] = up - lo;
    } else {
      real lo = robust_hinge_floor_zero(sr - rest, legacy_floor);
      real up = robust_hinge_min([sc, sr]', legacy_min);
      out[1] = lo;
      out[2] = robust_hinge_floor_zero(up - lo, legacy_min);
      out[3] = up - lo;
    }
    return out;
  }
  // share of the TRUE range [max(0, sr - rest), min(sc, sr)] that the smoothed cell cannot reach (0 if sharp)
  real seq_cell_range_loss(real sr, real sc, real rest, real width) {
    real tw = fmin(sc, sr) - fmax(0, sr - rest);
    return tw > 0 ? fmax(1 - width / tw, 0) : 0;
  }


  // ---------------------------------------------------------------------------------------------
  // ODDS-RATIO placement of one cell in a sequential allocation (kink free, no smoothing setting).
  // The cell x sits in the 2x2 table   [ x        sr - x          ]   this row
  //                                    [ sc - x   rest - sr + x   ]   the remaining rows
  //                                      this col   later cols
  // and lam is its log odds ratio:  lam = log x + log(rest - sr + x) - log(sr - x) - log(sc - x).
  // x is the root of a quadratic; it is always strictly inside max(0, sr - rest) < x < min(sc, sr),
  // covers that whole range as lam runs over the real line, and is smooth in sr, sc, rest.
  // lam = 0 gives independence, x = sr * sc / (sc + rest).
  real seq_or_cell(real sr, real sc, real rest, real lam) {
    real psi = exp(lam);
    real A = psi - 1;
    real B = -(psi * (sr + sc) + rest - sr);
    real Cq = psi * sr * sc;
    real disc = sqrt(fmax(B * B - 4 * A * Cq, 0));
    if (B <= 0) return 2 * Cq / (-B + disc);       // stable form of the root inside the bounds
    return (-B - disc) / (2 * A);
  }
  // log |dx / dlam| for that cell
  real seq_or_logjac(real x, real sr, real sc, real rest) {
    return -log(inv(x) + inv(sr - x) + inv(sc - x) + inv(rest - sr + x));
  }


  // ---------------------------------------------------------------------------------------------
  // ROW-WISE allocation: the common shift tau that makes a row add up.
  // The row takes cap[c] * inv_logit(eta[c] + tau) from column c; tau solves sum_c (that) = total.
  // The left side is increasing in tau, so this is a one-dimensional monotone solve (safeguarded Newton).
  // log of the integral over x >= 0 of the continuous Poisson density exp(x * log(mu) - mu - lgamma(x + 1)), as a function of
  // t = log(mu): cubic Hermite interpolation on a uniform grid (values v, derivatives d). The integral is 1 for large means
  // (zero returned above the grid) and falls below 1 for small ones (0.83 at mu = 1, 0.41 at mu = 0.1).
  real cpois_lognorm(real t, real t0, real h, vector v, vector d) {
    int n = rows(v);
    int lo = 1;
    int hi = n;
    real s;
    if (t >= t0 + h * (n - 1)) return 0;
    if (t <= t0) return v[1] + d[1] * (t - t0);
    while (hi - lo > 1) {
      int mid = (lo + hi) %/% 2;
      if (t >= t0 + h * (mid - 1)) lo = mid; else hi = mid;
    }
    s = (t - (t0 + h * (lo - 1))) / h;
    return (2 * s^3 - 3 * s^2 + 1) * v[lo] + (s^3 - 2 * s^2 + s) * h * d[lo] + (-2 * s^3 + 3 * s^2) * v[lo + 1] + (s^3 - s^2) * h * d[lo + 1];
  }
  real seq_row_shift(vector cap, vector eta, real total) {
    real base = logit(total / sum(cap));
    real lo = base - max(eta);                       // f(lo) <= 0
    real hi = base - min(eta);                       // f(hi) >= 0
    real tau = fmin(fmax(log(total) - log_sum_exp(log(cap) + eta), lo), hi);   // exact when the fractions are small
    for (it in 1:80) {
      vector[rows(cap)] p = inv_logit(eta + tau);
      real f = dot_product(cap, p) - total;
      real fp = dot_product(cap, p .* (1 - p));
      real nt = tau - f / fp;
      // converged to 1e-8: the Newton step from here is accurate to rounding error (quadratic convergence) and gives the gradient
      if (fabs(f) < 1e-8 * total) { tau = nt; break; }
      if (f > 0) hi = tau; else lo = tau;
      tau = (nt > lo && nt < hi) ? nt : 0.5 * (lo + hi);
    }
    return tau;
  }
  // ---------------------------------------------------------------------------------------------
  // Sequential allocation of the expected table with a chosen ORDER.
  //   row_order[1:(R-1)] = rows allocated, in order; row_order[R] = the remainder row.
  //   rem_col[r]         = the column of row r that is found by subtraction (cell-wise) / used as reference (row-wise).
  //   mode 2 = cell by cell, each cell placed by its 2x2 log odds ratio (seq_or_cell)
  //   mode 3 = row by row: row r takes cap[c] * inv_logit(lam_c + tau) from each column, tau from seq_row_shift
  // lam is laid out row by row in allocation order, C-1 values per row (columns other than rem_col, in natural order).
  // Returns (R+1) x C: rows 1:R = table, [R+1, 1] = log|d free cells / d lam|.
  matrix seq_alloc_ordered(vector w, vector m, vector lam, int mode, int[] row_order, int[] rem_col) {
    int R = rows(w);
    int C = rows(m);
    vector[R] sr = w;
    vector[C] sc = m;
    matrix[R + 1, C] out = rep_matrix(0, R + 1, C);
    real lj = 0;
    for (k in 1:(R - 1)) {
      int r = row_order[k];
      int last = rem_col[r];
      int i = 0;
      if (mode == 3) {
        vector[C] eta = rep_vector(0, C);
        vector[C] x;
        vector[C] v;
        real tau;
        for (c in 1:C) if (c != last) { i += 1; eta[c] = lam[(k - 1) * (C - 1) + i]; }
        tau = seq_row_shift(sc, eta, sr[r]);
        x = sc .* inv_logit(eta + tau);
        v = x .* (1 - x ./ sc);
        for (c in 1:C) out[r, c] = x[c];
        lj += sum(log(v)) - log(sum(v));
        sc -= x;
      } else {
        real rest = sum(sc);                         // capacity of the columns this row has not visited yet
        for (c in 1:C) if (c != last) {
          real x;
          i += 1;
          rest -= sc[c];
          x = seq_or_cell(sr[r], sc[c], rest, lam[(k - 1) * (C - 1) + i]);
          out[r, c] = x;
          lj += seq_or_logjac(x, sr[r], sc[c], rest);
          sc[c] -= x;
          sr[r] -= x;
        }
        out[r, last] = sr[r];
        sc[last] -= sr[r];
      }
    }
    for (c in 1:C) out[row_order[R], c] = sc[c];
    out[R + 1, 1] = lj;
    return out;
  }
  // the lam that reproduce table T under seq_alloc_ordered (T's own margins)
  vector seq_inv_ordered(matrix T, int mode, int[] row_order, int[] rem_col) {
    int R = rows(T);
    int C = cols(T);
    vector[R] sr;
    vector[C] sc;
    vector[(R - 1) * (C - 1)] lam;
    for (r in 1:R) sr[r] = sum(T[r]);
    for (c in 1:C) sc[c] = sum(col(T, c));
    for (k in 1:(R - 1)) {
      int r = row_order[k];
      int last = rem_col[r];
      int i = 0;
      if (mode == 3) {
        for (c in 1:C) if (c != last) {
          i += 1;
          lam[(k - 1) * (C - 1) + i] = logit(T[r, c] / sc[c]) - logit(T[r, last] / sc[last]);
        }
        sc -= T[r]';
      } else {
        real rest = sum(sc);
        for (c in 1:C) if (c != last) {
          real x = T[r, c];
          i += 1;
          rest -= sc[c];
          lam[(k - 1) * (C - 1) + i] = log(x) + log(rest - sr[r] + x) - log(sr[r] - x) - log(sc[c] - x);
          sc[c] -= x;
          sr[r] -= x;
        }
        sc[last] -= sr[r];
      }
    }
    return lam;
  }

      real[,,] ss_assign_cvals_lp (
        int n_areas, int R, int C,
        matrix row_margins, matrix col_margins,
        real[,,] lambda,
        real delta_floor, real delta_min, real slack_tol, real delta_rel){
    // constrains using sequential sampling approach described in Chen et. al 2005
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure. The function returns a flattened version of this information, together with an additional value, which is the contains the log determinant of the transform from the unconstrained lambda parameter to the cell values.
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;
     real cell_value[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;

     // slack_row=row_margins;
     // slack_col=col_margins;

     lower_pos[1]=0.0;
     log_det_J = 0;
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);

       free_R = 0;
       free_C = 0;

       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }
       array[free_R] int active_row_map;
       array[free_C] int active_col_map;
       { int t = 0; for (r in 1:R) if (row_margins[j,r]>0) { t += 1; active_row_map[t] = r; } }
       { int t = 0; for (c in 1:C) if (col_margins[j,c]>0) { t += 1; active_col_map[t] = c; } }


       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];

      rt=sum(slack_row);
      matrix[free_R, free_C] tmp_cell_value;

       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C -1 )){
           vector[3] bw = seq_cell_bounds(slack_row[r], slack_col[c], sum(tail(slack_col, free_C-c)), delta_rel, delta_floor, delta_min);
           lower_bound = bw[1];
           if (bw[3] < -slack_tol) reject("Negative cell width:", bw[3]);

           real x = lambda[j, r, c];
           real width = bw[2];
           this_inv_logit = inv_logit(x);

           tmp_cell_value[r,c]= lower_bound + this_inv_logit*width;

           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           if (slack_col[c] < -slack_tol) reject("Substantial negative column slack:", slack_col[c]);
           slack_col[c] = fmax(slack_col[c], 0.0);

           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           if (slack_row[r] < -slack_tol)  reject("Substantial negative row slack:", slack_row[r]);
           slack_row[r] = fmax(slack_row[r], 0.0);

           rt = rt - tmp_cell_value[r, c];
           if (rt < -slack_tol)  reject("Substantial negative total slack:", rt);
           rt = fmax(rt, 0.0);

           log_det_J += log(width) + log_inv_logit(x) + log1m_inv_logit(x);

           // tmp_cell_value[r,c]= fmax(0.0, lower_bound + this_inv_logit*(upper_bound-lower_bound));
           // slack_col[c]=fmax(0.0, slack_col[c] - tmp_cell_value[r,c]);
           // slack_row[r]=fmax(0.0, slack_row[r] - tmp_cell_value[r,c]);
           // rt = fmax(0.0, rt - tmp_cell_value[r, c]);
           // log_det_J += log(fmax(upper_bound - lower_bound, 1e-15)) + log_inv_logit(x) + log1m_inv_logit(x);
           // log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }

         tmp_cell_value[r, free_C]=slack_row[r];

         rt = rt - tmp_cell_value[r, free_C];
         if(rt < -slack_tol) reject("Substantial negative total slack:", rt);
         rt = fmax(rt, 0.0);

         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         if(slack_col[free_C] < -slack_tol) reject("Substantial negative column residual:", slack_col[free_C]);
         slack_col[free_C] = fmax(slack_col[free_C], 0.0);

         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
         if (slack_row[r] < -slack_tol) reject("Substantial negative final row slack:", slack_row[r]);
         if (slack_row[r] > slack_tol) reject("Substantial positive row slack after allocation should be complete:", slack_row[r]);
         slack_row[r] = fmax(slack_row[r], 0.0);
       }

       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt- tmp_cell_value[free_R, c];
         if (rt < -slack_tol) reject("Substantial negative total in final cell:", rt);
         rt = fmax(rt, 0.0);

         slack_col[c] = fmax(slack_col[c] - tmp_cell_value[free_R, c], 0.0);
         slack_row[free_R] = fmax(slack_row[free_R] - tmp_cell_value[free_R, c], 0.0);

         // tmp_cell_value[free_R, c] = slack_col[c];
         // rt = fmax(0.0, rt- tmp_cell_value[free_R, c]);
         // slack_col[c] = fmax(0.0, slack_col[c] - tmp_cell_value[free_R, c]);
         // slack_row[free_R] = fmax(0.0, slack_row[free_R] - tmp_cell_value[free_R, c]);
       }
       tmp_cell_value[free_R, free_C]=fmax(0.0, rt);

      int fr = 0;
      for (r in 1:R){
       if(row_margins[j, r]>0){
         fr += 1;
       }
       int fc = 0;
       for (c in 1:C){
         if(col_margins[j, c]>0 && row_margins[j, r]>0){
           fc += 1;
           cell_value[j, r, c] = tmp_cell_value[fr, fc];
         } else {
           cell_value[j, r, c] = 0;
         }
       }

     }
     }

    // Jacobian
    target += log_det_J;
    return cell_value;

  }


matrix inverse_alloc(matrix T, real dfloor, real dmin, real delta_rel, int use_or) {
  int R = rows(T);
  int C = cols(T);
  row_vector[R] sr;
  row_vector[C] sc;
  for (r in 1:R) sr[r] = sum(T[r]);
  for (c in 1:C) sc[c] = sum(col(T, c));

  matrix[R - 1, C - 1] lam;
  for (r in 1:(R - 1)) {
    for (c in 1:(C - 1)) {
      vector[3] bw = seq_cell_bounds(sr[r], sc[c], sum(sc[(c + 1):C]), delta_rel, dfloor, dmin);
      real p = fmin(fmax((T[r, c] - bw[1]) / fmax(bw[2], 1e-12), 1e-8), 1 - 1e-8);
      lam[r, c] = logit(p);
      if (use_or == 1) {
        real rest = sum(sc[(c + 1):C]);
        lam[r, c] = log(fmax(T[r, c], 1e-300)) + log(fmax(rest - sr[r] + T[r, c], 1e-300))
                    - log(fmax(sr[r] - T[r, c], 1e-300)) - log(fmax(sc[c] - T[r, c], 1e-300));
      }
      sc[c] = fmax(sc[c] - T[r, c], 0);
      sr[r] = fmax(sr[r] - T[r, c], 0);
    }
    sc[C] = fmax(sc[C] - sr[r], 0);
    sr[r] = 0;
  }
  return lam;
}
