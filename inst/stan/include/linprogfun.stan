
matrix lpivot (matrix m, int pivrow, int pivcol){

  real pivval = m[pivrow, pivcol];
  real pivremval;

  m[pivrow] = m[pivrow]/pivval;

  for (r in 1:rows(m)){
    if (r != pivrow){
      pivremval = mat[r, pivcol];
      m[r] = m[r] - pivremval * m[r];
    }
  }
  return(m);
}

matrix addart(matrix m){
  matrix mout;
  row_vector objfun;
  rm = rows(mat);
  norig = cols(mat);
  art = diag_matrix(rm);
  mout = append_col(m[1:rm, 1:(norig - 1)], art);
  mout = append_col(mout, m[1:rm, norig]);

  // add row with phase 1 objective function
  objfun = append_col(rep_row_vector(0, norig - 1), rep_row_vector(1, rm));
  objfun = append_col(objfun, 0);
  mout = append_row(mout, objfun);

  return(mout);
}

matrix dropart(matrix m){
  matrix mout;
  int n;
  // final row is objective function so drop that
  mout = m[1:(rows(m) - 1), 1:cols(m)];

  rm = rows(mout);
  n = ncol(mout) - rm;

  // the n non-artifical variables are in colums 1:(n - 1) and the final column

  mout = append_col(mout[1:rm, 1:(n - 1)], mout[1:rm, cols(mout)]);
  return(mout);
}

int mrtidx(vector avec, bvec){
  // returns the index of the row which meets the minimum ratio test
  // if there are only negative values it returns 0 (this would mean that the problem is unbounded)
  // if there are multiple equal max negative values it returns the index of the first occurance n the row_vector
  int idx = 0;
  real curmrt = 0;
  real thismrt;
  for (i in in 1:num_elements(avec)){
    if(avec[i] > 0){
      thismrt = bvec[i]/avec[i];
      if(thismrt < curmrt || curmrt == -1){
        curmrt = thismrt;
        idx = i;
      }
    }
  }
  return(idx);
}

int minnegidx(row_vector row_vec){
  // returns the index of the minimum negative value of a row_vector (excluding zero)
  // if there are no negative values it returns 0
  // if there are multiple equal minimum negative values it returns the index of the first occurance n the row_vector
  int idx = 0;
  real curmin = 0;
  for (i in in 1:num_elements(row_vec)){
    if(row_vec[i] < curmin){
      curmin = row_vec[i];
      idx = i;
    }
  }
  return(idx);
}

matrix run_simplex_alg(matrix m){
  int rpiv;
  int cpiv;

  while(cpiv > 0){
    cpiv = minnegidx(row(m, rows(m)));
    rpiv = mrtidx(col(m, cpiv), col(m, cols(m)));
    if(rpiv == 0) {reject("minimum ratio test found no positive values. Indicates an illegal unbounded maximization problem. Current matrix values are :"), mat};
    m = lpivot(m, rpiv, cpiv);
  }
  return(m);
}

int which_basic(vector vec){
  // function to find the index of the row that gives the value of a column which is in the basis of a linear programming solution
  // returns 0 if the column does not belong to a variable which is in the basis
  int not_zero[num_elements(vec)];
  int is_basic;
  int idx;

  not_zero = vec != 0;
  idx == 0;
  is_basic = sum(not_zero) == 1;
  if(is_basic == 1){
    for(i in 1:num_elements(not_zero)){
      if(not_zero[i] != 0){
        idx = i
      }
    }
  }
  return(idx);
}

vector get_var_lims(matrix m, int cm){
  // find limits for one variable from a matrix m giving a feasible starting point. Uses simplex algorithm with objective fuction which contains just one variable where the function to be maximized in f(x) =x (gives max) and f(x) = -x (gives min)
  int nvars = cols(m) - 1;
  int is_basic;
  int not_zero;
  vector var_lims;
  matrix mat_min;
  matrix mat_max;
  row_vector obj_max = rep_row_vector(0, nvars + 1);
  row_vector obj_min = rep_row_vector(0, nvars + 1);

  obj_max[cm] = -1;
  obj_min[cm] = 1;

  use_row = which_basic(col(m, cm));
  if(use_row > 0){
    obj_max = obj_max + m[use_row];
    obj_min = obj_min - m[use_row];
  }
  mat_min = run_simplex_alg(append_row(m, obj_min));
  mat_min = run_simplex_alg(append_row(m, obj_max));

  use_row = which_basic(col(mat_min, cm));
  if(use_row >0){
    var_lims[1] = mat_min[use_row, ncol(mat_min)]/mat_min[use_row, cm];
  } else{
    var_lims[1] = 0;
  }

  use_row = which_basic(col(mat_max, cm));
  if(use_row >0){
    var_lims[2] = mat_max[use_row, ncol(mat_max)]/mat_max[use_row, cm];
  } else{
    var_lims[1] = 0;
  }

  return(var_lims);
}
