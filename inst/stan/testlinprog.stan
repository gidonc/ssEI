
functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan
  #include include/lkjonionfun.stan
  #include include/linprogfun.stan

}
data{
 int<lower=0> J;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 int<lower=0, upper=1> lflag_has_target;
 matrix<lower=0>[J, R] row_margins; // the row margins in each area
 matrix<lower=0>[J, C] col_margins; // the column margins in each area
 matrix<lower=0>[R, C] overall_tots; // the overall total target distribution (if known - used only if lflag_has_target==1)

}
transformed data{
  int nconst = J*R + J*C + R*C;
  int ncells = J*R*C;
  array[J, R, C] int cellpos;
  matrix[nconst, ncells + 1] bfs = rep_matrix(0, nconst, ncells + 1);
  int cellmap[3, ncells];
  matrix[nconst + 1, ncells + nconst + 1] bfsart = rep_matrix(0, nconst + 1, ncells + nconst + 1);
  int bfs_basis[ncells];
  int nweights = 0;
  int curcellpos = 0;
  int the_row = 0;
  for (j in 1:J){
    for (r in 1:R){
      for(c in 1:C){
        curcellpos += 1;
        cellpos[j, r, c] = curcellpos;
        cellmap[1, curcellpos] = j;
        cellmap[2, curcellpos] = r;
        cellmap[3, curcellpos] = c;
      }
    }
  }

  for(j in 1:J){
    for(r in 1:R){
      the_row += 1;
      bfs[the_row, ncells + 1] = row_margins[j, r];
      for(c in 1:C){
        bfs[the_row, cellpos[j, r, c]] = 1;
      }
    }
  }
  for(j in 1:J){
    for(c in 1:C){
      the_row += 1;
      bfs[the_row, ncells + 1] = col_margins[j, c];
      for(r in 1:R){
        bfs[the_row, cellpos[j, r, c]] = 1;
      }
    }
  }
  for(r in 1:R){
    for(c in 1:C){
      the_row += 1;
      bfs[the_row, ncells + 1] = overall_tots[r, c];
      for(j in 1:J){
        bfs[the_row, cellpos[j, r, c]] = 1;
      }
    }
  }
  print("starting bfs calc");
  bfsart = addart(bfs);
  bfsart = run_simplex_alg(bfsart);
  bfs = dropart(bfsart);
  print("finished bfs calc");

  for (i in 1:ncells){
    bfs_basis[i] = which_basic(col(bfs, i));
    if(bfs_basis[i]==0){
      nweights += 1;
    }
  }


}
parameters{
  vector[nweights] w;
}
transformed parameters{
  real cell_values[J, R, C];

  cell_values = ss_asign_from_bfs(bfs, bfs_basis, w, cellmap, J, R, C);

}

model{
  w ~ normal(0, 1);
}
generated quantities{


}
