#bw_cv----
#f is the number of folds
#func is either SEC_only_calibration or ESC_only_calibration
bw_cv = function(bw_min, bw_max, step, f, func, method, Xc, Xe, Xs, y,intercept, utm_ev_sp, utm_st_sp){
  bandwidths = seq(bw_min, bw_max, step)
  if (bandwidths[end(bandwidths)[1]]<bw_max-0.0001){
    bandwidths = c(bandwidths, bw_max)
  }
  cvss = matrix(0,length(bandwidths), f)
  batch = floor(length(y)/f)
  rest = length(y) - batch * f
  rest_const = rest
  for (w in 1:length(bandwidths)){
    i = bandwidths[w]
    print(c("Bandwidth", i))
    for (b in 1:f){
      if (rest >0){
        batch_t = batch + 1
        beginning = batch_t*(b-1)+1
        end = b*batch_t
        rest = rest - 1
      } else {
        beginning = rest_const + batch*(b-1) + 1 #(batch+1)*rest_const + batch*(b-1-rest_const)+1 
        end = rest_const + batch*b
      }
      #create subsets
      Xc_temp = as.matrix(Xc)[-(beginning:end),]
      Xe_temp = as.matrix(Xe)[-(beginning:end),]
      Xs_temp = as.matrix(Xs)[-(beginning:end),]
      y_temp = y[-(beginning:end)]
      coords_e_temp = utm_ev_sp[-(beginning:end),]
      coords_s_temp = utm_st_sp[-(beginning:end),]
      
      Xc_test = as.matrix(Xc)[(beginning:end),]
      Xe_test = as.matrix(Xe)[(beginning:end),]
      Xs_test = as.matrix(Xs)[(beginning:end),]
      y_test = y[(beginning:end)]
      coords_e_test = utm_ev_sp[(beginning:end),]
      coords_s_test = utm_st_sp[(beginning:end),]
      
      if (batch == 1){
        pcoords = c(coords_e_test, coords_s_test)
      } else {
        pcoords = cbind(coords_e_test, coords_s_test)
      }
      
      temp = func(Xc_temp, Xe_temp, Xs_temp, y_temp, intercept,
                  i, i, coords_e_temp, coords_s_temp)
      if (method == "sec"){
        prediction = emgwr_prediction_points_no_param(Xc_temp, Xe_temp, Xs_temp, y_temp, temp, "sec", i, i,
                                             coords_e_temp, coords_s_temp, Xc_test, Xe_test, Xs_test,
                                             pcoords, 0.05)
      } else {
        prediction = emgwr_prediction_points_no_param(Xc_temp, Xe_temp, Xs_temp, y_temp, temp, "esc", i, i,
                                             coords_e_temp, coords_s_temp, Xc_test, Xe_test, Xs_test,
                                             pcoords, 0.05)
      }
      
      
      
      #compute residuals
      y_res = rep(0,length(y_test))
      y_hat = prediction$y0
      y_res = (y_hat - y_test)
      cvss[w,b] = Norm(y_res,2)^2
    }
  }
  
  cvsum = matrix(0, length(bandwidths), 2)
  for (w in 1:length(bandwidths)){
    i = bandwidths[w]
    cvsum[w,] = c(i,sum(cvss[w,]))
  }
  best = min(cvsum[,2])
  best_bw = cvsum[(cvsum[,2]==best),1]
  
  solution = list("best_bw" = best_bw, 
                  "cvsum" = cvsum
  )
  return(solution)
}

#bw_gcv----
#GCV for ESC and SEC in EMGWR with constant intercept
#func is either SEC_only_calibration or ESC_only_calibration
bw_gcv = function(bw_min, bw_max, step, func, Xc, Xe, Xs, y,intercept, utm_ev_sp, utm_st_sp){
  bandwidths = seq(bw_min, bw_max, step)
  if (bandwidths[end(bandwidths)[1]]<bw_max-0.0001){
    bandwidths = c(bandwidths, bw_max)
  }
  
  gcv = rep(0,length(bandwidths))
  for (w in 1:length(bandwidths)){
    i = bandwidths[w]
    print(c("Bandwidth", i))
    temp = func(Xc, Xe, Xs, y, intercept, i, i, utm_ev_sp, utm_st_sp)
    B = temp$B
    n = length(y)
    I = diag(1,n)
    Xcc = cbind(rep(1,n), Xc)
    H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
    res = (I-H)%*%y
    rss = sum(res^2)
    delta1 = n-2*tr(H)+tr(t(H)%*%H)
    gcv[w] = n*rss/delta1^2
  }
  best = min(gcv)
  best_bw = bandwidths[gcv==best]
  
  solution = list("best_bw" = best_bw, 
                  "gcv" = gcv)
  return(solution)
}
#bw_gcv_mei----
#GCV for ESC and SEC in EMGWR with constant intercept
#func is either SEC_only_calibration or ESC_oly_calibration
bw_gcv_mei = function(bw_min, bw_max, step, func, Xc, Xe, Xs, y,intercept, utm_ev_sp, utm_st_sp){
  bandwidths = seq(bw_min, bw_max, step)
  if (bandwidths[end(bandwidths)[1]]<bw_max-0.0001){
    bandwidths = c(bandwidths, bw_max)
  }
  
  gcv = rep(0,length(bandwidths))
  for (w in 1:length(bandwidths)){
    i = bandwidths[w]
    print(c("Bandwidth", i))
    temp = func(Xc, Xe, Xs, y, intercept, i, i, utm_ev_sp, utm_st_sp)
    B = temp$B
    n = length(y)
    I = diag(1,n)
    Xcc = cbind(rep(1,n), Xc)
    H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
    res = (I-H)%*%y
    gcv[w] = sum(res^2/(rep(1, length(diag(H)))-diag(H))^2)
  }
  best = min(gcv)
  best_bw = bandwidths[gcv==best]
  
  solution = list("best_bw" = best_bw, 
                  "gcv" = gcv)
  return(solution)
}
#gcv_mei_only_one----
#GCV for ESC and SEC in EMGWR with constant intercept
#func is either SEC_only_calibration or ESC_only_calibration
gcv_mei_only_one = function(bwe, bws, func, Xc, Xe, Xs, y,intercept, utm_ev_sp, utm_st_sp){
  temp = func(Xc, Xe, Xs, y, intercept, bwe, bws, utm_ev_sp, utm_st_sp)
  B = temp$B
  n = length(y)
  I = diag(1,n)
  Xcc = cbind(rep(1,n), Xc)
  H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
  res = (I-H)%*%y
  gcv = sum(res^2/(rep(1, length(diag(H)))-diag(H))^2)
  
  solution = list("bwe" = bwe,
                  "bws" = bws,
                  "gcv" = gcv)
  return(solution)
}
#SEC_general----
SEC_general = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
    #print(c("Hs",i))
  }
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
    #print(c("He",i))
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_e for whole grid
  for (i in 1:L){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    beta_e[,i] = Ae %*% y_tilde
    #print(c("beta_e",i))
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-He)%*%y_tilde
  
  #5)compute beta_s for whole grid
  for (i in 1:L){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    beta_s[,i] = As %*% y_tilde_s
    #print(c("beta_s",i))
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}
#ESC_general----
ESC_general = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
    #print(c("He",i))
  }
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    Hs[i,] = Xs[i,] %*% As
    #print(c("Hs",i))
  }
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_s for whole grid
  for (i in 1:L){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    beta_s[,i] = As %*% y_tilde
    #print(c("beta_s",i))
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-Hs)%*%y_tilde
  
  #5)compute beta_e for whole grid
  for (i in 1:L){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    beta_e[,i] = Ae %*% y_tilde_s
    #print(c("beta_e",i))
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}
#SCE_general----
SCE_general = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  Hc = matrix(0,N,N)
  Hs = matrix(0,N,N)
  He = matrix(0,N,N)
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
  }
  #create Hc
  Ac = (solve((t(Xc)%*%t(I-Hs)%*%(I-Hs)%*%Xc))) %*% t(Xc) %*% t(I-Hs)
  Hc = Xc %*% Ac
  
  #create B
  B = I - Hc - Hs + Hs %*% Hc
  
  #1)compute beta_e for whole grid
  for (i in 1:(length(e1)*length(e2))){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve(t(Xe)%*%t(B)%*%We%*%B%*%Xe)) %*% t(Xe) %*% t(B) %*% We %*% B
    beta_e[,i] = Ae %*% y
  }
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve(t(Xe)%*%t(B)%*%We%*%B%*%Xe)) %*% t(Xe) %*% t(B) %*% We %*% B
    He[i,] = Xe[i,] %*% Ae
  }
  
  #2)compute y_tilde
  y_tilde = (I-He) %*% y
  
  #3)compute beta_c
  beta_c = Ac %*% y_tilde
  
  #4)compute y_tilde_s
  y_tilde_s = (I-Hc)%*%y_tilde
  
  #5)compute beta_s for whole grid
  for (i in 1:(length(s1)*length(s2))){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*%  Ws
    beta_s[,i] = As %*% y_tilde_s
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "Hc" = Hc,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}
#ECS_general----
ECS_general = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hc = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
  }
  #create Hc
  Ac = (solve((t(Xc)%*%t(I-He)%*%(I-He)%*%Xc))) %*% t(Xc) %*% t(I-He)
  Hc = Xc %*% Ac
  
  #create B
  B = I - Hc - He + He %*% Hc
  
  #1)compute beta_s for whole grid
  for (i in 1:(length(s1)*length(s2))){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve(t(Xs)%*%t(B)%*%Ws%*%B%*%Xs)) %*% t(Xs) %*% t(B) %*% Ws %*% B
    beta_s[,i] = As %*% y
  }
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve(t(Xs)%*%t(B)%*%Ws%*%B%*%Xs)) %*% t(Xs) %*% t(B) %*% Ws %*% B
    Hs[i,] = Xs[i,] %*% As
  }
  
  #2)compute y_tilde
  y_tilde = (I-Hs) %*% y
  
  #3)compute beta_c
  beta_c = Ac %*% y_tilde
  
  #4)compute y_tilde_s
  y_tilde_s = (I-Hc)%*%y_tilde
  
  #5)compute beta_e for whole grid
  for (i in 1:(length(e1)*length(e2))){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*%  We
    beta_e[,i] = Ae %*% y_tilde_s
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "Hc" = Hc,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}
#CSE_general----
CSE_general = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  Hc = matrix(0,N,N)
  Hs = matrix(0,N,N)
  He = matrix(0,N,N)
  
  #create Hc
  Ac = (solve((t(Xc)%*%Xc))) %*% t(Xc)
  Hc = Xc %*% Ac
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-Hc)))%*%Ws%*%(I-Hc)%*%Xs))) %*% t(Xs) %*% (t((I-Hc))) %*% Ws %*% (I-Hc)
    Hs[i,] = Xs[i,] %*% As
  }
  
  #create B
  B = I - Hc - Hs + Hc %*% Hs
  
  #1)compute beta_e for whole grid
  for (i in 1:(length(e1)*length(e2))){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve(t(Xe)%*%t(B)%*%We%*%B%*%Xe)) %*% t(Xe) %*% t(B) %*% We %*% B
    beta_e[,i] = Ae %*% y
  }
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve(t(Xe)%*%t(B)%*%We%*%B%*%Xe)) %*% t(Xe) %*% t(B) %*% We %*% B
    He[i,] = Xe[i,] %*% Ae
  }
  
  #2)compute y_tilde
  y_tilde = (I-He) %*% y
  
  #3)compute beta_s for whole grid
  for (i in 1:(length(s1)*length(s2))){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%(t((I-Hc)))%*%Ws%*%(I-Hc)%*%Xs))) %*% t(Xs) %*% (t((I-Hc))) %*% Ws %*% (I-Hc)
    beta_s[,i] = As %*% y_tilde
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-Hs)%*%y_tilde
  
  #5)compute beta_c
  beta_c = Ac %*% y_tilde_s
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "Hc" = Hc,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}
#CES_general----
CES_general = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hc = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hc
  Ac = (solve((t(Xc)%*%Xc))) %*% t(Xc)
  Hc = Xc %*% Ac
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hc)))%*%We%*%(I-Hc)%*%Xe))) %*% t(Xe) %*% (t((I-Hc))) %*% We %*% (I-Hc)
    He[i,] = Xe[i,] %*% Ae
  }
  
  #create B
  B = I - Hc - He + Hc %*% He
  
  #1)compute beta_s for whole grid
  for (i in 1:(length(s1)*length(s2))){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve(t(Xs)%*%t(B)%*%Ws%*%B%*%Xs)) %*% t(Xs) %*% t(B) %*% Ws %*% B
    beta_s[,i] = As %*% y
  }
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve(t(Xs)%*%t(B)%*%Ws%*%B%*%Xs)) %*% t(Xs) %*% t(B) %*% Ws %*% B
    Hs[i,] = Xs[i,] %*% As
  }
  
  #2)compute y_tilde
  y_tilde = (I-Hs) %*% y
  
  #3)compute beta_e for whole grid
  for (i in 1:(length(e1)*length(e2))){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hc)))%*%We%*%(I-Hc)%*%Xe))) %*% t(Xe) %*% (t((I-Hc))) %*% We %*% (I-Hc)
    beta_e[,i] = Ae %*% y_tilde
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-He)%*%y_tilde
  
  #5)compute beta_c
  beta_c = Ac %*% y_tilde_s
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "Hc" = Hc,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}

#emgwr_prediction_points----
#prediction for points with their coordinates
#Xc, Xe, Xs, y are data used for calibration
#emgwr is output of a calibration method is either esc or sec
#bwe, bws, utm_ev_sp, utm_st_sp are calibration coordinates
#pc, pe, ps prediction covariates
#pcoords are prediction coordinates (pass as matrix)
#alfa is confidence level for prediction interval
emgwr_prediction_points = function(Xc, Xe, Xs, y, emgwr, H_hat, delta1, delta2, method, bwe, bws, utm_ev_sp, utm_st_sp, pc, pe, ps, pcoords, alfa){
  
  N = length(y)
  K = dim(pcoords)[1]
  if (is.null(K)){
    K = 1
  }
  Xc = cbind(rep(1,N), Xc)
  if (K==1 | ((is.null(dim(pc)) & length(pc)!= K & length(pe)!= K & length(ps)!= K))){ #only one point to predict
    x0 = matrix(c(1, pc, pe, ps), ncol = 1)
  } else {
    x0 = cbind(1, pc, pe, ps)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  He = emgwr$He
  Hs = emgwr$Hs
  B = emgwr$B
  
  # H_hat = I - B + B %*% Xc %*% solve(t(Xc)%*%t(B)%*%B%*%Xc) %*% t(Xc) %*% t(B)%*% B
  res = (I-H_hat)%*%y
  # delta1 = n_sample-2*tr(H_hat)+tr(t(H_hat)%*%H_hat)
  # delta2 = tr((t(I-H_hat)%*%(I-H_hat)) %*% (t(I-H_hat)%*%(I-H_hat)))
  rss = sum(res^2)
  sigma2hat = rss/delta1
  sigmahat = sqrt(sigma2hat)
  
  #1)compute beta_c
  Ac = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B
  beta_c = Ac%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  if (K==1){
    dist_e_sim = gw.dist(utm_ev_sp, t(as.matrix(pcoords)[1:2]), focus=0, p=2, theta=0, longlat=F)
    dist_s_sim = gw.dist(utm_st_sp, t(as.matrix(pcoords)[3:4]), focus=0, p=2, theta=0, longlat=F)
    
    if (method == "esc") {
      Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
      As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
      beta_s = As %*% y_tilde
      print(c("beta_s",i))
      
      y_tilde_s = (I-Hs)%*%y_tilde
      
      We = diag(gauss_kernel(dist_e_sim[,1],bwe))
      Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
      beta_e = Ae %*% y_tilde_s
      print(c("beta_e",i))
      
      Q = rbind(Ac, Ae%*%(I-Hs)%*%(I-Xc%*%Ac), As%*%(I-Xc%*%Ac))
      
      s0 = t(x0)%*%Q%*%t(Q)%*%x0
      
      betas = c(beta_c, beta_e, beta_s)
      y0 = t(x0)%*%betas
      t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
      var0 = sigma2hat*s0
      interval = c(y0-sigmahat*sqrt(s0+1)*t_alfa2, y0, y0+sigmahat*sqrt(s0+1)*t_alfa2)
      
    } else {
      We = diag(gauss_kernel(dist_e_sim[,1],bwe))
      Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
      beta_e = Ae %*% y_tilde
      print(c("beta_e",i))
      
      y_tilde_s = (I-He)%*%y_tilde
      
      Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      beta_s = As %*% y_tilde_s
      print(c("beta_s",i))
      
      Q = rbind(Ac, Ae%*%(I-Xc%*%Ac), As%*%(I-He)%*%(I-Xc%*%Ac))
      
      s0 = t(x0)%*%Q%*%t(Q)%*%x0
      
      betas = c(beta_c, beta_e, beta_s)
      y0 = t(x0)%*%betas
      t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
      var0 = sigma2hat*s0
      interval = c(y0-sigmahat*sqrt(s0+1)*t_alfa2, y0, y0+sigmahat*sqrt(s0+1)*t_alfa2)
    }
  }
  
  else if (K > 1){
    betas = matrix(0,K,n_c+n_e+n_s)
    y0 = rep(0,K)
    var0 = rep(0,K)
    s0 = rep(0,K)
    interval = matrix(0,K,3)
    
    if (method == "esc"){
      for (k in 1:K){
        dist_e_sim = gw.dist(utm_ev_sp, coordinates(pcoords[k,1:2]), focus=0, p=2, theta=0, longlat=F)
        dist_s_sim = gw.dist(utm_st_sp, coordinates(pcoords[k,3:4]), focus=0, p=2, theta=0, longlat=F)
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
        As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
        beta_s = As %*% y_tilde
        print(c("beta_s",k))
        
        y_tilde_s = (I-Hs)%*%y_tilde
        
        We = diag(gauss_kernel(dist_e_sim[,1],bwe))
        Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
        beta_e = Ae %*% y_tilde_s
        print(c("beta_e",k))
        
        Q = rbind(Ac, Ae%*%(I-Hs)%*%(I-Xc%*%Ac),As%*%(I-Xc%*%Ac))
        
        s0[k] = t(x0[k,])%*%Q%*%t(Q)%*%x0[k,]
        
        betas[k,] = c(beta_c, beta_e, beta_s)
        y0[k] = x0[k,]%*%betas[k,]
        t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
        var0[k] = sigma2hat*s0[k]
        interval[k,] = c(y0[k]-sigmahat*sqrt(s0[k]+1)*t_alfa2, y0[k], y0[k]+sigmahat*sqrt(s0[k]+1)*t_alfa2)
      }
    } else {
      for (k in 1:K){
        dist_e_sim = gw.dist(utm_ev_sp, coordinates(pcoords[k,1:2]), focus=0, p=2, theta=0, longlat=F)
        dist_s_sim = gw.dist(utm_st_sp, coordinates(pcoords[k,3:4]), focus=0, p=2, theta=0, longlat=F)
        
        We = diag(gauss_kernel(dist_e_sim[,1],bwe))
        Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
        beta_e = Ae %*% y_tilde
        print(c("beta_e",k))
        
        y_tilde_s = (I-He)%*%y_tilde
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
        As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
        beta_s = As %*% y_tilde_s
        print(c("beta_s",k))
        
        Q = rbind(Ac, Ae%*%(I-Xc%*%Ac), As%*%(I-He)%*%(I-Xc%*%Ac))
        s0[k] = t(x0[k,])%*%Q%*%t(Q)%*%x0[k,]
        
        betas[k,] = c(beta_c, beta_e, beta_s)
        y0[k] = x0[k,]%*%betas[k,]
        t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
        var0[k] = sigma2hat*s0[k]
        interval[k,] = c(y0[k]-sigmahat*sqrt(s0[k]+1)*t_alfa2, y0[k], y0[k]+sigmahat*sqrt(s0[k]+1)*t_alfa2)
      }
    }
  }
  
  
  result <- list("s0" = s0, 
                 "y0" = y0,
                 "betas" = betas,
                 "interval" = interval,
                 "var0" = var0,
                 "sigma2hat" = sigma2hat)
  return(result)
}
#emgwr_prediction_points_no_param----
#prediction for points with their coordinates without already computed delta1, delta2 and H_hat
emgwr_prediction_points_no_param = function(Xc, Xe, Xs, y, emgwr, method, bwe, bws, utm_ev_sp, utm_st_sp, pc, pe, ps, pcoords, alfa){
  
  N = length(y)
  K = dim(pcoords)[1]
  if (is.null(K)){
    K = 1
  }
  Xc = cbind(rep(1,N), Xc)
  if (K==1 | ((is.null(dim(pc)) & length(pc)!= K & length(pe)!= K & length(ps)!= K))){ #only one point to predict
    x0 = matrix(c(1, pc, pe, ps), ncol = 1)
  } else {
    x0 = cbind(1, pc, pe, ps)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  He = emgwr$He
  Hs = emgwr$Hs
  B = emgwr$B
  
  H_hat = I - B + B %*% Xc %*% solve(t(Xc)%*%t(B)%*%B%*%Xc) %*% t(Xc) %*% t(B)%*% B
  res = (I-H_hat)%*%y
  delta1 = n_sample-2*tr(H_hat)+tr(t(H_hat)%*%H_hat)
  delta2 = tr((t(I-H_hat)%*%(I-H_hat)) %*% (t(I-H_hat)%*%(I-H_hat)))
  rss = sum(res^2)
  sigma2hat = rss/delta1
  sigmahat = sqrt(sigma2hat)
  
  #1)compute beta_c
  Ac = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B
  beta_c = Ac%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  if (K==1){
    dist_e_sim = gw.dist(utm_ev_sp, t(as.matrix(pcoords)[1:2]), focus=0, p=2, theta=0, longlat=F)
    dist_s_sim = gw.dist(utm_st_sp, t(as.matrix(pcoords)[3:4]), focus=0, p=2, theta=0, longlat=F)
    
    if (method == "esc") {
      Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
      As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
      beta_s = As %*% y_tilde
      print(c("beta_s",i))
      
      y_tilde_s = (I-Hs)%*%y_tilde
      
      We = diag(gauss_kernel(dist_e_sim[,1],bwe))
      Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
      beta_e = Ae %*% y_tilde_s
      print(c("beta_e",i))
      
      Q = rbind(Ac, Ae%*%(I-Hs)%*%(I-Xc%*%Ac), As%*%(I-Xc%*%Ac))
      
      s0 = t(x0)%*%Q%*%t(Q)%*%x0
      
      betas = c(beta_c, beta_e, beta_s)
      y0 = t(x0)%*%betas
      t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
      var0 = sigma2hat*s0
      interval = c(y0-sigmahat*sqrt(s0+1)*t_alfa2, y0, y0+sigmahat*sqrt(s0+1)*t_alfa2)
      
    } else {
      We = diag(gauss_kernel(dist_e_sim[,1],bwe))
      Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
      beta_e = Ae %*% y_tilde
      print(c("beta_e",i))
      
      y_tilde_s = (I-He)%*%y_tilde
      
      Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      beta_s = As %*% y_tilde_s
      print(c("beta_s",i))
      
      Q = rbind(Ac, Ae%*%(I-Xc%*%Ac), As%*%(I-He)%*%(I-Xc%*%Ac))
      
      s0 = t(x0)%*%Q%*%t(Q)%*%x0
      
      betas = c(beta_c, beta_e, beta_s)
      y0 = t(x0)%*%betas
      t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
      var0 = sigma2hat*s0
      interval = c(y0-sigmahat*sqrt(s0+1)*t_alfa2, y0, y0+sigmahat*sqrt(s0+1)*t_alfa2)
    }
  }
  
  else if (K > 1){
    betas = matrix(0,K,n_c+n_e+n_s)
    y0 = rep(0,K)
    var0 = rep(0,K)
    s0 = rep(0,K)
    interval = matrix(0,K,3)
    
    if (method == "esc"){
      for (k in 1:K){
        dist_e_sim = gw.dist(utm_ev_sp, pcoords[k,1:2], focus=0, p=2, theta=0, longlat=F)
        dist_s_sim = gw.dist(utm_st_sp, pcoords[k,3:4], focus=0, p=2, theta=0, longlat=F)
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
        As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
        beta_s = As %*% y_tilde
        print(c("beta_s",k))
        
        y_tilde_s = (I-Hs)%*%y_tilde
        
        We = diag(gauss_kernel(dist_e_sim[,1],bwe))
        Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
        beta_e = Ae %*% y_tilde_s
        print(c("beta_e",k))
        
        Q = rbind(Ac, Ae%*%(I-Hs)%*%(I-Xc%*%Ac),As%*%(I-Xc%*%Ac))
        
        s0[k] = t(x0[k,])%*%Q%*%t(Q)%*%x0[k,]
        
        betas[k,] = c(beta_c, beta_e, beta_s)
        y0[k] = x0[k,]%*%betas[k,]
        t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
        var0[k] = sigma2hat*s0[k]
        interval[k,] = c(y0[k]-sigmahat*sqrt(s0[k]+1)*t_alfa2, y0[k], y0[k]+sigmahat*sqrt(s0[k]+1)*t_alfa2)
      }
    } else {
      for (k in 1:K){
        dist_e_sim = gw.dist(utm_ev_sp, pcoords[k,1:2], focus=0, p=2, theta=0, longlat=F)
        dist_s_sim = gw.dist(utm_st_sp, pcoords[k,3:4], focus=0, p=2, theta=0, longlat=F)
        
        We = diag(gauss_kernel(dist_e_sim[,1],bwe))
        Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
        beta_e = Ae %*% y_tilde
        print(c("beta_e",k))
        
        y_tilde_s = (I-He)%*%y_tilde
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
        As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
        beta_s = As %*% y_tilde_s
        print(c("beta_s",k))
        
        Q = rbind(Ac, Ae%*%(I-Xc%*%Ac), As%*%(I-He)%*%(I-Xc%*%Ac))
        s0[k] = t(x0[k,])%*%Q%*%t(Q)%*%x0[k,]
        
        betas[k,] = c(beta_c, beta_e, beta_s)
        y0[k] = x0[k,]%*%betas[k,]
        t_alfa2 = qt(1-alfa/2, delta1^2/delta2)
        var0[k] = sigma2hat*s0[k]
        interval[k,] = c(y0[k]-sigmahat*sqrt(s0[k]+1)*t_alfa2, y0[k], y0[k]+sigmahat*sqrt(s0[k]+1)*t_alfa2)
      }
    }
  }
  
  
  result <- list("s0" = s0, 
                 "y0" = y0,
                 "betas" = betas,
                 "interval" = interval,
                 "var0" = var0,
                 "sigma2hat" = sigma2hat)
  return(result)
}
#ESC_only_calibration----
ESC_only_calibration = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
    print(c("He",i))
  }
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    Hs[i,] = Xs[i,] %*% As
    print(c("Hs",i))
  }
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
#ESC_only_constant_intercept_calibration----
#ESC, only constant intercept, only calibration
ESC_only_constant_intercept_calibration = function(Xe, Xs, y, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  Xc = rep(1,N)
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = 1
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
    print(c("He",i))
  }
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    Hs[i,] = Xs[i,] %*% As
    print(c("Hs",i))
  }
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  return(calibration)
}
#ESC_no_intercept_calibration----
ESC_no_intercept_calibration = function(Xc, Xe, Xs, y, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
    print(c("He",i))
  }
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    Hs[i,] = Xs[i,] %*% As
    print(c("Hs",i))
  }
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
#ESC_grid_creation----
ESC_grid_creation = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid, emgwr){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  B = emgwr$B
  He = emgwr$He
  Hs = emgwr$Hs
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_s for whole grid
  for (i in 1:L){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    beta_s[,i] = As %*% y_tilde
    print(c("beta_s",i))
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-Hs)%*%y_tilde
  
  #5)compute beta_e for whole grid
  for (i in 1:L){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    beta_e[,i] = Ae %*% y_tilde_s
    print(c("beta_e",i))
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}

#SEC_only_calibration----
SEC_only_calibration = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
    print(c("Hs",i))
    #print(c("Hs",i))
  }
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
    print(c("He",i))
    #print(c("He",i))
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
#SEC_only_constant_intercept_calibration----
#SEC, only constant intercept, only calibration
SEC_only_constant_intercept_calibration = function(Xe, Xs, y, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  Xc = rep(1,N)
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = 1
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
    print(c("Hs",i))
  }
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
    print(c("He",i))
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  return(calibration)
}
#SEC_no_intercept_calibration----
SEC_no_intercept_calibration = function(Xc, Xe, Xs, y, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
    print(c("Hs",i))
  }
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
    print(c("He",i))
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
#SEC_grid_creation----
SEC_grid_creation = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid, emgwr){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  B = emgwr$B
  He = emgwr$He
  Hs = emgwr$Hs
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_e for whole grid
  for (i in 1:L){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    beta_e[,i] = Ae %*% y_tilde
    print(c("beta_e",i))
    #print(c("beta_e",i))
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-He)%*%y_tilde
  
  #5)compute beta_s for whole grid
  for (i in 1:L){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    beta_s[,i] = As %*% y_tilde_s
    print(c("beta_s",i))
    #print(c("beta_s",i))
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}
#gauss_kernel----
#Gaussian kernel
#d  a vector of distances (one column of a distance matrix)
#h is the  bandwidth
#return a vector of weight

gauss_kernel = function(d,h){
  exp(-0.5*(d/h)^2)
}
#mixed_SC_calibration----
mixed_SC_calibration = function(Xc, Xs, y, bw, intercept, coords_s_sim){
  
  N = length(y) #y vector of responses
  
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }

  I = diag(rep(1,N))
  n_c = dim(Xc)[2] #number constant of covariates
  if (is.null(n_c)) {n_c = 1}
  n_s = dim(Xs)[2] #number of site-dependent covariates
  if (is.null(n_s)) {n_s = 1}
  Hs = matrix(0,N,N)
  dist_s_sim_cal = gw.dist(coords_s_sim, coords_s_sim, focus=0, p=2, theta=0, longlat=F)
  
  #create Hs
  if (n_s > 1){
    for (i in 1:N){
      Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      Hs[i,] = Xs[i,] %*% As
      print(c("H_s",i))
    }
  }
  
  else if (n_s == 1){
    for (i in 1:N){
      Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      Hs[i,] = Xs[i] %*% As
      print(c("H_s",i))
    }
  }
  
  betas <- list("Hs" = Hs)
  return(betas)
}
#mixed_SC_no_intercept_calibration----
mixed_SC_no_intercept_calibration = function(Xc, Xs, y, bw, coords_s_sim){
  
  N = length(y) #y vector of responses
  I = diag(rep(1,N))
  n_c = dim(Xc)[2] #number constant of covariates
  if (is.null(n_c)) {n_c = 1}
  n_s = dim(Xs)[2] #number of site-dependent covariates
  if (is.null(n_s)) {n_s = 1}
  Hs = matrix(0,N,N)
  dist_s_sim_cal = gw.dist(coords_s_sim, coords_s_sim, focus=0, p=2, theta=0, longlat=F)
  
  #create Hs
  if (n_s > 1){
    for (i in 1:N){
      Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      Hs[i,] = Xs[i,] %*% As
      print(c("H_s",i))
    }
  }
  
  else if (n_s == 1){
    for (i in 1:N){
      Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      Hs[i,] = Xs[i] %*% As
      print(c("H_s",i))
    }
  }
  
  betas <- list("Hs" = Hs)
  return(betas)
}
