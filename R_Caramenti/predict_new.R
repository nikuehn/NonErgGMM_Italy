#emgwr_prediction_points_no_param----
#prediction for points with their coordinates without already computed delta1, delta2 and H_hat
emgwr_prediction_points_no_param2 = function(Xc, Xe, Xs, y, emgwr, method, bwe, bws, utm_ev_sp, utm_st_sp, pc, pe, ps, pcoords, alfa){
  
  N = length(y)
  n_sample <- N
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
    dist_e_sim = gw.dist(utm_ev_sp, coordinates(t(as.matrix(pcoords)[1:2])), focus=0, p=2, theta=0, longlat=F)
    dist_s_sim = gw.dist(utm_st_sp, coodinates(t(as.matrix(pcoords)[3:4])), focus=0, p=2, theta=0, longlat=F)
    
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
        #print(c("beta_s",k))
        
        y_tilde_s = (I-Hs)%*%y_tilde
        
        We = diag(gauss_kernel(dist_e_sim[,1],bwe))
        Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
        beta_e = Ae %*% y_tilde_s
        #print(c("beta_e",k))
        
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
        #print(c("beta_e",k))
        
        y_tilde_s = (I-He)%*%y_tilde
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bws))
        As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
        beta_s = As %*% y_tilde_s
        #print(c("beta_s",k))
        
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
