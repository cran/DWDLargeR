##**********************************************************************
## DWD: an SGS based ADMM method for solving 
##  min sum_i (1/ri^q) + C*<e,xi> 
##  s.t  r = ZT*w + beta*y + xi, r > 0, xi>=0, norm(w) <= 1
##
##  INPUT:
##   X = dim x n, n = sample size, dim = feature dimension
##   y = nx1 (+1,-1 labels)
##   C = penalty parameter
##   expon = exponent q
##   tol = stopping tolerance (default:1e-5)
##   maxIter = maximum iterations allowed (default:2000)
##   method = 1 (default): sGS-ADMM
##          = 2 directly extended ADMM
##   printDetails = switch on printer for the iteration details
##  OUTPUT:
##   w,beta,xi,r = variables for primal problem
##   alpha = variable for dual problem
##   info = run information
##  --------------------------------------------------------------------
##  Copyright (c) 2016 by
##  Xin-Yee Lam, J.S. Marron, Defeng Sun,
##  Kim-Chuan Toh (corresponding author)
##**********************************************************************

genDWD = function(X,y,C,expon, tol = 1e-5, maxIter = 2000, method = 1, printDetails = 0,
                  rmzeroFea = 1, scaleFea = 1){
  
  tau = 1.618

  dim = nrow(X)
  n = ncol(X)
  idxpos = which(y>0)
  idxneg = which(y<0)
  np = length(idxpos)
  nn = length(idxneg)
  tstart = proc.time()

  nnz = sum(X@ra!=0)
  if (nnz > 0.4*dim*n && dim <= 5000){
    X = as.matrix(X)
  }
    
  ##
  ## remove zero features
  ##
  if (rmzeroFea!=0){
    normX = sqrt(rowSums(as(X*X,"dgCMatrix")))
    nzrow = which(normX>0)
    if (length(nzrow) < length(normX)){
      if (is.matrix.csr(X)){
        X = rbind(X[nzrow,1:n], 0*as.matrix.csr(1,1,n))
      }
      else{
        X = rbind(X[nzrow,1:n], 0*rep(0,n))
      }
      dim = nrow(X)
    }
  }
  
  ##
  ## scale features (to have roughly same magnitude)
  ##
  if (scaleFea!=0){
    DD = 1
    if(dim > 0.5*n){
      normX = sqrt(rowSums(as(X*X,"dgCMatrix")))
      if (max(normX) > 2*min(normX)){
        if (dim > 3*n){
          DD = new("matrix.csr", ra = 1/pmax(1,sqrt(normX)), ja = 1:dim, ia = 1:(dim+1), dimension = c(dim,dim))
        }
        else{
          DD = new("matrix.csr", ra = 1/pmax(1,normX), ja = 1:dim, ia = 1:(dim+1), dimension = c(dim,dim))
        }
        if (is.matrix.csr(X)){
          X = DD %*% X
        } else{
          X = as.matrix(DD) %*% X
        }
      }
    }
  }
  
  ##
  use_balanced = TRUE;
  if (use_balanced){
    K = n/log(n);
    tmpvec = matrix(1,n,1);
    tmpvec[idxpos] = matrix(nn/K,np,1)
    tmpvec[idxneg] = matrix(np/K,nn,1)
    weightoptions = 2
    if(weightoptions == 0){
      resweight = 1
      penweight = 1
    }
    else if(weightoptions == 1){ 
      resweight = tmpvec^(1/(2+expon))
      penweight = tmpvec^(1/(2+expon))
    }
    else if(weightoptions == 2){ 
      resweight = tmpvec^(1/(1+expon))
      penweight = 1
    }
    else if(weightoptions == 3){ 
      resweight = 1
      penweight = tmpvec^(1/(2+expon))
    }
    resweight = resweight/max(resweight)
    y = y/resweight
    Cvec = C*(penweight*resweight)
  }
  else{
    Cvec = C*matrix(1,n,1)
  }
  maxC = max(Cvec);
  
  ##
  if(is.matrix.csr(X)){
    Z = X %*% new("matrix.csr", ra = as.numeric(y), ja = 1:n, ia = 1:(n+1), dimension = c(n,n))
  }else{
    Z = X %*% diag(as.vector(y))
  }
  scale_data = 1
  if (scale_data==1){
    Zscale = sqrt(fnorm(X))
    Z = Z/Zscale
    sigma = min(10*C,1*n) ## important
  } else {
    Zscale = 1
    sigma = max(1,log(n/dim)*fnorm(X)) 
  }
  sigma = sigma^expon
  sigmastart = sigma
  normZ = 1+sqrt(max(colSums(as(Z*Z,"dgCMatrix"))))

  ##  
  ## initial iterate
  ##
  r = matrix(1,n,1) 
  wbeta = matrix(0,dim+1,1)
  u = matrix(0,dim,1) 
  xi = matrix(0,n,1)
  alpha = matrix(0,n,1)
  p = matrix(0,dim,1)
  
  ##  
  ## chol of M
  ##  
  const = 1
  if (dim>5000){
    if (n<0.2*dim && n<=2500){
      Solver = 'SMW'
    } else{
      Solver = 'iterative'
    }
  } else{
    Solver = 'direct'
  }
  
  ZT = t(Z)
  if (Solver == 'direct'){ 
    if (is.matrix.csr(X)){
      M1 = Z %*% ZT + const*new("matrix.csr", ra = rep(1,dim), ja = 1:dim, ia = 1:(dim+1), dimension = c(dim,dim))
      M4 = as.matrix.csr(t(y) %*% y)
    } else{
      M1 = Z %*% ZT + const*diag(dim)
      M4 = t(y) %*% y
    }
    M2 = Z %*% y
    M3 = t(M2)
    M = rbind(cbind(M1,M2),cbind(M3,M4))
    if (dim > 4000){
      R = chol(M,tmpmax = 1000*dim)
    }else{
      R = chol(M)
    }
  }
  else if (Solver == 'SMW'){
    normy = fnorm(y)
    yunit = y/normy    
    H11 = new("matrix.csr", ra = rep(1,n), ja = 1:n, ia = 1:(n+1), dimension = c(n,n)) + (1/const)*(ZT %*% Z)
    R = chol(H11)
    invH11yunit = linsysolve(R,yunit)
    schurmat = t(yunit) %*% invH11yunit
    schurvec = c(invH11yunit, -1)
  }
  else if (Solver == 'iterative'){
    ff = list("Z"=Z,"ZT"=ZT,"y"=y,"const"=const)
    diagM = rowSums(as(Z*Z,"dgCMatrix")) + const
    L = list()
    L[["invdiagM"]] = 1/c(diagM,1)
    L[["precond"]] = 1
  }
  if (printDetails!=0){
    cat('\n------------------------------------------------------')
    cat('--------------------------------------------------------')
    cat('\n iter  time    sigma    primfeas  dualfeas')
    cat(' relcomp   primobj    dualobj    relgap'); 
    cat('   errNewton  min-r   dc_meas  dcmpute   psqmriter   trainerr')
    cat('\n------------------------------------------------------')
    cat('---------------------------------------------------------')
    cat(sprintf('\n sample size = %3.0f, feature dimension = %3.0f',n,dim))
    cat(sprintf('\n expon = %2.1f, penalty constant C = %4.2e',expon,C))
    cat(sprintf('\n initial sigma = %3.2e',sigma))
    cat(sprintf('\n Zscale = %3.2e',Zscale))
    cat(sprintf('\n norm(scaled Z) = %3.2e',fnorm(Z)))
    cat(sprintf('\n Linear system solver = %s',Solver))
    cat('\n------------------------------------------------------')
    cat('--------------------------------------------------------')
  }
  ##
  ## Set up run history  
  ##
  primfeas_his = c(); dualfeas_his = c()
  relcomp_his  = c(); trainerr_his = c()
  primobj_his  = c(); dualobj_his  = c()
  relgap_his   = c(); psqmr_his    = c()
  Newton_his   = c(); doublecompute_his = c()
  
  ##
  ## Main Loop
  ##  
  breakyes = FALSE
  doublecompute = 0
  psqmriter = 0
  rhs = matrix(0,dim+1,1)
  rhsnew = matrix(0,dim+1,1)

  for (iter in 1:maxIter){
    rold = r; wbetaold = wbeta; uold = u
    xiold = xi; alphaold = alpha; pold = p
    
    # update w, beta
    tmp = rold - xiold + alphaold/sigma
    rhs1 = as.matrix(Z %*% tmp) + const*uold + pold/sigma
    rhs = rbind(rhs1,t(y) %*% tmp)
    if (Solver == 'iterative'){
      psqmrTol = max(min(5e-2,1/(iter^2)),1e-8)*max(1,sqrt(fnorm(rhs)))
      psqmrMaxiter = 100
      runpsqmr = psqmr(ff,rhs,L,wbetaold,psqmrTol,psqmrMaxiter)
      wbeta = runpsqmr$x
      resnrm = runpsqmr$resnrm
      solve_ok = runpsqmr$solve_ok
      psqmriter = psqmriter + length(resnrm) - 1
      if (solve_ok != 1 && printDetails!=0){
        cat(sprintf('\n iter=%2.0f: PSQMR not successful,num iter=%3.0f,residual=%3.2e',iter,length(resnrm)-1,resnrm[length(resnrm)]/fnorm(rhs)))
      }
    }
    else if (Solver == 'direct'){
      wbeta = linsysolve(R,rhs)
    }
    else if (Solver == 'SMW'){
      wbeta = smw(R,Z,ZT,yunit,schurmat,schurvec,normy,const,rhs)
    }

    # update r
    w = wbeta[1:dim]
    beta = wbeta[dim+1]
    ZTwpbetay = as.matrix(ZT %*% w) + beta*y
    cc = ZTwpbetay + xiold - alphaold/sigma;
    runNewton = polyRootsNewton(cc,expon,sigma,rold)
    r = runNewton[1:n]
    errNewton = runNewton[n+1]
    iterNewton = runNewton[n+2]
    r = pmax(r,0)

    # update w, beta again
    if(method==1){
      doublecompute_measure = normZ*fnorm(r-rold)*iter^1.5
      if ((doublecompute_measure > 10) || (iter < 50)){
        doublecompute = doublecompute + 1
        tmpnew = r - xiold + alphaold/sigma 
        rhsnew1 = as.matrix(Z %*% tmpnew) + const*uold + pold/sigma
        rhsnew = rbind(rhsnew1,t(y) %*% tmpnew)
        if (Solver == 'iterative'){
          runpsqmr = psqmr(ff,rhsnew,L,wbeta,psqmrTol,psqmrMaxiter)
          wbeta = runpsqmr$x
          resnrm = runpsqmr$resnrm
          solve_ok = runpsqmr$solve_ok
          psqmriter = psqmriter + length(resnrm) - 1
          if (solve_ok != 1 && printDetails!=0){
            cat(sprintf('\n iter=%2.0f: PSQMR not successful,num iter=%3.0f,residual=%3.2e',iter,length(resnrm)-1,resnrm[length(resnrm)]/fnorm(rhs)))
          }
        }
        else if (Solver == 'direct'){
          wbeta = linsysolve(R,rhsnew)
        }
        else if (Solver == 'SMW'){
          wbeta = smw(R,Z,ZT,yunit,schurmat,schurvec,normy,const,rhsnew)
        }
        w = wbeta[1:dim]
        beta = wbeta[dim+1]
        ZTwpbetay = as.matrix(ZT %*% w) + beta*y
      }
    }
    else{
      doublecompute_measure = 0;
      doublecompute = 0;
    }
    
    # update u, xi
    uinput = w -pold/(const*sigma);
    u = Zscale*uinput/max(Zscale,fnorm(uinput));     
    xiinput = r - ZTwpbetay + (alphaold-Cvec)/sigma;
    xi = pmax(as.vector(xiinput),0);  

    # update alpha, p
    Rp = ZTwpbetay +xi-r; 
    alpha = alphaold -tau*sigma*Rp;
    p = pold -(tau*sigma*const)*(w-u);
    
    # check for termination
    rexpon1 = r^(expon+1); 
    comp1 = abs(t(y) %*% alpha)
    comp2 = abs(t(xi) %*% (Cvec-alpha));      
    comp3 = min(fnorm(alpha*rexpon1-expon),fnorm(alpha-expon/rexpon1)^2); 
    relcomp = max(comp1,comp2,comp3)/(1+maxC); 
    primfeas = max(fnorm(Rp),fnorm(w-u),max(fnorm(w)-Zscale,0))/(1+maxC);      
    dualfeas = max(fnorm(pmin(0,alpha)),fnorm(pmax(alpha-Cvec,0)))/(1+maxC);   
    trainerr = length(which(ZTwpbetay<=0))/n*100; 
    if (((max(primfeas,dualfeas) < tol) && (iter%%20==1)) || (iter%%100==1)){
      primobj1 = sum(r/rexpon1)
      primobj2 = sum(Cvec*xi)+1e-8
      primobj = primobj1 + primobj2
      kappa = ((expon+1)/expon)*expon^(1/(expon+1));          
      dualobj = kappa*sum(pmax(0,alpha)^(expon/(expon+1)))-Zscale*fnorm(as.matrix(Z %*% alpha))
      relgap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj))
    }
    tol2 = 0.1/2
    if ((iter > 50) && (max(primfeas,dualfeas) < tol) && (min(relcomp,relgap) < sqrt(tol)) && (((relcomp < tol2) && (relgap < sqrt(tol))) || (((relcomp < sqrt(tol)) && (relgap < tol2))))){
      KKTerr = max(max(primfeas,dualfeas),min(relcomp,relgap))
      breakyes = 1
      cat(sprintf('\n Algorithm stopped with error %3.2e',KKTerr))
    }
    if ((iter > 50) && (max(primfeas,dualfeas) < 5*tol) && (min(relcomp,relgap) < 10*tol) && (((relcomp < tol2) && (relgap < sqrt(tol))) || (((relcomp < sqrt(tol)) && (relgap < tol2))))){
      KKTerr = max(max(primfeas,dualfeas),min(relcomp,relgap))
      breakyes = 2
      cat(sprintf('\n Algorithm stopped with error %3.2e',KKTerr))
    }
    if ((iter > 50) && (max(primfeas,dualfeas) < tol) && (fnorm(alpha)/(1+maxC) < 1e-3)){
      KKTerr = max(max(primfeas,dualfeas),min(relcomp,relgap))
      breakyes = 3
      cat(sprintf('\n Algorithm stopped with error %3.2e',KKTerr))
    }
    if (iter <= 100){
      print_iter = 20
    }else{
      print_iter = 100
    }     
    if ((iter%%print_iter==1) || (breakyes>0)){        
      ttime = as.numeric((proc.time()-tstart)[3])
      if (printDetails!=0){
        cat(sprintf('\n%4.0f| %6.2f| %3.2e| %3.2e %3.2e %3.2e| %5.4e %5.4e %3.2e| %3.2e %3.2e| %3.2e %3.0f| %5.0f %4.2f|',iter,ttime,sigma,primfeas,dualfeas,relcomp,primobj,dualobj,relgap,errNewton,min(r),doublecompute_measure,doublecompute,psqmriter,trainerr));
        cat(sprintf(' %3.2e',fnorm(alpha)/(1+maxC)))
      }
    }
    primfeas_his = c(primfeas_his,primfeas)
    dualfeas_his = c(dualfeas_his,dualfeas)
    relcomp_his = c(relcomp_his,relcomp)
    trainerr_his = c(trainerr_his,trainerr)
    primobj_his = c(primobj_his,primobj)
    dualobj_his = c(dualobj_his,dualobj)
    relgap_his = c(relgap_his,relgap)
    psqmr_his = c(psqmr_his,psqmr)
    Newton_his = c(Newton_his,iterNewton)
    doublecompute_his = c(doublecompute_his,doublecompute)
    
    # adjust sigma
    sigma_update_iter = sigma_update(iter);
    if (iter%%sigma_update_iter==0){
      primfeas2 = max(primfeas,0.2*tol);
      dualfeas2 = max(dualfeas,0.2*tol);
      ratio = primfeas2/dualfeas2;
      const2 = 1.1; 
      if (max(ratio,1/ratio) > 500){
        const2 = const2*2;
      }
      else if (max(ratio,1/ratio) > 50){
        const2 = const2*1.5;
      }
      if (ratio > 5){
        sigma = min(sigma*const2,1e6);
      }
      else if (1/ratio > 5){
        sigma = max(sigma/const2,1e-3);
      }
    }
    if (breakyes){
      cat('\n')
      break;
    } 
  }
  
  ##
  ## End of main loop
  ##
  
  Zalpha = as.matrix(Z %*% alpha)
  w = w/Zscale;
  
  # Calculate train error
  res = t(X) %*% w + beta*as.matrix.csr(1,n,1)
  error = length(which(y*sign(res@ra)<=0))/n*100
  
  cat(sprintf('\n sample size = %3.0f, feature dimension = %3.0f',n,dim));
  cat(sprintf('\n positve sample = %3.0f, negative sample = %3.0f',np,nn));
  cat(sprintf('\n number of iterations = %3.0f',iter));
  cat(sprintf('\n time taken = %3.2f',ttime));
  cat(sprintf('\n error of classification (training) = %3.2f (%%)',error))
  cat(sprintf('\n primfeas = %3.2e',primfeas));
  cat(sprintf('\n dualfeas = %3.2e',dualfeas));
  cat(sprintf('\n relative gap = %3.2e\n',relgap));
  
  runhist = list("primfeas"=primfeas_his, "dualfeas"=dualfeas_his, "relcomp" = relcomp_his, "trainerr"=trainerr_his, "primobj"=primobj_his, "dualobj"=dualobj_his, "relgap"=relgap_his, "Newton"=Newton_his, "doublecompute"=doublecompute_his)  
  info = list("iter"=iter,"time"=ttime,"penaltyParameter"=C,"sampsize"=n,"np" = np, "nn"=nn, "dimension"=dim, "sigmastart"=sigmastart, "primfeas"=primfeas, "dualfeas"=dualfeas, "relcomp"=relcomp, "relgap"=relgap, "primobj"=primobj, "dualobj"=dualobj, "trainerr"=trainerr, "psqmr"=psqmriter,"doublecompute"=doublecompute);
  return(list("w"=w,"beta"=beta,"xi"=xi,"r"=r,"alpha"=alpha,"info"=info,"runhist"=runhist))
}

##**********************************************************************
##**********************************************************************
sigma_update = function(iter){
  const = 0.5
  if (iter <= 25){
    sigma_update_iter = 10*const
  }
  else if (iter <= 50){
    sigma_update_iter = 20*const
  }
  else if (iter <= 100){
    sigma_update_iter = 40*const
  }
  else if (iter <= 500){
    sigma_update_iter = 60*const
  }
  else if (iter <= 1000){
    sigma_update_iter = 80*const
  }
  else{
    sigma_update_iter = 100
  }
  return(sigma_update_iter)
}

linsysolve = function(R,r){
  if(is.matrix(R)){
    x = backsolve(R,forwardsolve(t(R),r))
  }else{
    x = backsolve(R,r)
  }
  
  return(x)
}

##**********************************************************************
##**********************************************************************
smw = function(R,Z,ZT,yunit,schurmat,schurvec,normy,const,r){
  n = length(yunit)
  dim = length(r) - 1
  b1 = (1/const)*as.matrix(ZT %*% r[1:dim]) + (r[dim+1]/normy)*yunit
  b2 = r[dim+1]/normy
  b = rbind(b1,b2)
  
  ##
  tmpvec = (t(schurvec) %*% b/schurmat) %*% schurvec
  a1 = backsolve(R,b1) - tmpvec[1:n]
  a2 = -b2 - tmpvec[n+1]
  q = matrix(0,dim+1,1)
  q[1:dim] = (1/const)*(as.matrix(Z %*% (-a1)) + r[1:dim])
  q[dim+1] = (1/normy)*(r[dim+1]/normy - t(yunit)%*%a1 - a2)
  
  return (q)
}
##**********************************************************************
##**********************************************************************
fnorm = function(x){
    if (typeof(x) == "S4"){
      xentry = x@ra
      return (sqrt(sum(xentry * xentry)))
    }else{
      return (sqrt(sum(x*x)))
    }

}
##**********************************************************************
polyRootsNewton = function(c,q,sigma,x0){
  tol = 1e-12
  x = x0
  d = q/sigma; cq1 = c*(q+1); q2 = q+2; 
  maxiter = 50
  
  if (q == 1){
    for (iter in 1:maxiter){
      idx = which(x<=0);
      if (length(idx)!=0){
        x[idx]=max(0,0.5+c[idx])
      }
      xq1 = x*x
      xq = x
      grad = xq1*(x-c)-d
      hess = xq*(q2*x - cq1)
      x = x-grad/hess
      err = max(abs(grad))
      if (err<tol){
        break
      }
    }
  }
  else if (q == 2){
    for (iter in 1:maxiter){
      idx = which(x<=0);
      if (length(idx)!=0){
        x[idx]=max(0,0.5+c[idx])
      }
      xq = x*x
      xq1 = xq*x
      grad = xq1*(x-c)-d
      hess = xq*(q2*x-cq1)
      x = x - grad/hess
      err = max(abs(grad))
      if (err<tol){
        break
      }
    }
  }
  else if (q == 3){
    for (iter in 1:maxiter){
      idx = which(x<=0);
      if (length(idx)!=0){
        x[idx]=max(0,0.5+c[idx])
      }
      x2 = x*x
      xq = x2*x
      xq1 = x2*x2
      grad = xq1*(x-c)-d
      hess = xq*(q2*x-cq1)
      x = x - grad/hess
      err = max(abs(grad))
      if (err<tol){
        break
      }
    }
  }
  else if (q == 4){
    for (iter in 1:maxiter){
      idx = which(x<=0);
      if (length(idx)!=0){
        x[idx]=max(0,0.5+c[idx])
      }
      x2 = x*x
      xq = x2*x2
      xq1 = xq*x
      grad = xq1*(x-c)-d
      hess = xq*(q2*x-cq1)
      x = x - grad/hess
      err = max(abs(grad))
      if (err<tol){
        break
      }
    }
  }
  
  return(rbind(x,err,iter))
}

##**********************************************************************
##**********************************************************************

psqmr = function(ff,b,L,x0,tol,maxit){
  
  N = length(b)
  if (!exists('maxit')) maxit = max(1000,sqrt(length(b)))
  if (!exists('tol')) tol = 1e-8*fnorm(b)
  if (!exists('L')) L[["precond"]] = 0
  if (!exists('x0')) x0 = matrix(0,N,1)
  
  solve_ok = 1
  stagnate_check = 20*1
  miniter = 0
  printlevel = FALSE
  
  ##
  x = x0
  if (fnorm(x) > 0){
    Aq = vecMultiply(ff,x)
  }else{
    Aq = matrix(0,N,1)
  }
  r = b - Aq
  err = fnorm(r)
  resnrm = err
  minres = err
  
  ##
  q = precondfun(L,r)
  tau_old = fnorm(q)
  rho_old = as.numeric(t(r) %*% q)
  theta_old = 0
  d = matrix(0,N,1)
  res = r
  Ad = matrix(0,N,1)
  
  ##
  ## Main Loop
  ##
  tiny = 1e-30
  for (iter in 1:maxit){
    Aq = vecMultiply(ff,q)
    sigma = as.numeric(t(q) %*% Aq)
    if (abs(sigma) < tiny){
      solve_ok = 2
      if (printlevel) cat('s1')
      break
    }else{
      alpha = rho_old/sigma
      r = r - alpha*Aq
    }
    u = precondfun(L,r)
    
    ##
    theta = fnorm(u)/tau_old
    c = 1/sqrt(1+theta^2)
    tau = tau_old*theta*c
    gam = (c^2*theta_old^2)
    eta = c^2*alpha
    d = gam*d + eta*q
    x = x + d
    ##-------------stopping conditions---------
    Ad = gam*Ad + eta*Aq
    res = res - Ad
    err = fnorm(res)
    resnrm = c(resnrm,err)
    if (err < minres) minres = err
    if ((err < tol) && (iter > miniter) && (t(b) %*% x > 0)) break
    if ((iter > stagnate_check) && (iter > 10)){
      ratio = resnrm[(iter-9):(iter+1)]/resnrm[(iter-10):iter]
      if ((min(ratio) > 0.997) && (max(ratio) < 1.003)){
        if (printlevel) cat('s')
        solve_ok = -1
        break
      }
    }
    ##--------------------------------------------
    if (abs(rho_old) < tiny){
      solve_ok = 2
      cat('s2')
      break
    }else{
      rho = as.numeric(t(r) %*% u)
      beta = rho/rho_old
      q = u + beta*q
    }
    rho_old = rho
    tau_old = tau
    theta_old = theta
  }
  if (iter == maxit) solve_ok = -2
  if (solve_ok != -1){
    if (printlevel) cat(' ')
  }
  
  return (list("x"=x,"resnrm"=resnrm,"solve_ok"=solve_ok))
}

precondfun = function(L,r){
  
  precond = L$precond
  
  if (precond == 0){
    q = r
  }
  else if (precond == 1){
    q = L$invdiagM * r
  }
  else if (precond == 2){
    q = backsolve(L$R,r)
  }
  
  return(q)
}

vecMultiply = function(ff,x){
  
  Z = ff$Z
  ZT = ff$ZT
  y = ff$y
  const = ff$const
  d = length(x) - 1
  w = x[1:d]; beta = x[d+1]
  tmp = as.matrix(ZT %*% w) + beta*y
  Aq = matrix(0,d+1,1)
  Aq[1:d] = as.matrix(Z %*% tmp) + const*w
  Aq[d+1] = t(y) %*% tmp
  
  return (Aq)
}


