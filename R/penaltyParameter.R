##**************************************************************
## find a good choice of C by the data X
## each column of X is a sample
##**************************************************************

penaltyParameter = function(X,y,expon,rmzeroFea = 1, scaleFea = 1){
  set.seed(0)
  
  dim = nrow(X)
  sampsize = ncol(X)
  
  ##
  ## remove zero features
  ##
  if (rmzeroFea!=0){
    normX = sqrt(rowSums(as(X*X,"dgCMatrix")))
    nzrow = which(normX>0)
    if (length(nzrow) < length(normX)){
      X = rbind(X[nzrow,1:sampsize], 0*as.matrix.csr(1,1,sampsize))
      dim = nrow(X)
    }
  }
  
  ##
  ## scale features (to have roughly same magnitude)
  ##
  if (scaleFea!=0){
    DD = 1
    if(dim > 0.5*sampsize){
      normX = sqrt(rowSums(as(X*X,"dgCMatrix")))
      cat(sprintf('\n max-normX, min-normX = %3.2e, %3.2e',max(normX),min(normX)))
      if (max(normX) > 2*min(normX)){
        if (dim > 3*sampsize){
          DD = new("matrix.csr", ra = 1/pmax(1,sqrt(normX)), ja = 1:dim, ia = 1:(dim+1), dimension = c(dim,dim))
        }
        else{
          DD = new("matrix.csr", ra = 1/pmax(1,normX), ja = 1:dim, ia = 1:(dim+1), dimension = c(dim,dim))
        }
        X = DD %*% X
      }
    }
  }
  
  positive=which(y==1)
  negative=which(y==-1)
  if ((dim > 1e4) && (sampsize > 1e4)){
    len = 100
  }
  else{
    len = 200
  }
  if (length(positive) > len){
    idx = sample(1:length(positive))
    positive = positive[idx[1:len]]
  }
  if (length(negative) > len){
    idx = sample(1:length(negative))
    negative = negative[idx[1:len]]
  }   
  ##    
  posX = as.matrix(X[1:dim,positive])
  n1 = ncol(posX)
  d1 = nrow(posX)
  negX = as.matrix(X[1:dim,negative])
  n2 = ncol(negX)
  d2 = nrow(negX)
  ddist = matrix(0,n2,n1)

  for (i in 1:n2){
    for (j in 1:n1){
      Xtmp = posX[1:d1,j] - negX[1:d2,i]
      ddist[i,j] = sqrt(sum(Xtmp*Xtmp))
    }
  }
  
  ddist = as.vector(ddist)   
  dd = median(ddist)
  if (expon==1){     
    const = log(sampsize)*max(1000,dim)^(1/3)
  }
  else{
    const = 10*log(sampsize)*max(1000,dim)^(1/3)      
  }
  C = 10^(expon+1)*max(1,const/dd^(expon+1))
  return(C)
}
