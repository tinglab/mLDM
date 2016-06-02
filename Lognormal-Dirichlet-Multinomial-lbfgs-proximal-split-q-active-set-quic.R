# compute the det via cholesky decomposition
computeLogDet <- function(X){
  eigen_value <- eigen(X, symmetric=TRUE, only.values=TRUE)$values

  return(sum(log(eigen_value)))
}

# check if alpha will overflow and save in log_file
checkAlpha <- function(alpha_i, mu_i, z_i){
  if(length(which(alpha_i <= 0)) > 0){
    print('error alpha_i:')
    print(alpha_i[alpha_i <= 0])
    print('error mu_i:')
    print(mu_i[alpha_i <= 0])
    print('error z_i:')
    print(z_i[alpha_i <= 0])
  }
}

# when alpha i is too small, set to a not zero threshold
filterAlphaZero <- function(x){
    min_integer <- -30
    x[x<exp(min_integer)] <- exp(min_integer)
    return(x)
}

# compute alpha_i = exp(u_i + z_i)
getAlpha <- function(mu_i, z_i){
  alpha_i <- exp(mu_i + z_i)
  checkAlpha(alpha_i, mu_i, z_i)
  #alpha_i <- filterAlphaZero(alpha_i)
  return(alpha_i)
}

# compute LnT
LnT <- function(x) {
  return(sum(lgamma(x)) - lgamma(sum(x)))
}

# compute LnTT test
LnTT <- function(x) {
  print('lg sum x:')
  print(lgamma(sum(x)))
  return(sum(lgamma(x)) - lgamma(sum(x)))
}

# compute Ln'T
DerLnT <- function(x) {
  return(digamma(x) - digamma(sum(x)))
}

# compute lgamma part of objective funciton for Z_i
computeAlphaVector <- function(alpha_ix, alpha_i){
  obj <- -sum(lgamma(alpha_ix) - lgamma(alpha_i)) + (lgamma(sum(alpha_ix)) - lgamma(sum(alpha_i)))
  return(obj)
}

# compute the objective function for the z_i
objZi <- function(xv, mu_i, x_i, Theta, n, B0) {
  z_i <- xv
  alpha_i <- getAlpha(mu_i, z_i)
  alpha_ix <- alpha_i + x_i
  obj <- computeAlphaVector(alpha_ix, alpha_i)

  obj <- obj + t(z_i - B0)%*%Theta%*%(z_i - B0)/2
  return(as.numeric(obj/n))
}

# compute digamma part of objective funciton for gradient of Z_i
computeAlphaDerVector <- function(alpha_ix, alpha_i){
  der <- -(digamma(alpha_ix) - digamma(alpha_i)) + (digamma(sum(alpha_ix)) - digamma(sum(alpha_i)))
  return(der)
}

# compute the derivate of the objective function for the  Z_i
derObjZi <- function(xv, mu_i, x_i, Theta, n, B0) {
  z_i <- xv
  alpha_i <- getAlpha(mu_i, z_i)
  alpha_ix <- alpha_i + x_i

  der <- computeAlphaDerVector(alpha_ix, alpha_i)*alpha_i + Theta%*%(z_i - B0)

  return(as.numeric(der/n))
}

# compute the second derivate of the objective funciton for the Z_i
derObjZi2 <- function(z_i, mu_i, x_i, Theta, n, B0){
  alpha_i <- getAlpha(mu_i, z_i)
  p <- length(alpha_i)
  alpha_ix <- alpha_i + x_i
  HessZi <- Theta
  cons1 <- digamma(sum(alpha_ix)) - digamma(sum(alpha_i))
  cons2 <- trigamma(sum(alpha_ix)) - trigamma(sum(alpha_i))
  HessZi <- HessZi + matrix(alpha_i,p,1)%*%matrix(alpha_i,1,p)*cons2

  der1 <- digamma(alpha_ix) - digamma(alpha_i)
  der2 <- trigamma(alpha_ix) - trigamma(alpha_i)
  HessZi <- HessZi + diag(-der1*alpha_i) + diag(-der2*alpha_i^2)

  HessZi <- HessZi + diag(alpha_i)*cons1

  return(HessZi/n)
}

# compute lgamma matrix part of objective funciton for matrix B
computeAlhpaMatrix <- function(AlphaMX, AlphaM){
  obj <- -sum(lgamma(AlphaMX) - lgamma(AlphaM)) + sum(lgamma(colSums(AlphaMX)) - lgamma(colSums(AlphaM)))
  return(obj)
}

# compute  the objective function for the matrix B
objB <- function(Bv, ZB0, X, M, lambda2, q, p, n){
  B <- matrix(Bv, q, p)
  AlphaM <- exp(t(M%*%B) + ZB0)
  AlphaMX <- AlphaM + t(X)

  obj <- computeAlhpaMatrix(AlphaMX, AlphaM)

  return(obj/n)
}

# compute digamma matrix part of objective funciton for the matrix B
computeAlphaDerMatrix <- function(AlphaMX, AlphaM){
  der <- -t(digamma(AlphaMX) - digamma(AlphaM)) + (digamma(colSums(AlphaMX)) - digamma(colSums(AlphaM)))
  der <- t(der)
  return(der)
}

# compute the derivate of the objective function for the  matrix B
derObjB <- function(Bv, ZB0, X, M, lambda2, q, p, n){
  B <- matrix(Bv, q, p)
  AlphaM <- exp(t(M%*%B) + ZB0)
  AlphaMX <- AlphaM + t(X)
  der <- (computeAlphaDerMatrix(AlphaMX, AlphaM)*AlphaM)%*%M
  der <- t(der)/n

  return(rep(der))
}

# compute the diag of the second derivate of the obj funciton for the B
derObjB2D <- function(w, ZB0, X, M, lambda2, q, p, n){
  B <- matrix(w, q, p)

  AlphaM <- exp(t(M%*%B) + ZB0)
  AlphaMX <- AlphaM + t(X)

  der1 <- -t(digamma(AlphaMX) - digamma(AlphaM)) + (digamma(colSums(AlphaMX)) - digamma(colSums(AlphaM)))
  der1 <- t(der1)

  der2 <- -t(trigamma(AlphaMX) - trigamma(AlphaM)) + (trigamma(colSums(AlphaMX)) - trigamma(colSums(AlphaM)))
  der2 <- t(der2)

  der <- (der1*AlphaM + der2*AlphaM^2)%*%M^2
  der <- t(der)/n

  return(rep(der))
}

# compute diag of Bt
computeDiag <- function(gammat, Q, Qh, d){
  ret <- rep(0,d)
  for(i in 1:d){
    ret[i] <- Q[i,]%*%Qh[,i]
  }

  return(gammat - ret)
}

# compute Bt for proximal qusi newton method
computeBtQ <- function(gammat, St, Yt, approx_num, approx_count, approx_index){
  d <- nrow(St)
  Q <- matrix(0,d,2)
  Qh <- matrix(0,2,d)
  BtQ <- list()

  if(approx_count != 0){
    real_index <- rep(0, approx_count)
    for(i in 1:approx_count){
      delta <- approx_count - i
      real_index[i] <- (approx_index - delta + approx_num - 1)%%approx_num + 1
    }
    
    Stp <- St[,real_index]
    Ytp <- Yt[,real_index]

    Q <- cbind(gammat*Stp, Ytp)

    SY <- t(Stp)%*%Ytp

    if(approx_count!=1){
      Dt <- diag(diag(SY))
    }else{
      Dt <- diag(SY)
    }
    Dt <- as.matrix(Dt)
    BtQ$Dt <- Dt

    Lt <- SY
    Lt[upper.tri(Lt, diag=TRUE)] <- 0

    Rt <- rbind(cbind(gammat*t(Stp)%*%Stp, Lt), cbind(t(Lt),-Dt))
    Rt <- solve(Rt)

    Qh <- Rt%*%t(Q)
  }

  BtQ$Q <- Q
  BtQ$Qh <- Qh

  return(BtQ)
}

# check and guarantee Bt positive definite
checkBtPD <- function(gammat, BtQ, round, d){
  if(round > 1){
    Dt <- diag(BtQ$Dt)
    if(min(Dt) < 0){
      Q <- BtQ$Q
      Qh <- BtQ$Qh
      QQh <- Q%*%Qh

      while(TRUE){
        temp <- diag(gammat*rep(1,d)) - QQh

        if(det(temp) > 1 && min(diag(temp)) > 0){
          break
        }else if(gammat > 0){
          gammat <- gammat*2
        }else{
          gammat <- 1
        }
      }

    }
  }

  return(gammat)
}

softthrehold <- function(a, b){
  return(sign(a)*max(abs(a) - b, 0))
}

# set near zero value to zero~
filter_dir <- function(direction){
  #threshold <- 1e-6
  threshold <- 1e-10
  direction[abs(direction) < threshold] <- 0
  return(direction)
}

scale_dir <- function(direction){
  direction <- direction / sqrt(sum(direction*direction))
  return(direction)
}

# select the not zero elements in w and gradient of w
# as the active set
getActiveSet <- function(w, gt, lambda, d){
  active <- list()
  threshold <- 1e-6
  act <- c()
  count <- 0

  for(i in 1:d){
    a <- w[i]
    b <- gt[i]
    add_flag <- TRUE
    if(abs(a) < threshold){
      subgt <- max(0, b - lambda, -(b + lambda))
      if(abs(subgt) < threshold){
        add_flag <- FALSE
      }
    }

    if(add_flag){
      count <- count + 1
      act[count] <- i
    }
  }

  active$act <- act
  active$num <- count

  return(active)
}

# get the direciton via coordinate descent
coordinateDescentD <- function(w, gammat, Q, Qh, gt, lambda, max_iteration, threshold){
  d <- nrow(Q)
  wt <- w
  activeSet <- getActiveSet(w, gt, lambda, d)
  act <- activeSet$act
  act_size <- activeSet$num
  #act <- c(1:d)
  #act_size <- d
  #print('active set size:')
  #print(act_size)

  # set the number of iteration to the half of size of active set
  #max_iteration <- round(activeSet$num*0.5)
  round <- 0
  Bt <- diag(gammat*rep(1,d)) - Q%*%Qh
  delta <- 1

  while(delta > threshold && round < max_iteration){
    round <- round + 1
    #print('round coor:')
    #print(round)

    wt_old <- wt
    for(index in 1:act_size){
      i <- act[index]
      a <- Bt[i,i]
      b <- gt[i] + Bt[i,]%*%(wt-w) - a*wt[i]

      wt[i] <- softthrehold(-b, lambda)/a
    }

    delta <- sum(abs(wt - wt_old)) / d

    if(is.infinite(delta) || is.na(delta)){
        print('coordinate error!')
        print(delta)
        print('alpha:')
        print(alpha)
    }
    #print('delta:')
    #print(delta)
  }

  direction <- wt - w
 # direction <- filter_dir(direction)
  direction <- scale_dir(direction)
 # print('direction scale not zero number:')
 # print(sum(abs(direction)>1e-10))

  return(direction)
}

L1_norm <- function(x){
  return(sum(abs(x)))
}

# compute objective funciton for index-th row of the matrix B
objbp <- function(xv, index, BZ, X, M, q, p, n){
  bv <- xv
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  #AlphaM <- filterAlphaZero(AlphaM)
  AlphaMX <- AlphaM + t(X)
  obj <- computeAlhpaMatrix(AlphaMX, AlphaM)

  return(obj/n)
}

# compute derivate function for index-th row of the matrix B
derObjbp <- function(xv, index, BZ, X, M, q, p, n){
  bv <- xv
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  #AlphaM <- filterAlphaZero(AlphaM)
  AlphaMX <- AlphaM + t(X)
  der <- (computeAlphaDerMatrix(AlphaMX, AlphaM)*AlphaM)%*%M[,index]
  der <- der / n

  return(as.vector(der))
}

# find new w via linesearch based on strong wolfe condition
linesearch <- function(objFunc, derObjFunc, w, direction, lambda, max_linesearch, f0, g0, delta1, delta2, ...){
  beta <- 0.5
  alpha <- 1
  k <- 0
  exist <- FALSE
  wt <- w
  ret <- list()
  f1 <- f0
  g1 <- g0
  if(sum(is.na(direction))!=0){
      print('linesearch direction NaN!!!!')
      print('NaN direction:')
      print(direction)
      ret$exist <- FALSE
      ret$wt <- w
      ret$value <- f0
      ret$grad <- g0
      ret$k <- 0
      return(ret)
  }
  delta1 <- delta1
  delta2 <- delta2
  d1 <- as.numeric(t(g0)%*%direction)

  if(is.infinite(d1) || is.na(d1)){
      print('linesearch direction d1 error!')
      print('NaN d1:')
      print(d1)
      ret$exist <- FALSE
      ret$wt <- w
      ret$value <- f0
      ret$grad <- g0
      ret$k <- 0
      return(ret)
  }
 # if(d1 > 0){
 #     print('linesearch is not decrease !!!!')
 #     print('increase d1:')
 #     print(d1)
 #     ret$exist <- FALSE
 #     ret$wt <- w
 #     ret$value <- f0
 #     ret$grad <- g0
 #     ret$k <- 0
 #     return(ret)
 # }
  dg0 <- d1 + lambda*L1_norm(w)
  #print('dg0:')
  #print(dg0)
  #print('lambda*L1_norm(w):')
  #print(lambda*L1_norm(w))
  d1 <- d1 + lambda*(L1_norm(w+direction) - L1_norm(w))
  d1 <- d1*delta1
  
  while(k <= max_linesearch){
    alpha <- beta^k
    f1 <- objFunc(xv=w+alpha*direction, ...)
    f1 <- f1 + lambda*L1_norm(w+alpha*direction)
    if(is.infinite(f1) || is.na(f1)){
        print('line search error!')
        print(f1)
        print('alpha:')
        print(alpha)
    }
    part <- alpha*d1

    k <- k + 1
   # print('f0:')
   # print(f0)
   # print('f1:')
   # print(f1)
   # print('part:')
   # print(part)
    if(f1 <= f0 + part){
      g1 <- derObjFunc(xv=w+alpha*direction, ...)
      dg1 <- as.numeric(t(g1)%*%direction) + alpha*lambda*L1_norm(w+direction)
      #print('dg1:')
      #print(dg1)
      #print('lambda*L1_norm(w+alpha*direction):')
      #print(alpha*lambda*L1_norm(w+direction))
      if(abs(dg1/dg0) <= delta2){
       # print('line search ok~')
        exist <- TRUE
        wt <- w+alpha*direction
        break
      }
    }
  }

  wt <- filter_dir(wt)

  ret$wt <- wt
  ret$exist <- exist
  ret$value <- f1
  ret$grad <- g1
  ret$k <- k

  return(ret)
}

computeL2 <- function(x){
  return(sqrt(sum(x^2)))
}

# compute second derivate of objective funciton of index-th row of the matrix B
derObjbp2 <- function(bv, index, BZ, X, M, lambda2, q, p, n){
  AlphaM <- exp(BZ + matrix(bv,p,1)%*%matrix(M[,index],1,n))
  AlphaMX <- AlphaM + t(X)

  der1 <- -t(digamma(AlphaMX) - digamma(AlphaM)) + (digamma(colSums(AlphaMX)) - digamma(colSums(AlphaM)))
  der1 <- t(der1)

  der2 <- -t(trigamma(AlphaMX) - trigamma(AlphaM)) + (trigamma(colSums(AlphaMX)) - trigamma(colSums(AlphaM)))
  der2 <- t(der2)

  der <- (der1*AlphaM + der2*AlphaM^2)%*%matrix(M[,index]^2,n,1)
  der <- der/n

  return(as.vector(der))
}

# estimate w via proximal quasi-newton method
proximalQusiNewtonB <- function(objFunc, derObjFunc, w, lambda, approx_num, max_linesearch, max_iteration, threshold, delta1_threshold, delta2_threshold, sy_threshold, max_iteration_coor, threshold_coor, ...){
  approx_count <- 0
  approx_index <- 0
  gammat <- 1

  d <- length(w)
  St <- matrix(0,d,approx_num)
  Yt <- matrix(0,d,approx_num)

  round <- 0
  # when SYdot < sy_threshold, change to steepest descent
  steepest_active <- FALSE
  steepest_active_height <- 1
  steepest_count <- 0

  objNew <- 0
  objOld <- objFunc(xv=w, ...) + lambda*L1_norm(w)
  gt0 <- derObjFunc(xv=w, ...)

  ret <- list()

  while(TRUE){
    round <- round + 1

   # print('Round B:')
   # print(round)
    if(!steepest_active){
      BtQ <- computeBtQ(gammat, St, Yt, approx_num, approx_count, approx_index)
      Q <- BtQ$Q
      Qh <- BtQ$Qh
    # adjust gammat and guarantee Bt positive definite
    #gammat <- checkBtPD(gammat, BtQ, round, d)
    #print('gammat adjust:')
   # print(gammat)

    #Bt <- diag(gammat*rep(1,d)) - Q%*%Qh

      direction <- coordinateDescentD(w, gammat, Q, Qh, gt0, lambda, max_iteration = max_iteration_coor, threshold = threshold_coor)
    }else{
      direction <- coordinateDescentD(w, 1, matrix(0,d,2), matrix(0,2,d), gt0, lambda, max_iteration = max_iteration_coor, threshold = threshold_coor)
      steepest_count <- steepest_count + 1
      if(steepest_count >= steepest_active_height){
        steepest_active <- FALSE
      }
    }
    #print('direciton norm1:')
    #print(sum(abs(direction)))
    line <- linesearch(objFunc=objFunc, derObjFunc=derObjFunc, w=w, direction=direction, lambda=lambda, max_linesearch=max_linesearch, f0=objOld, g0=gt0, delta1=delta1_threshold, delta2=delta2_threshold, ...)
    
    exist <- line$exist
    if(!exist){
     # print('line search failed!!!')
      break
    }
    wt <- line$wt
    gt1 <- line$grad
    objNew <- line$value

   # print('line search k:')
   # print(line$k)

   # print('obj new:')
   # print(objNew)

    delta <- objNew - objOld
   # print('delta:')
   # print(delta)

    # update BFGS
    St_per <- wt - w
    Yt_per <- gt1 - gt0
    SYdot <- as.numeric(t(St_per)%*%Yt_per)

    if(SYdot > sy_threshold){
      approx_count <- approx_count + 1
      if(approx_count > approx_num){
        approx_count <- approx_num
      }
    #  print('approx_count:')
    #  print(approx_count)
      approx_index <- approx_index %% approx_num + 1
      St[,approx_index] <- St_per
      Yt[,approx_index] <- Yt_per
      gammat <- SYdot/(t(St_per)%*%St_per)
      gammat <- as.numeric(gammat)
    }else{
      #print(paste('bi sydot <= 0 when round:', as.character(round)))
      steepest_active <- TRUE
      steepest_count <- 0
      print('steepest descent active!')
    }
   # print('gammat:')
   # print(gammat)

    objOld <- objNew
    gt0 <- gt1
    w <- wt

    if(abs(delta) < threshold || round > max_iteration){
     # print('proximal B stop!')
      break
    }
  }

  ret$par <- w
  ret$value <- objOld
  ret$gradient <- gt0
  return(ret)
}

# compute the current objective function for f
computeObjf <- function(X, M, B, Theta, Z, lambda1, lambda2){
  n <- nrow(X)
  B0 <- colMeans(Z)
  Z_center <- t(t(Z) - B0)
  AlphaM <- exp(t(M%*%B) + t(Z))
  #AlphaM <- filterAlphaZero(AlphaM)
  AlphaMX <- AlphaM + t(X)
  S <- t(Z_center)%*%Z_center/n

  obj <- computeAlhpaMatrix(AlphaMX, AlphaM)/n
  obj <- obj - computeLogDet(Theta)/2 + sum(diag(S%*%Theta))/2 + lambda1*sum(abs(Theta))/2 + lambda2*sum(abs(B))
  return(obj)
}

# compute the EBIC
computeEBIC <- function(objf, B, Theta, p, q, n, lambda1, lambda2){
  g <- 0.5
  E1 <- (sum(Theta!=0) - p)/2
  E2 <- sum(B!=0)
  EBIC <- 2*n*(objf - lambda2*sum(abs(B)) - lambda1*sum(abs(Theta))/2) + (E1 + E2)*log(n) + 4*g*E1*log(p) + 2*g*E2*log(p*q)
  return(EBIC)
}

getSparsity <- function(X){
  n <- ncol(X)
  not_zero <- sum(X!=0) - n
  return(not_zero / (n*(n-1)))
}

computeEdgesVary <- function(Theta, Theta_old, B, B_old){
  not1 <- (Theta != 0)
  not2 <- (B != 0)
  not1_old <- (Theta_old != 0)
  not2_old <- (B_old != 0)
  notd1 <- sum(abs(not1_old - not1))/2
  notd2 <- sum(abs(not2_old - not2))
  res <- list()
  res$notd1 <- notd1
  res$notd2 <- notd2
  return(res)
}

# estimate parameters B,b0,Theta,Z on the lambda1 and lambda2 for the
# Lognormal-Dirichlet-Multinomial model
LDM <- function(X, M, n, p, q, B, B0, Theta, Z, lambda1, lambda2, max_iteration, threshold, approx_num_Z, max_linesearch_Z, debug, approx_num_B, max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, threshold_B_coor){
  # the value of objective function on the last step
  objOld <- 0
  # the value of objective function on the current step
  objNew <- 0
  # the round of the iterations
  round <- 0

  objList <- list()
  edges1_list <- list()
  edges2_list <- list()
  edges1_vary_list <- list()
  edges2_vary_list <- list()
  derZ <- matrix(0, n, p)
  B0 <- colMeans(Z)

  while(TRUE){
    round <- round + 1
    print('Round:')
    print(round)

    objOld <- objNew

    B_old <- B
    Theta_old <- Theta

    # estimate the Z
    for(i in 1:n){
      z_i <- Z[i,]
      x_i <- X[i,]
      mu_i <- t(B)%*%M[i,]

      z_i_result <- lbfgs(call_eval=objZi, call_grad=derObjZi, vars=z_i, mu_i=mu_i, x_i=x_i, Theta=Theta, n=n, B0=B0, invisible=1, m=approx_num_Z, max_linesearch=max_linesearch_Z)
      
      Z[i,] <- z_i_result$par
      derZ[i,] <- derObjZi(Z[i,], mu_i, x_i, Theta, n, B0)

      if(i==1){
        print('obj z_i = 1:')
        print(z_i_result$value)
        
      }
    }
    print('norm2 Z:')
    print(computeL2(rep(Z)))
    print('norm2 der Z:')
    print(computeL2(rep(derZ)))
    B0 <- colMeans(Z)

    print('Z max min mean:')
    print(max(Z))
    print(min(Z))
    print(mean(abs(Z)))

    # estimate the B
    Bv <- rep(B,1)
    ZB0 <- t(Z)

    for(i in 1:q){
      print('q:')
      print(i)
      bv <- B[i,]
      Bp <- B
      Bp[i,] <- 0
      BMp <- M%*%Bp
      BZ <- t(BMp) + ZB0

      b_result_per <- proximalQusiNewtonB(objFunc=objbp, derObjFunc=derObjbp, w=bv, lambda=lambda2, approx_num = approx_num_B, max_linesearch = max_linesearch_B, max_iteration = max_iteration_B, threshold = threshold_B, delta1_threshold = delta1_threshold_B, delta2_threshold = delta2_threshold_B, sy_threshold = sy_threshold_B, max_iteration_coor = max_iteration_B_coor, threshold_coor = threshold_B_coor, BZ=BZ, index=i, X=X, M=M, q=q, p=p, n=n)
      B[i,] <- b_result_per$par
      print('B per obj value:')
      print(b_result_per$value)
    }

    objBvalue <- objB(rep(B), ZB0, X, M, lambda2, q, p, n)

    derB <- derObjB(rep(B), ZB0, X, M, lambda2, q, p, n)
    print('norm2 B:')
    print(computeL2(rep(B)))
    print('norm2 derB:')
    print(computeL2(derB))
    print('B max min mean mean_abs:')
    print(max(B))
    print(min(B))
    print(mean(abs(B)))
    print(min(abs(B)))

    # estimate the Theta
    Z_center <- t(t(Z) - B0)
    S <- t(Z_center)%*%Z_center/n
    quic_res <- QUIC(S, rho=lambda1)
    Theta <- quic_res$X

    edges1 <- (sum(Theta!=0)-p)/2
    edges2 <- sum(B!=0)

    der_edges <- computeEdgesVary(Theta, Theta_old, B, B_old)
    print('vary of edges1:')
    print(der_edges$notd1)
    print('vary of edges2:')
    print(der_edges$notd2)
    edges1_vary_list[[round]] <- der_edges$notd1
    edges2_vary_list[[round]] <- der_edges$notd2

    edges1_list[[round]] <- edges1
    edges2_list[[round]] <- edges2

    print('Edges1:')
    print(edges1)
    print('Edges 2:')
    print(edges2)
    # compute the current obj
    objNew <- computeObjf(X, M, B, Theta, Z, lambda1, lambda2)

    print('Objective Function:')
    print(objNew)
    objList[[round]] <- objNew

    delta = abs(objNew - objOld)
    print('delta:')
    print(objNew - objOld)

    print('lambda1:')
    print(lambda1)
    print('lambda2:')
    print(lambda2)

    if(delta < threshold || round >= max_iteration){
      print('iteration stop~')
      break
    }
  }
  EBIC <- computeEBIC(objNew, B, Theta, p, q, n, lambda1, lambda2)

  optimalSolution <- list(B, B0, Theta, Z, lambda1, lambda2, objList, EBIC, edges1_list, edges2_list, edges1_vary_list, edges2_vary_list)
  return(optimalSolution)
}
