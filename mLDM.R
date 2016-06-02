####################################################################
# mLDM -- metagenomic Lognormal-Dirichlet-Multinomial Model
# ------------------------------------------------------------------
# return association networks among OTUs and between EFs and OTUs
####################################################################
# Input
# ----------------------------------------
# n -- the number of samples
# p -- the number of OTUs
# q -- the number of environmental factors (EFs)
# X -- n*p matrix, OTU data 
# M -- n*q matrix, Meta data
# max_iteration -- the number of max iterations
# threshold -- the threshold for termination
# approx_num -- the number of gradient vector to approximate the hessian matrix for Z
# max_linesearch -- the number of line search of lbfgs for Z
# model_selection_num -- the number of different lambda to select model, 
#                        the model_selection_num*model_selection_num combinations of lambda1 and lambda2 will be tested
# debug -- true / false, whether print intermediate information
# approx_num_B -- the number of gradient vector to approximate the hessian matrix for B
# max_linesearch_B -- the number of line search of proximal method for B
# max_iteration_B -- the max iterations for B
# threshold_B -- the threshold for termination for B
# delta1_threshold_B and delta2_threshold_B -- the parameters of line search based on strong wolfe condition for B
# sy_threshold_B -- test to maintain positive definite for hessian approximation for B, when < sy_threshold_B, choose steepest descent method
##########################################
# Output
# ----------------------------------------
#
##########################################
library('lbfgs')
library('QUIC')
library('dirmult')
library("psych")
library("MASS")

mLDM <- function(X, M, Z_mean, max_iteration = 1500, threshold = 1e-5, approx_num_Z = 10, max_linesearch_Z = 30, model_selection_num = 6, debug = FALSE, approx_num_B = 10, max_linesearch_B = 30, max_iteration_B = 1000, threshold_B = 1e-5, delta1_threshold_B = 1e-4, delta2_threshold_B = 0.9, sy_threshold_B = 1e-6, max_iteration_B_coor = 20, threshold_B_coor = 1e-6){
#  source('./Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic.R')
  M <- scale(M, center = TRUE, scale = TRUE)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(M)
  
  # set initial values for every variable: Z, B, Theta, B0
  Xr <- X / rowSums(X)
  spearman_ef_otu <- corr.test(x = data.frame(Xr), y = data.frame(M), method="spearman")
  cor_x <- cor(X)
  cor_m <- spearman_ef_otu$r
  # for Z
  #fit <- dirmult(X, epsilon = threshold)
  #print('init Z mean:')
  #print(mean(fit$gamma))
  #Z_init <- mvrnorm(n, mu = rep(mean(fit$gamma),p), Sigma = diag(rep(1,p)))
  #Z_init <- matrix(4,n,p)
  Z_init <- log(X+1) + Z_mean
  print('Z init max min mean:')
  print(max(Z_init))
  print(min(Z_init))
  print(mean(Z_init))
  # for B
  B_init <- t(cor_m)
  # for Theta
  Theta_init <- cor_x
  diag_value <- 1
  while(det(Theta_init) < p){
    diag(Theta_init) <- diag(Theta_init) + rep(diag_value, p)
  }
  Theta_init <- solve(Theta_init)
  # for B0
  B0_init <- colMeans(Z)
  
  # set combinations of lambda1 and lambda2
  cor_x_list <- unique(rep(abs(cor_x)))
  cor_m_list <- unique(rep(abs(cor_m)))
  lambda1_left <- quantile(cor_x, 0.6)
  lambda1_right <- quantile(cor_x, 0.9)
  delta1 <- (lambda1_right - lambda1_left)/model_selection_num
  # lambda1 for Theta
  lambda1_list <- seq(lambda1_left, lambda1_right, delta1)
  
  lambda2_left <- quantile(cor_m, 0.6)
  lambda2_right <- quantile(cor_m, 0.9)
  delta2 <- (lambda2_right - lambda2_left)/model_selection_num
  # lambda2 for B
  lambda2_list <- seq(lambda2_left, lambda2_right, delta2)
  
  # record the minimum of EBIC
  EBIC_min <- 1e+20
  
  LDM_result <- list()
  LDM_result_all <- list()
  count <- 0
  
  length1 <- length(lambda1_list)
  length2 <- length(lambda2_list)
  
  exist_best <- FALSE
  exist_per <- rep(0,length1)
  best_i <- 0
  
  time0 <- proc.time()
  
  for(j in 1:length2){
    for(i in 1:length1){
      solution <- list()
      
      if(!exist_best || exist_per[i]){
        lambda1 <- lambda1_list[i]
        lambda2 <- lambda2_list[j]
        
        solution_EBIC <- 1e+30
        
        solution <- LDM(X, M, n, p, q, B_init, B0_init, Theta_init, Z_init, lambda1, lambda2, max_iteration, threshold, approx_num_Z, max_linesearch_Z, debug, approx_num_B, max_linesearch_B, max_iteration_B, threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B, max_iteration_B_coor, threshold_B_coor)
        
        if(length(solution)!=0){
          solution_EBIC <- solution[[8]]
        }
        
        if(solution_EBIC > 0 && solution_EBIC < EBIC_min){
          EBIC_min <- solution[[8]]
          LDM_result <- solution
          best_i <- i
          if(j == 1 && i == 1){
            B_init <- LDM_result[[1]]
            B0_init <- LDM_result[[2]]
            Theta_init <- LDM_result[[3]]
            Z_init <- LDM_result[[4]]
          }
        }
        print('memory used:')
        print(memory.profile())
        
      }
      
      count <- count + 1
      LDM_result_all[[count]] <- solution
    }
    
    if(best_i!=0){
      exist_per[best_i] <- TRUE
      exist_best <- TRUE
    }
  }
  time1 <- proc.time()
  print('used time for LDM:')
  print(time1 - time0)

  # record all results for LDM model
  LDM_record <- list(LDM_result,LDM_result_all, lambda1_list, lambda2_list)
  
  warnings()
  
  return(LDM_record)
}
