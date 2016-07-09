# mLDM
mLDM: a new hierarchical Bayesian statistical model for sparse microbial association discovery

Usage:
  1. R packages: lbfgs, QUIC, dirmult, psych, MASS should be installed first!
  2. Download the mLDM.R and Lognormal-Dirichlet-Multinomial-lbfgs-proximal-split-q-active-set-quic.R
  3. Source these two files before running mLDM

Input for mLDM.R:
   n -- the number of samples
   p -- the number of OTUs
  q -- the number of environmental factors (EFs)
  X -- n*p matrix, OTU data 
  M -- n*q matrix, Meta data 
  Z_mean -- a positive integer for initalization for latent variable Z
            default is 1, but need to set Z_mean a little bit large when 
            the biggest OTU is >> the smallest OTU, try to maintain the 
            minimum of latent variable Z >= 0
  max_iteration -- the number of max iterations
  threshold -- the threshold for termination
  approx_num -- the number of gradient vector to approximate the hessian matrix for Z
  max_linesearch -- the number of line search of lbfgs for Z
  model_selection_num -- the number of different lambda to select model, 
                         the model_selection_num*model_selection_num combinations of lambda1 and lambda2 will be tested
  debug -- true / false, whether print intermediate information
  approx_num_B -- the number of gradient vector to approximate the hessian matrix for B
  max_linesearch_B -- the number of line search of proximal method for B
  max_iteration_B -- the max iterations for B
  threshold_B -- the threshold for termination for B
  delta1_threshold_B and delta2_threshold_B -- the parameters of line search based on strong wolfe condition for B
  sy_threshold_B -- test to maintain positive definite for hessian approximation for B, when < sy_threshold_B, choose steepest descent method

## TARA
### TARA_Validation_dataset.RData
 * Selected subset data (67 OTUs, 17 EFs, 221 Samples and 28 known genus-level associations) of original TARA data for validation

            -- X   N*P matrix  N samples, P otus  OTU table (N = 221, P = 67)
            -- M   N*Q matrix  N samples, P environmental factors (Q = 17)
            -- otus_ids  P*1 vector  Selected otus' rank in original TARA OTU table
            -- otus_annotations  P*1 list  Selected otus' annotations
            -- genus_associations_ids  G*1 vector  Selected genus-level associations IDs in original TARA known associations table (G = 28)
            -- genus_associations_names  G*2 list  Two genus list of selected genus-level associations
            -- genus_associations_pair  G*3 matrix  
                -- 1-th column  association IDs
                -- 2-th column  OTU-1 id
                -- 3-th column  OTU-2 id
