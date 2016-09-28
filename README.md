# mLDM
mLDM: a new hierarchical Bayesian statistical model for sparse microbial association discovery

##USAGE:
```
  1. R packages: lbfgs, QUIC, Rcpp, RcppEigen should be installed first!
  2. two approaches to use mLDM: 
    2.1  a. Download the mLDM.cpp
         b. library(Rcpp), library(RcppEigen), library(lbfgs), library(QUIC), and sourceCpp("mLDM.cpp")

    2.2  a. Download the R package mLDM_1.0.tar.gz and install.packages('PATH/mLDM_1.0.tar.gz', repos=NULL, type="source")
         b. library(mLDM)
         c. help(mLDM)
```

### R package:
```
  mLDM_1.0.tar.gz
```

###Input Parameters for mLDM:
```
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
  approx_num_B -- the number of gradient vector to approximate the hessian matrix for B 
  max_linesearch_B -- the number of line search of proximal method for B
  max_iteration_B -- the max iterations for B 
  threshold_B -- the threshold for termination for B 
  delta1_threshold_B and delta2_threshold_B -- the parameters of line search based on strong wolfe condition for B 
  sy_threshold_B -- test to maintain positive definite for hessian approximation for B, when < sy_threshold_B, choose steepest  descent method
  max_iteration_B_coor -- The max iteration of coordinate descent when optimize B
  threshold_B_coor -- Stop the coordinate descent when the the variation of the direction samll than the threshold
  ratio1 -- Set ratio1 to control minimum values of lambda1 and lambda2.
  ratio2 -- Set ratio2 to control maximum values of lambda1 and lambda2.
  verbose -- Set FALSE to run mLDM at a silent mode; set TRUE to see the debug information from the program.
```
###Output Parameters for mLDM:
return a list consists of optimal and all results used in model selection
```
  a list consists of optimal and all results from mLDM are returned:
  list$optimal -- the optimal result via the model selection
  list$all -- all results corresponding to different lambda1 and lambda2
  list$lambda1 -- the list of all lambda1
  list$lambda2 -- the list of all lambda2


```
### TARA
#####TARA_Validation_dataset.RData
  Selected subset data (67 OTUs, 17 EFs, 221 Samples and 28 known genus-level associations) of original TARA data for validation
```
   X   N*P matrix  N samples, P otus  OTU table (N = 221, P = 67)
   M   N*Q matrix  N samples, P environmental factors (Q = 17)
   otus_ids  P*1 vector  Selected otus' rank in original TARA OTU table
   otus_annotations  P*1 list  Selected otus' annotations
   genus_associations_ids  G*1 vector  Selected genus-level associations IDs in original TARA known associations table (G = 28)
   genus_associations_names  G*2 list  Two genus list of selected genus-level associations
   genus_associations_pair  G*3 matrix  
        -- 1-th column  association IDs
        -- 2-th column  OTU-1 id
        -- 3-th column  OTU-2 id
```
#####Suplemental Tabel 4.csv
Selected 28 known genus-level interactions and corresponding 67 OTUs. This table is subset of the original provided by TARA OCEANS project website (http://www.raeslab.org/companion/ocean-interactome.html).
##### result-67-17-221-3-1-1-TARA-67-filter.RData
Estimated associations by mLDM on TARA Oceans dataset.

###Colorectal Cancer
#####Baxter_CRC.RData
```
   X   N*P matrix  N samples, P otus  OTU table (N = 490, P = 117)
   M   N*Q matrix  N samples, P environmental factors (Q = 13)
       We encode 'site' (4 cities), 'Dx_bin' (5 diagnosis states) and 'Gender' (male and female) with 0-1 coding and get 4, 5, 2 features respectively.
   X_name  OTU numbers for 117 OTUs
   X_tax   OTU taxonomy for 117 OTUs
```
#####glne007.csv
OTU table used by mLDM
#####metadata.csv
Meta data used by mLDM

###Human Microbiome Project
#####HMP-All.RData
```
   X   N*P matrix  N samples, P otus  OTU table (N = 112, P = 110)
   M   N*Q matrix  N samples, P environmental factors (Q = 3)
       We encode 'Gender' with 0-1 coding and get 2 features.
```
#####HMP_2_groups.RData
```
   X_two   a list consist of two subsets of 112 samples
           we distribute samples according to their visit number of the same subject.
           every subset has 95 samples
   M_two   meta data related to two subsets
```
#####v13_map_uniquebyPDN_stool.csv
   Human stool samples we selected

###West English Channel dataset 
######48 not zero data.RData
Subset of West English Channel dataset selected by mLDM
```
   X  N*P matrix N samples, P otus   OTU table (N=47, P=48)
   M  N*Q matrix N samples, Q environmental factors (Q=8)
   X_class  OTUs' annotation on phylum level
   X_name   OTU labels for P=48 OTUs
   X_tax    OTUs' annotation in detail
   M_name   names of environmental factors
```
