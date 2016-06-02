# mLDM
mLDM: a new hierarchical Bayesian statistical model for sparse microbial association discovery

## TARA
### TARA_Validation_dataset.RData
            > Selected subset data (67 OTUs, 17 EFs, 221 Samples and 28 known genus-level associations) of original TARA data for validation
            -- X   N*P matrix  N samples, P otus  OTU table (N = 221, P = 67)
            -- M   N*Q matrix  N samples, P environmental factors (Q = 17)
            -- otus_ids  P*1 vector  Selected otus' rank in original TARA OTU table
            -- otus_annotations  P*1 list  Selected otus' annotations
            -- genus_associations_ids  G*1 vector  Selected genus-level associations IDs in original TARA known associations table (G = 28)
            -- genus_associations_names  G*2 list  Two genus list of selected genus-level associations
            -- genus_associations_pair  G*3 matrix  
            -- first column  association IDs
            -- second column  OTU-1 id
            -- third column  OTU-2 id
