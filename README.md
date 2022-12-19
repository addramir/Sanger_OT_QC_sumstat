# QC_sumstat

Quality control step to preform on GWAS-catalog summary statistics data.


## Steps:

### 1. AF-AF_check.py

Here first we estimates the differences in allele frequencies higher than 0.2 using gnomad nfe as reference. then wetae calculate the proportion of of allele for each study that pass this threshold. 

#### Studies that fail this step need to be removed

### 2. mean_Beta.py

Work in progress !

#### Studies that fail this step need to be removed



### 3 and 4. PZ_N_priors.py

Here first we estimate the priors which are the gold standard thresholds that needs to be respected using the UK Biobank dataset.

#### Studies that fail this step need to replace their p-value with the one calculated using the betas and se from the study

#### P-Z test priors values pre-computed:

3rd_quartile_beta_prior_pz = 1.4662742614746094e-05

3rd_quartile_beta_se_prior_pz = 1.1408809896806815e-07

3rd_quartile_intercept_prior_pz = 2.4025500351854134e-05

3rd_quartile_se_prior_pz = 6.564901955385949e-08


### SE-N test:

#### Studies that fail this step needs to be discarded from further analyses

3rd_quartile_beta_prior = *****

3rd_quartile_beta_se_prior = *****

3rd_quartile_intercept_prior = *****

3rd_quartile_intercept_se_prior = *****


#### To run the PZ_N_priors test we need a bigger machine:

gcloud dataproc clusters create QC_gwas \
        --image-version=2.0 \
        --project=open-targets-genetics-dev \
        --region= europe-west1-b \
	--master-machine-type=n1-highmem-96 \
        --enable-component-gateway \
        --single-node \
        --max-idle=10m
	
	
gcloud dataproc jobs submit pyspark --cluster=QC_gwas --project=open-targets-genetics-dev --region=europe-west1-b PZ_N_priors.py

