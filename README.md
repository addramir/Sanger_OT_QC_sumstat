# QC_sumstat
Quality control tests on summary statistics data

Quality control step to do GWAS summary statistics databases such as GWAS-catalog, FINNGEN, etc..

The first part consist in estimating the priors for each test. I have already calculated these for UK Biobanks and are reported below. The others need to be calucated and scripts will soon be ready to do that


#### P-Z test:
3rd_quartile_beta_prior = 1.4662742614746094e-05

3rd_quartile_beta_se_prior = 1.1408809896806815e-07

3rd_quartile_intercept_prior = 2.4025500351854134e-05

3rd_quartile_se_prior = 6.564901955385949e-08


#### SE-N test:
3rd_quartile_beta_prior = *****

3rd_quartile_beta_se_prior = *****

3rd_quartile_intercept_prior = *****

3rd_quartile_intercept_se_prior = *****


This test was run using a VM with 16 cores.

To make another test run:

gcloud dataproc clusters create ${CLUSTER_NAME} \
        --image-version=2.0 \
        --project=open-targets-genetics-dev \
        --region=${REGION} \
		--master-machine-type=n1-highmem-96 \
        --enable-component-gateway \
        --initialization-actions=gs://genetics_etl_python_playground/initialisation/initialise_cluster.sh \
        --single-node \
        --max-idle=10m

