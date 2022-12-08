# This script calculate the beta and se using the estimated and reported pvalues in uk biobank studies. 
# These are used as priors to do QC on other GWAS studies
# Approx 2K UK Biobank studies are considered

from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
import scipy as sc
from scipy import stats
import numpy as np

#Spark initialization and configuration
spark = SparkSession.builder.master("yarn").getOrCreate()

#functions definition
def calculate_pval(z_score):
       return(float(sc.stats.chi2.sf((z_score), 1)))

calculate_pval_udf = f.udf(calculate_pval, t.DoubleType())

def calculate_lin_reg(y, x):
    lin_reg = sc.stats.linregress(y,x)
    return [float(lin_reg.slope), float(lin_reg.stderr), float(lin_reg.intercept), float(lin_reg.intercept_stderr)]

lin_udf = f.udf(calculate_lin_reg,  linear_reg_Schema)

def calculate_iqr(spark_dataframe):
    bounds = {
        c: dict(
            zip(["q1", "q3"], df.approxQuantile(c, [0.25, 0.75], 0))
        )
        for c in df.columns
    }
    for c in bounds:
        iqr = bounds[c]['q3'] - bounds[c]['q1']
        bounds[c]['3iqr'] = 3*iqr
    return bounds

 
#main script

path_ukbio = "gs://genetics-portal-dev-sumstats/filtered/significant_window_2mb/gwas/NEALE2_*"
NEALE_studies = spark.read.parquet(path_ukbio)

NEALE_studies = (
    NEALE_studies
        .withColumn("zscore", f.col("beta")/f.col("se"))
        .withColumn("new_pval", calculate_pval_udf(f.col("zscore")**2))
        .select("study_id","pval", "new_pval")
            )

linear_reg_Schema = t.StructType([
    t.StructField("beta", t.FloatType(), False),
    t.StructField("beta_stderr", t.FloatType(), False),
    t.StructField("intercept", t.FloatType(), False),
    t.StructField("intercept_stderr", t.FloatType(), False)])

final_data_priors = (NEALE_studies
              .groupBy("study_id")
              .agg(f.collect_list("pval").alias("pval_vector"),
                   f.collect_list("new_pval").alias("new_pval_vector")
                   )
              .withColumn("result_lin_reg", lin_udf(f.col("pval_vector"), f.col("new_pval_vector")))
            )


#Given that this script takes 2 hours to run it is advised to save the results of this run
#final_data.write.parquet("/home/ba13/prior_ukbio_pz")

prior_uk_bio_lr = prior_uk_bio.select("result_lin_reg.*")
prior_uk_bio_bounds = calculate_bounds(prior_uk_bio_lr)
print(prior_uk_bio_bounds)

#{'beta': {'q1': 0.9999961853027344, 'q3': 1.000001072883606, '3iqr': 1.4662742614746094e-05}, 'beta_stderr': {'q1': 1.3771888518476771e-08, 'q3': 5.1801254841166156e-08, '3iqr': 1.1408809896806815e-07}, 'intercept': {'q1': -8.286592674267013e-07, 'q3': 7.17984084985801e-06, '3iqr': 2.4025500351854134e-05}, 'intercept_stderr': {'q1': 7.51434381385252e-09, 'q3': 2.9397350331805683e-08, '3iqr': 6.564901955385949e-08}}
