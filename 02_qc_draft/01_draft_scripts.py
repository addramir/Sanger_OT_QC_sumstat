# This script calculate the priors for QC GWAS data ingested 
# Approx 2K UK Biobank studies are considered

from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window
import scipy as sc
from scipy import stats
import numpy as np

#Spark initialization and configuration
spark = SparkSession.builder.master("yarn").getOrCreate()

#functions definition
def calculate_logpval(z2):
        logpval=-np.log10(sc.stats.chi2.sf((z2), 1))
        return(float(logpval))

calculate_logpval_udf = f.udf(calculate_logpval, t.DoubleType())


def logpval(pval):
        logpval=-np.log10(pval)
        return(float(logpval))

logpval_udf = f.udf(logpval, t.DoubleType())


linear_reg_Schema = t.StructType([
    t.StructField("beta", t.FloatType(), False),
    t.StructField("beta_stderr", t.FloatType(), False),
    t.StructField("intercept", t.FloatType(), False),
    t.StructField("intercept_stderr", t.FloatType(), False)])

def calculate_lin_reg(y, x):
    lin_reg = sc.stats.linregress(y,x)
    return [float(lin_reg.slope), float(lin_reg.stderr), float(lin_reg.intercept), float(lin_reg.intercept_stderr)]

lin_udf = f.udf(calculate_lin_reg,  linear_reg_Schema)

def calculate_iqr(spark_dataframe):
    bounds = {
        c: dict(
            zip(["q1", "q3"], spark_dataframe.approxQuantile(c, [0.25, 0.75], 0))
        )
        for c in spark_dataframe.columns
    }
    for c in bounds:
        iqr = bounds[c]['q3'] - bounds[c]['q1']
        bounds[c]['3iqr'] = 3*iqr
    return bounds

 
#main script

path_ukbio = "gs://genetics-portal-dev-sumstats/filtered/significant_window_2mb/gwas/NEALE2_100001_raw.parquet"
NEALE_studies = spark.read.parquet(path_ukbio)

grp_window = Window.partitionBy('study_id')

NEALE_studies_columns = (NEALE_studies
    .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
    .withColumn("zscore", f.col("beta")/f.col("se"))
    .withColumn("new_logpval", calculate_logpval_udf(f.col("zscore")**2))
    .withColumn("logpval",logpval_udf(f.col("pval")))
    .withColumn("var_af", 2*(f.col("eaf") * (1-f.col("eaf"))))
    .withColumn("pheno_var", ((f.col("se")**2) * (f.col("n_total") * f.col("var_af"))) + ((f.col("beta")**2) * f.col("var_af")))
    .withColumn("total_SNP", f.count("id").over(grp_window))
    )

PHENO_MEDIAN=float(NEALE_studies_columns.agg(f.percentile_approx("pheno_var",0.5)).toPandas().iloc[0,0])
total_SNP=int(NEALE_studies_columns.agg(f.percentile_approx("total_SNP",0.5)).toPandas().iloc[0,0])
n_total=int(NEALE_studies_columns.agg(f.percentile_approx("n_total",0.5)).toPandas().iloc[0,0])

NEALE_studies_columns = (NEALE_studies_columns  
    .withColumn("N_hat",(PHENO_MEDIAN - ((f.col("beta")**2) * f.col("var_af"))/((f.col("se")**2) * f.col("var_af")))))


#NEALE_studies_columns.agg("result_lin_reg", lin_udf(f.col("logpval"), f.col("new_logpval")))

NEALE_studies_res = (NEALE_studies_columns 
    .agg(
       f.collect_list("logpval").alias("pval_vector"),
       f.collect_list("new_logpval").alias("new_pval_vector"),
       f.percentile_approx( f.col("N_hat")/f.col("n_total"), 0.5).alias("median_N"),
       f.stddev(f.col("N_hat")/f.col("n_total")).alias("se_N")
       )
    .withColumn("result_lin_reg", lin_udf(f.col("pval_vector"), f.col("new_pval_vector"))))


NEALE_studies_nfe = (NEALE_studies
    .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
    .select("study_id", "type", "id", "beta", "eaf")
    )

NEALE_studies_nfe = NEALE_studies_nfe.join(variant_annotation, on = "id", how = "inner")

grp_window = Window.partitionBy('study_id')
NEALE_studies_nfe = (NEALE_studies_nfe
                    .withColumn("delta_freq", f.abs(f.col("nfe_freq") - f.col("eaf")))
                    .withColumn("total_SNP", f.count("id").over(grp_window))
                    .filter(f.col("delta_freq")> 0.2)
                    .withColumn("over_threshold", f.count("id").over(grp_window))
                    .withColumn("delta_prop", f.col("over_threshold")/f.col("total_SNP"))
                    .select("study_id", "delta_prop", "total_SNP", "over_threshold")
                    .distinct()
                    )

NEALE_studies_nfe = (NEALE_studies_nfe
                    .withColumn("total_SNP", f.count("id").over(grp_window)))






NEALE_studies_beta = float(NEALE_studies_columns.select(f.mean("beta")).toPandas().iloc[0,0])

