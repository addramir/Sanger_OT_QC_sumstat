from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window
import scipy as sc
from scipy import stats
import numpy as np
import pandas as pd

gwas_list=!gsutil ls gs://genetics-portal-dev-sumstats/unfiltered/gwas
gwas_list=gwas_list[1:]

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

 
#main script

gw=gwas_list[0]
i=0
out=pd.DataFrame()
for i,gw in enumerate(gwas_list):
	print("N "+str(i)+": "+gw)
	GWAS=spark.read.parquet(gw)
	grp_window = Window.partitionBy('study_id')
	GWAS_columns = (GWAS
	    .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
	    .withColumn("zscore", f.col("beta")/f.col("se"))
	    .withColumn("new_logpval", calculate_logpval_udf(f.col("zscore")**2))
	    .withColumn("logpval",logpval_udf(f.col("pval")))
	    .withColumn("var_af", 2*(f.col("eaf") * (1-f.col("eaf"))))
	    .withColumn("pheno_var", ((f.col("se")**2) * (f.col("n_total") * f.col("var_af"))) + ((f.col("beta")**2) * f.col("var_af")))
	    .withColumn("total_SNP", f.count("id").over(grp_window))
	    )
	PHENO_MEDIAN=float(GWAS_columns.agg(f.percentile_approx("pheno_var",0.5)).toPandas().iloc[0,0])
	total_SNP=int(GWAS_columns.agg(f.percentile_approx("total_SNP",0.5)).toPandas().iloc[0,0])
	n_total=int(GWAS_columns.agg(f.percentile_approx("n_total",0.5)).toPandas().iloc[0,0])
	GWAS_columns = (GWAS_columns  
	    .withColumn("N_hat",((PHENO_MEDIAN - ((f.col("beta")**2) * f.col("var_af")))/((f.col("se")**2) * f.col("var_af")))))
	GWAS_columns_res = (GWAS_columns 
	    .agg(
	       f.collect_list("logpval").alias("pval_vector"),
	       f.collect_list("new_logpval").alias("new_pval_vector"),
	       f.percentile_approx( f.col("N_hat")/f.col("n_total"), 0.5).alias("median_N"),
	       f.stddev(f.col("N_hat")/f.col("n_total")).alias("se_N")
	       )
	    .withColumn("result_lin_reg", lin_udf(f.col("pval_vector"), f.col("new_pval_vector")))
	    .select('median_N', 'se_N', 'result_lin_reg')
	    )
	GWAS_columns_res=GWAS_columns_res.toPandas()
	GWAS_beta = float(GWAS_columns.select(f.mean("beta")).toPandas().iloc[0,0])
	out_r=pd.DataFrame({"gw":gw,"total_SNP":total_SNP,"n_total":n_total,"GWAS_beta":GWAS_beta},index=[i])
	out_r=pd.concat([out_r,GWAS_columns_res],axis=1)
	out=out.append(out_r)


