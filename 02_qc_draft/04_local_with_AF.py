from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window
import scipy as sc
from scipy import stats
import numpy as np


#Spark initialization and configuration

global spark
spark = (
    SparkSession.builder
    .master('local[*]')
    .config('spark.driver.memory', '15g')
    .appName('spark')
    .getOrCreate()
)
#spark = SparkSe

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

#variant_annotation = (spark.read.parquet("/mnt/disks/gwas/variant-index/"))
#variant_annotation = (variant_annotation
#                .withColumn("id", f.concat_ws("_", f.col("chr_id"), f.col("position"), f.col("ref_allele"), f.col("alt_allele")))
#                )
#variant_annotation = (variant_annotation
#	.filter(f.col("gnomad_nfe")>=0.0001)
#	.filter(f.col("gnomad_nfe")<=0.9999)
#	.select("id","gnomad_nfe")
#	)
#variant_annotation.write.parquet("/mnt/disks/gwas/variant_annotation_gnomad_nfe")
variant_annotation = spark.read.parquet("/mnt/disks/gwas/variant_annotation_gnomad_nfe")

 
#main script

gwas_list=!ls /mnt/disks/gwas/raw_20230120
exist_list=!ls /mnt/disks/gwas/SS_QC/
ll=set(gwas_list)-set(exist_list)

#gw=ll[0]
#i=0
for i,gw in enumerate(ll):
	print("N "+str(i)+": "+gw)
	GWAS=spark.read.parquet("/mnt/disks/gwas/raw_20230120/"+gw)
	grp_window = Window.partitionBy('study_id')
	GWAS_columns = (GWAS
		.filter(f.col("beta")!=0)
		.filter(f.col("se")>0)
		.withColumn("zscore", f.col("beta")/f.col("se"))
		.filter((f.col("zscore")**2)<=500)
		.withColumn("new_logpval", calculate_logpval_udf(f.col("zscore")**2))
		.withColumn("logpval",logpval_udf(f.col("pval")))
		.withColumn("var_af", 2*(f.col("eaf") * (1-f.col("eaf"))))
		.withColumn("pheno_var", ((f.col("se")**2) * (f.col("n_total") * f.col("var_af"))) + ((f.col("beta")**2) * f.col("var_af")))
		.withColumn("pheno_median", f.percentile_approx( "pheno_var", 0.5).over(grp_window))
		.withColumn("N_hat",((f.col("pheno_median") - ((f.col("beta")**2) * f.col("var_af")))/((f.col("se")**2) * f.col("var_af"))))
		.groupBy("study_id")
		.agg(
			f.collect_list("logpval").alias("pval_vector"),
			f.collect_list("new_logpval").alias("new_pval_vector"),
			f.percentile_approx( f.col("N_hat")/f.col("n_total"), 0.5).alias("median_N"),
			f.stddev(f.col("N_hat")/f.col("n_total")).alias("se_N"),
			f.mean("n_total",).alias("N"),
			f.count("n_total").alias("total_SNP"),
			f.mean("beta").alias("mean_beta")
		)
		.withColumn("result_lin_reg", lin_udf(f.col("pval_vector"), f.col("new_pval_vector")))
		.select('study_id','median_N','se_N','N','total_SNP','mean_beta','result_lin_reg')
	)
	GWAS_eaf = (GWAS
				.filter(f.col("beta")!=0)
				.filter(f.col("se")>0)
				.withColumn("zscore", f.col("beta")/f.col("se"))
				.filter((f.col("zscore")**2)<=500)
                .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
                .select("study_id", "type", "id", "beta", "eaf")
                )
	GWAS_eaf = GWAS_eaf.join(variant_annotation, on = "id", how = "inner")
	GWAS_nfe = (GWAS_eaf
                    .withColumn("delta_freq", f.abs(f.col("gnomad_nfe") - f.col("eaf")))
                    .withColumn("total_SNP", f.count("id").over(grp_window))
                    .filter(f.col("delta_freq")> 0.2)
                    .withColumn("over_threshold", f.count("id").over(grp_window))
                    .withColumn("delta_prop", f.col("over_threshold")/f.col("total_SNP"))
                    .select("delta_prop","over_threshold")
                    .distinct()
                    )
	GWAS_columns=GWAS_columns.join(GWAS_nfe)
	GWAS_columns.write.parquet("/mnt/disks/gwas/SS_QC/"+gw)
