from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window
import scipy as sc
from scipy import stats
import numpy as np

#Spark initialization and configuration
spark = SparkSession.builder.master("yarn").getOrCreate()

gwas_cat_data = (spark.read.parquet("gs://genetics-portal-dev-sumstats/filtered/significant_window_2mb/gwas/GCST000*.parquet")
                .withColumn("id", f.concat_ws("_", f.col("chrom"), f.col("pos"), f.col("ref"), f.col("alt")))
                .select("study_id", "type", "id", "beta", "eaf")
                )

variant_annotation = (spark.read.parquet("gs://genetics_etl_python_playground/XX.XX/output/python_etl/parquet/variant_annotation")
                        .select("id", f.col("alleleFrequencies.alleleFrequency")[5].alias("nfe_freq"))
                        .filter(~f.col("nfe_freq").isNull())
                        )


grp_window = Window.partitionBy('study_id')


gwas_cat_data_nfe = gwas_cat_data.join(variant_annotation, on = "id", how = "inner")
gwas_cat_data_nfe = (gwas_cat_data_nfe
                    .withColumn("delta_freq", f.abs(f.col("nfe_freq") - f.col("eaf")))
                    .withColumn("total_SNP", f.count("id").over(grp_window))
                    .filter(f.col("delta_freq")> 0.2)
                    .withColumn("over_threshold", f.count("id").over(grp_window))
                    .withColumn("delta_prop", f.col("over_threshold")/f.col("total_SNP"))
                    .select("study_id", "delta_prop", "total_SNP", "over_threshold")
                    .distinct()
                    .filter(f.col("delta_prop") > 0.2)
                    )

gwas_cat_data_nfe.write.parquet("/home/ba13/AF-AF_filter")
