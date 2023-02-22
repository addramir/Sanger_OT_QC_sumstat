from pyspark.sql import SparkSession
import pyspark.sql.functions as f
import pyspark.sql.types as t
from pyspark.sql import Window
import scipy as sc
from scipy import stats
import numpy as np
import pandas as pd

#Spark initialization and configuration

global spark
spark = (
    SparkSession.builder
    .master('local[*]')
    .config('spark.driver.memory', '5g')
    .config("spark.driver.cores", "1")
    .config("spark.executor.cores", "1") 
    .appName('spark')
    .getOrCreate()
)

gwas_list=!ls /mnt/disks/gwas/SS_QC/

i=0
gw=gwas_list[i]
print("N "+str(i)+": "+gw)
GWAS=spark.read.parquet("/mnt/disks/gwas/SS_QC/"+gw)
out=GWAS.toPandas()

for i,gw in enumerate(gwas_list[1:]):
    print("N "+str(i)+": "+gw)
    GWAS=spark.read.parquet("/mnt/disks/gwas/SS_QC/"+gw)
    GWAS=GWAS.toPandas()
    out=pd.concat([out,GWAS])

out["beta_PZ"]=[con[0] for con in (out.iloc[:,6])]
out["se_beta_PZ"]=[con[1] for con in (out.iloc[:,6])]
out["intercept_PZ"]=[con[2] for con in (out.iloc[:,6])]
out["se_intercept_PZ"]=[con[3] for con in (out.iloc[:,6])]

out.to_csv("~/projects/Sanger_OT_QC_sumstat/02_qc_draft/QC_results.csv",index=False,sep=";")