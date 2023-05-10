import hail as hl
import pandas as pd

# module load Anaconda3/2022.05
# module load OpenBLAS/0.3.8-GCC-9.2.0
# export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.8-GCC-9.2.0/lib/libopenblas.so

hl.init(spark_conf={'spark.driver.memory': '165g'},tmp_dir="/well/lindgren/dpalmer/tmp/")

import os
cpu_count = os.cpu_count()

bed = "/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr22_v2.bed"
bim = "/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr22_v2.bim"
fam = "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
mt = hl.import_plink(bed=bed, bim=bim, fam=fam, reference_genome='GRCh37')

# Import phenotype data
chr_PRS_file = "/well/lindgren/dpalmer/test_off_chr22.tsv"
phenotype_file = "/well/lindgren/dpalmer/test_BMI_imp.tsv"

ht = hl.import_table(chr_PRS_file, impute=True, no_header=True, types={"f0": hl.tstr}, missing=['', 'NA'], key='f0')
ht = ht.transmute(PRS=ht.f2).select("PRS")

# Import and join phenotype data
ht_phe = hl.import_table(phenotype_file, impute=True, no_header=True, types={"f0": hl.tstr}, missing=['', 'NA'], key='f0')
ht_phe = ht_phe.transmute(BMI_imp=hl.float(ht_phe.f2)).select("BMI_imp")

# Normalize genotypes in the MatrixTable
mt = hl.experimental.ldscsim.normalize_genotypes(mt.GT)

# Standardize y
stats = ht_phe.aggregate(hl.agg.stats(ht_phe["BMI_imp"]))
ht_phe = ht_phe.annotate(BMI_imp = (ht_phe["BMI_imp"] - stats.mean) / stats.stdev)

mt = mt.annotate_cols(PRS = ht[mt.s]["PRS"])
mt = mt.annotate_cols(BMI_imp = ht_phe[mt.s]["BMI_imp"])

# Compute gene-environment interaction term by multiplying the normalized genotypes and off-chromosome polygenic scores
mt = mt.select_entries(GxE = mt.PRS * mt.norm_gt)
covariates_to_use = ["PRS"]

# Perform linear regression using the gene-environment interaction term as the predictor
gwas = hl.linear_regression_rows(
    y=mt.BMI_imp, 
    x=mt.GxE,
    covariates=[1] + [mt[x] for x in covariates_to_use])

gwas.export(f'/well/lindgren/barney/gxe/data/gwas_{chrom}.tsv.bgz')
