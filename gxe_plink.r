library(data.table)
library(dplyr)

# export PATH="/well/lindgren/dpalmer/:$PATH"

# Import phenotype data
chr_PRS_file <- "/well/lindgren/barney/gxe/data/BMI_imp_int_pgs_chrom.txt.gz"
dt_chr <- fread(cmd = paste("zcat", chr_PRS_file))
for (chr in seq(1,22)) {
	dt_chr[, c(paste0("off_chr", chr)) := Reduce(`+`, .SD), .SDcol = seq(2,23)[-(chr+1)]]
}

# Initial covariate file, write to disk
# FID and IID in first two columns, covariates in remaining columns
fwrite(dt_chr[, .(sid, sid, off_chr22)],
	file="/well/lindgren/dpalmer/test_off_chr22.tsv",
	sep='\t', col.names=FALSE)

# Filter phenotype file to a single phenotype
phenotype_file <- "/well/lindgren/barney/gxe/data/curated_phenotypes_cts.tsv"
dt_pheno <- fread(phenotype_file)

# The first and second columns of that file must contain family and within-family IDs, respectively.
fwrite((dt_pheno %>% filter(!is.na(BMI_imp)))[, .(eid, eid, BMI_imp)],
	file="/well/lindgren/dpalmer/test_BMI_imp.tsv",
	sep='\t', col.names=FALSE)

# The first and second columns of that file must contain family and within-family IDs, respectively.
bed <- "/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr22_v2.bed"
bim <- "/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr22_v2.bim"
fam <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"

covar <- "/well/lindgren/dpalmer/test_off_chr22.tsv"
pheno <- "/well/lindgren/dpalmer/test_BMI_imp.tsv"
out <- "testing_1_2_3"
# Run plink, filtering to these samples and writing the result
system(paste("plink",
	"--bim", bim,
	"--bed", bed,
	"--fam", fam,
	"--pheno", pheno,
	"--linear interaction standard-beta hide-covar",
	"--covar", covar,
	"--tests 3",
	"--out", out))

system(paste("plink2",
	"--bim", bim,
	"--bed", bed,
	"--fam", fam,
	"--pheno", pheno,
	"--linear interaction hide-covar",
	"--variance-standardize",
	"--covar", covar,
	"--tests 3",
	"--out", out))


# Good, now do the same thing using hail, and make sure it's the same.

# Now, include all the covariates and rerun, and ensure that it's the same.

# Prepare all of the plink files from the .bgens, filtering to the required samples for running plink.

# Run them