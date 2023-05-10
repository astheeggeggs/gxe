import hail as hl
import sys

# module load Anaconda3/2022.05
# module load OpenBLAS/0.3.8-GCC-9.2.0
# export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.8-GCC-9.2.0/lib/libopenblas.so
# module load java/1.8.0_latest

chrom = sys.argv[1]

hl.init(spark_conf={'spark.driver.memory': '165g'},tmp_dir="/well/lindgren/UKBIOBANK/dpalmer/tmp/")
mt = hl.read_matrix_table(f"/well/lindgren/barney/gxe/hail/ukb_imp_chr{chrom}_v3.mt")

# Import and filter samples based on superpopulation labels
ancs = hl.import_table("/well/lindgren/barney/gxe/data/superpopulation_labels.tsv", force=True, impute=True, missing=['', 'NA'])
ancs = ancs.key_by("sample.ID")
ancs = ancs.select(ancs.classification_strict)
mt = mt.annotate_cols(classification_strict = ancs[mt.s].classification_strict)
mt = mt.filter_cols(mt.classification_strict == "EUR")

# Perform sample quality control and filter based on call rate
mt = hl.sample_qc(mt)
mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.97)

# Import and annotate variant information
infos = hl.import_table(f"/well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_mfi_chr{chrom}_v3.txt", force=True, impute=True, no_header=True, missing=['', 'NA'])
infos = infos.drop(*["f1","f2","f3","f4"])
infos = infos.rename({"f0": "varid","f5": "MAF", "f6":"minor_allele", "f7": "info_score"})
infos = infos.key_by("varid")
mt = mt.annotate_rows(**infos[mt.varid])

# remove related samples
unrelated = hl.import_table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_non_finnish_eur_unrelated.txt", no_header=True, key="f0")

mt = mt.filter_cols(hl.is_defined(unrelated[mt.s]))

# Filter variants based on MAF and INFO score
mt = mt.filter_rows((mt.MAF > 0.001) & (mt.info_score > 0.8))

# Filter to the subset of samples for which you have PRS information.

# Perform variant quality control and filter based on allele frequency and HWE p-value
mt = hl.variant_qc(mt)
mt = mt.filter_rows(
    (mt.variant_qc.AF[0] > 0.005) & 
    (mt.variant_qc.AF[1] > 0.005) &
    (mt.variant_qc.p_value_hwe > 1e-10))

# Count the matrix table
mt.count()

# Write to plink
hl.export_plink(mt, "/well/lindgren/UKBIOBANK/dpalmer/gxe/ukb_chr_{chrom}", ind_id=mt.s)
