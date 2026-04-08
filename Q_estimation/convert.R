library(radiator)

setwd("/local/workdir/yc2644/CV_CT_array/centroid_sim")

genomic_converter(
  data="merged.gen",
  strata = NULL,
  output = c("vcf"),
  filename = "geno_n1601_in_silico_wild_AQ_fp_noinvers",
  parallel.core = 4, #parallel::detectCores() - 1,
  verbose = TRUE,
)

genomic_converter(
  data="test_fp_noinvers.gen",
  strata = NULL,
  output = c("vcf"),
  filename = "test_fp_noinvers",
  parallel.core = 4, #parallel::detectCores() - 1,
  verbose = TRUE,
)
