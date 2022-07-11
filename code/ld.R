library(LDlinkR)
library(vcfR)

genotypes <- read.vcfR("data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz")


LDproxy_batch(snp=random.hits,
              pop="YRI",
              r2d="r2",
              token="36e356ef3c1e")

