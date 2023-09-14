# Comparing Phasing Methods with Pseudo-diploids

This repository contains the code necessary to perform the analyses contained in *something*. 

## VCF Pre-processing

In our analyses, we make use of the [1000 Genomes Project (1kGP) 30x deep sequencing release](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) from the NYGC. In particular, we make use of the X chromosome VCF from the `20220422_3202_phased_SNV_INDEL_SV` directory. 

Our pre-processing consists of:

1. subsetting the dataset down to the 2,504 unrelated individuals from the phase 3 1kGP release
2. applying the 1kGP pilot accessibility mask
3. Removing PAR

Remove non-SNPs, remove PAR, apply mask:

```
bcftools view -Ov --types snps -S unrel_ids.txt 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz | bcftools view -Ov -t ^chrX:10001-2781479 | bcftools view -Ov -t ^chrX:155701382-156030895 | bedtools intersect -a - -b 20160622.chrX.pilot_mask.bed -wa -header | bcftools view -Ob > 1kGP_nonPAR_mask.bcf
```

Remove singletons:

```
bcftools view -Ov 1kGP_nonPAR_mask.bcf | vcftools --singletons --vcf -

awk '{if($1 == "chrX")print($1"\t"$2)}' out.singletons > exclude.txt

bcftools view -Ob -T ^exclude.txt 1kGP_nonPAR_mask.bcf > 1kGP_nonPAR_mask_noSing.bcf
```

## Other Data Prep

### Get maf for 1kGP variants used in PD construction

```
bcftools view /net/snowwhite/home/beckandy/research/phasing/data/1kgp/1kGP_nonPAR_mask_noSing.bcf |\
vcftools --vcf - --freq --out 1kGP_nonPAR_mask_noSing
```

## Sample Pairs for Pseudodiploid Construction

The file `sample_X_pairs.R` in the code directory randomly samples pairs of males from the same sub-population. This was modified slightly on 9 Aug 2023 to ensure that each sampled pair was unique. Also, the AMR, EAS, and SAS populations were added (100 pseudodiploids for each population).
