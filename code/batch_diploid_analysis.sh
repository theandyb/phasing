#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-100
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/slurm/pair.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/slurm/pair.%A.%a.out

# Read line from sample_pairs.csv
line=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing/data/sample_pairs.csv)
arrSub=(${line//,/ })
sub_a=${arrSub[1]}
sub_b=${arrSub[2]}

echo "Our subjects are $sub_a and $sub_b"

start_time=$(date +%s)

# Get list of sites for our two subjects
bcftools view -v snps -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR.bcf | \
  bcftools view -s $sub_a,$sub_b -c1 -Ou |\
  bcftools query -e 'GT="1|0" || GT="0|1"' -f '%CHROM\t%POS\n' | tr -s ' ' > /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.tsv
echo -e "Site list generated"

# Filter VCF to our subjects and their list of sites
bcftools view -v snps -s $sub_a,$sub_b -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR.bcf |\
	bcftools view -T /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.tsv -Ov |\
	bcftools view -c 1 -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Number of sites
n_sites=$(bcftools query -f '%POS\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz | wc -l)
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_sites" > /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_all.txt

# Generate pseudodiploid from two haplotypes
Rscript /net/snowwhite/home/beckandy/research/phasing/code/diploidifier.R /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}

# eagle documentation processing steps
bcftools view /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf
bcftools index /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf

# Save number of heterozygotes
n_hets=$(bcftools query -f '[%GT\t]\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf | awk '{if($1 == "1/0" || $1 == "0/1")count++}END{print(count)}')
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_hets" > /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_hets.txt

### REFERENCE
# Filter our subjects out of the reference
bcftools view -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR.bcf | \
  bcftools view -s^$sub_a,$sub_b -Ou | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf

n_ref=$(bcftools query -f '%POS\n' /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf | wc -l)
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_ref" > /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_ref.txt

### PHASING
eagle --vcfTarget /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf \
  --vcfRef /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --geneticMapFile=/net/snowwhite/home/beckandy/software/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
  --vcfOutFormat v \
  --chrom chrX \
  --outPrefix=/tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}

echo "EAGLE done"

/net/snowwhite/home/beckandy/software/bad_shapeit/shapeit4-4.2.2/bin/shapeit4.2 \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf \
  --map /net/snowwhite/home/beckandy/research/phasing/data/shapeit/chrX.b38.gmap.gz \
  --reference /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --region chrX \
  --thread 2 \
  --output /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz

#/net/snowwhite/home/beckandy/research/phasing/output/14mar22_switch_errors
# Generate switch error files

vcftools --vcf /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}

# Count number of heterozygous sites within each region
## Generate "bed"-like file with sites
bcftools query -i 'GT="1/0" || GT="0/1"' -f '%CHROM\t%POS\t%POS\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf  |\
  tr -s ' ' > /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed

## bedtools intersect
tail -n+2 /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
  bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/switch_errors/eagle/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed

tail -n+2 /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
  bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/switch_errors/shapeit/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed

# whatshap adds additional summary statistics regarding phase errors
# need gunzipped truth
gunzip /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz
whatshap compare --names truth,eagle --tsv-pairwise /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/whatshap/eagle/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf
# need gunzipped shapeit
gunzip /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,shapeit --tsv-pairwise /net/snowwhite/home/beckandy/research/phasing/output/18mar22_switch_errors/whatshap/shapeit/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf


# perform a task
end_time=$(date +%s)

# elapsed time with second resolution
elapsed=$(( end_time - start_time ))
echo "$elapsed"

# cleanup
rm /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.tsv
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.bcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed
