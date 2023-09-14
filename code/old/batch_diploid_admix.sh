#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --time=01:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-200
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/admix_switch_errors/slurm/pair.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing/output/admix_switch_errors/slurm/pair.%A.%a.out


out_dir="/net/snowwhite/home/beckandy/research/phasing/output/admix_switch_errors"
# Read line from sample_pairs_admx_8nov22.csv
line=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing/data/sample_pairs_admx_8nov22.csv)
arrSub=(${line//,/ })
sub_a=${arrSub[0]}
sub_b=${arrSub[1]}

echo "Our subjects are $sub_a and $sub_b"

# Get list of sites for our two subjects
bcftools view -v snps -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf | \
  bcftools view -s $sub_a,$sub_b -c1 -Ou |\
  bcftools query -e 'GT="1|0" || GT="0|1"' -f '%CHROM\t%POS\n' | tr -s ' ' > $out_dir/site_list/sites_${SLURM_ARRAY_TASK_ID}.tsv
echo -e "Site list generated"

# Filter VCF to our subjects and their list of sites
bcftools view -v snps -s $sub_a,$sub_b -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf |\
	bcftools view -T $out_dir/site_list/sites_${SLURM_ARRAY_TASK_ID}.tsv -Ov |\
	bcftools view -c 1 -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Generate pseudodiploid from two haplotypes
Rscript /net/snowwhite/home/beckandy/research/phasing/code/diploidifier.R /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID}

# eagle documentation processing steps
bcftools view /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ou -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | \
  bcftools view -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test2.vcf.gz
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
mv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test2.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
bcftools index /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz

# Save number of heterozygotes
n_hets=$(bcftools query -f '[%GT\t]\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz | awk '{if($1 == "1/0" || $1 == "0/1")count++}END{print(count)}')
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_hets" > $out_dir/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_hets.txt
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_hets"

# Get list of heterozygouse sites
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz |\
awk '{if($3 == "1/0" || $3 == "0/1")print($0)}' > $out_dir/het_loc/pair_${SLURM_ARRAY_TASK_ID}_het_loc.txt

# Reference VCF
bcftools view -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf | \
  bcftools view -s^$sub_a,$sub_b -Ou | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf

n_ref=$(bcftools query -f '%POS\n' /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf | wc -l)
echo -e "${SLURM_ARRAY_TASK_ID}\t$n_ref" > $out_dir/vcf_n_sites/pair_${SLURM_ARRAY_TASK_ID}_ref.txt

### PHASING
eagle --vcfTarget /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  --vcfRef /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --geneticMapFile=/net/snowwhite/home/beckandy/software/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
  --vcfOutFormat v \
  --chrom chrX \
  --numThreads 4 \
  --outPrefix=/tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}

echo "EAGLE done"

/net/snowwhite/home/beckandy/software/bad_shapeit/shapeit4-4.2.2/bin/shapeit4.2 \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  --map /net/snowwhite/home/beckandy/research/phasing/data/shapeit/chrX.b38.gmap.gz \
  --reference /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --region chrX \
  --thread 4 \
  --output /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "SHAPEIT done"

# BEAGLE

bcftools view /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf -Ov > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf

java -Xmx4g -jar /net/snowwhite/home/beckandy/bin/beagle.05May22.33a.jar \
  gt=/tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  ref=/tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf \
  map=/net/snowwhite/home/beckandy/research/phasing/data/ref/plink.chrX.map \
  nthreads=4 \
  impute=false \
  out=/tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}

# Generate switch error files

vcftools --vcf /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out $out_dir/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out $out_dir/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out  $out_dir/switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}

# Count number of heterozygous sites within each region
## Generate "bed"-like file with sites
bcftools query -i 'GT="1/0" || GT="0/1"' -f '%CHROM\t%POS\t%POS\n' /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz  |\
  tr -s ' ' |\
  awk '{print($1"\t"$2-1"\t"$3)}' > /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed

## bedtools intersect
tail -n+2 $out_dir/switch_errors/eagle/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
  bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > $out_dir/switch_errors/eagle/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed

tail -n+2 $out_dir/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
  bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > $out_dir/switch_errors/shapeit/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed

tail -n+2 $out_dir/switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}.diff.switch |\
  bedtools intersect -a - -b /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed -c > $out_dir/switch_errors/beagle/hets_per_region_${SLURM_ARRAY_TASK_ID}.bed

# whatshap adds additional summary statistics regarding phase errors
# need gunzipped truth
gunzip /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz
whatshap compare --names truth,eagle --tsv-pairwise $out_dir/whatshap/eagle/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf
# need gunzipped shapeit
gunzip /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,shapeit --tsv-pairwise $out_dir/whatshap/shapeit/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf

gunzip /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,beagle --tsv-pairwise $out_dir/whatshap/beagle/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf

# stash truth vcf
bcftools view -Ob /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf > $out_dir/vcf/pair_${SLURM_ARRAY_TASK_ID}.bcf

# stash test vcf
bcftools view -Oz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz > $out_dir/vcf/pair_test_${SLURM_ARRAY_TASK_ID}.vcf.gz

# clean up tmp files
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.log
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_eagle_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test*
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf*
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf*
rm /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.bed
