#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --time=01:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-400
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/27sept/slurm/pair.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing/output/27sept/slurm/pair.%A.%a.out

out_dir="/net/snowwhite/home/beckandy/research/phasing/output/27sept"
in_dir="/net/snowwhite/home/beckandy/research/phasing/output/9aug22_switch_errors"
# Read line from sample_pairs.csv
line=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing/data/sample_pairs_16aug2022.csv)
arrSub=(${line//,/ })
sub_a=${arrSub[1]}
sub_b=${arrSub[2]}

echo "Our subjects are $sub_a and $sub_b"

site_file=$in_dir/site_list/sites_${SLURM_ARRAY_TASK_ID}.tsv

# Filter VCF to our subjects and their list of sites
bcftools view -v snps -s $sub_a,$sub_b -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf |\
	bcftools view -T $site_file -Ov |\
	bcftools view -c 1 -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Generate pseudodiploid from two haplotypes
Rscript /net/snowwhite/home/beckandy/research/phasing/code/diploidifier.R /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}

# eagle documentation processing steps
bcftools view /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ou -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | \
  bcftools view -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test2.vcf.gz
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
mv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test2.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
bcftools index /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz

# Reference VCF
bcftools view -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf | \
  bcftools view -s^$sub_a,$sub_b -Ou | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf
bcftools index /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf

#SHAPEIT4 - use phasing option

/net/snowwhite/home/beckandy/software/bad_shapeit/shapeit4-4.2.2/bin/shapeit4.2 \
  --input /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  --map /net/snowwhite/home/beckandy/research/phasing/data/shapeit/chrX.b38.gmap.gz \
  --reference /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf \
  --region chrX \
  --thread 4 \
  --sequencing \
  --output /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "SHAPEIT done"

# BEAGLE - no genetic map

bcftools view /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf -Ov > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf

java -Xmx4g -jar /net/snowwhite/home/beckandy/bin/beagle.05May22.33a.jar \
  gt=/tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  ref=/tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf \
  nthreads=4 \
  out=/tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}

# Generate switch error files

vcftools --gzvcf /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out $out_dir/switch_errors/shapeit/error_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out  $out_dir/switch_errors/beagle/error_${SLURM_ARRAY_TASK_ID}

# whatshap adds additional summary statistics regarding phase errors
# need gunzipped truth
gunzip /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz
gunzip /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,shapeit --tsv-pairwise $out_dir/whatshap/shapeit/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf

gunzip /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
whatshap compare --names truth,beagle --tsv-pairwise $out_dir/whatshap/beagle/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf  /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf


# clean up tmp files
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.log
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test*
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.bcf*
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf*
rm /tmp/andy_shapeit_${SLURM_ARRAY_TASK_ID}.vcf
