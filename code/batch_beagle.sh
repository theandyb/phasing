#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-300
#SBATCH --partition=nomosix
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/19april22_switch_errors/slurm/beagle.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing/output/19april22_switch_errors/slurm/beagle.%A.%a.out

# Read line from sample_pairs.csv
line=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing/data/sample_pairs_all.csv)
arrSub=(${line//,/ })
sub_a=${arrSub[1]}
sub_b=${arrSub[2]}

echo "Our subjects are $sub_a and $sub_b"

start_time=$(date +%s)

# Get list of sites for our two subjects
bcftools view -v snps -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf | \
  bcftools view -s $sub_a,$sub_b -c1 -Ou |\
  bcftools query -e 'GT="1|0" || GT="0|1"' -f '%CHROM\t%POS\n' | tr -s ' ' > /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.tsv
echo -e "Site list generated"

# Filter VCF to our subjects and their list of sites
bcftools view -v snps -s $sub_a,$sub_b -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf |\
	bcftools view -T /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.tsv -Ov |\
	bcftools view -c 1 -Oz > /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Generate pseudodiploid from two haplotypes
Rscript /net/snowwhite/home/beckandy/research/phasing/code/diploidifier.R /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz /tmp/andy_input_${SLURM_ARRAY_TASK_ID}

### REFERENCE
# Filter our subjects out of the reference
bcftools view -c 2 -Ou /net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot.bcf | \
  bcftools view -s^$sub_a,$sub_b -Ou | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Oz -d none -f /net/snowwhite/home/beckandy/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf.gz

java -Xmx4g -jar /net/snowwhite/home/beckandy/bin/beagle.05May22.33a.jar \
  gt=/tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz \
  ref=/tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  map=/net/snowwhite/home/beckandy/research/phasing/data/ref/plink.chrX.GRCh38.map \
  nthreads=2 \
  out=/tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}

echo "BeAGLE DONE"

#/net/snowwhite/home/beckandy/research/phasing/output/14mar22_switch_errors
# Generate switch error files

vcftools --gzvcf /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz \
  --gzdiff /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz \
  --diff-switch-error \
  --out /net/snowwhite/home/beckandy/research/phasing/output/19april22_switch_errors/switch_errors/beagle_nomap/error_${SLURM_ARRAY_TASK_ID}
# perform a task
end_time=$(date +%s)

# elapsed time with second resolution
elapsed=$(( end_time - start_time ))
echo "$elapsed"

# cleanup
rm /tmp/andy_ref_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.tsv
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_test.vcf.gz
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}_truth.vcf.gz
rm /tmp/andy_input_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.vcf.gz
rm /tmp/andy_beagle_${SLURM_ARRAY_TASK_ID}.log

