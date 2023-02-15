#!/bin/bash
#
#SBATCH --job-name=countHets
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-1233
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/X_hets/slurm/count.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing/output/X_hets/slurm/count.%A_%a.out

sub=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'NR==id{ print; exit }' /net/snowwhite/home/beckandy/research/phasing/data/1kgp/male_ids.txt)
vcf="/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_nonPAR_pilot_red.bcf"
cVcf="/net/snowwhite/home/beckandy/research/phasing/data/ref/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
outdir="/net/snowwhite/home/beckandy/research/phasing/output/X_hets"

# Get list of het sites
bcftools view -s ${sub} -v snps $vcf |\
  bcftools query -e 'GT="hom"' -f '%POS\n' > /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.txt
sed -i -e 's/^/chrX\t/' /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.txt

n_hets=$(wc -l < /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.txt)
echo $n_hets > $outdir/nhet/${SLURM_ARRAY_TASK_ID}.txt

# Get GT call from cVcf
bcftools view -s $sub -R /tmp/andy_sites_1.txt $cVcf |\
  bcftools query -f'[%GT\n]' |\
  awk '{count[$0]++}END{for(key in count)print(key"\t"count[key])}' >\
  $outdir/call_count/${SLURM_ARRAY_TASK_ID}.tsv

# clean up
rm /tmp/andy_sites_${SLURM_ARRAY_TASK_ID}.txt
