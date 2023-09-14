#!/bin/bash
#
#SBATCH --job-name=tmErr
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --array=1-400
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/filter_switch_errors/slurm/topmed.%A_%a.err
#SBATCH --output=/net/snowwhite/home/beckandy/research/phasing/output/filter_switch_errors/slurm/topmed.%A.%a.out

in_dir="/net/snowwhite/home/beckandy/research/phasing/output/final_switch_errors/vcf"
what_dir="/net/snowwhite/home/beckandy/research/phasing/output/final_switch_errors/whatshap/topmed"
vt_dir="/net/snowwhite/home/beckandy/research/phasing/output/final_switch_errors/switch_errors/topmed"

bcftools view -Ov ${in_dir}/pair_${SLURM_ARRAY_TASK_ID}_true.bcf > /tmp/andy_true_${SLURM_ARRAY_TASK_ID}
bcftools view -s "FAKE${SLURM_ARRAY_TASK_ID}" -Ov ${in_dir}/topmed_res/topmed_all.bcf > /tmp/andy_topmed_${SLURM_ARRAY_TASK_ID}

whatshap compare --names truth,topmed --tsv-pairwise ${what_dir}/eval_${SLURM_ARRAY_TASK_ID}.tsv /tmp/andy_true_${SLURM_ARRAY_TASK_ID} /tmp/andy_topmed_${SLURM_ARRAY_TASK_ID}

vcftools --vcf /tmp/andy_topmed_${SLURM_ARRAY_TASK_ID} \
  --diff /tmp/andy_true_${SLURM_ARRAY_TASK_ID} \
  --diff-switch-error \
  --out $vt_dir/error_${SLURM_ARRAY_TASK_ID}

rm /tmp/andy_true_${SLURM_ARRAY_TASK_ID}
rm /tmp/andy_topmed_${SLURM_ARRAY_TASK_ID}
