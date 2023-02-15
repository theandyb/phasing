#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=meth
#SBATCH --array=1-300
#SBATCH --requeue
#SBATCH -e /net/snowwhite/home/beckandy/research/phasing/output/19april22_switch_errors/slurm/meth-e-%J.err
#SBATCH -o /net/snowwhite/home/beckandy/research/phasing/output/19april22_switch_errors/slurm/meth-e-%J.out

Rscript /net/snowwhite/home/beckandy/research/phasing/code/add_meth.R /net/snowwhite/home/beckandy/research/phasing/output/19april22_switch_errors/switch_errors/beagle/annotated ${SLURM_ARRAY_TASK_ID}
