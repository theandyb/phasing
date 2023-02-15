from pyfaidx import Fasta
import pandas as pd

ref_file = "/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_mask.fasta"

fasta_obj = Fasta(ref_file)
seq = fasta_obj["chrX"]
seqstr = seq[0:len(seq)].seq
cpg_pos = []

dimer_count = {
  "AA":0,
  "AC":0,
  "AG":0,
  "AT":0,
  "CA":0,
  "CC":0,
  "CG":0,
  "CT":0,
  "GA":0,
  "GC":0,
  "GG":0,
  "GT":0,
  "TA":0,
  "TC":0,
  "TG":0,
  "TT":0
}

for i in range(len(seqstr)-1):
  if(seqstr[i:(i+2)] == "CG"):
    cpg_pos.extend([i, (i+1)])

for i in range(len(seqstr) - 1):
  test = seqstr[i:(i+2)]
  if test in dimer_count:
    dimer_count[test] += 1
    
end_pos = [x + 1 for x in cpg_pos]
chr_list = ["chrX" for x in cpg_pos]

df_dict = {'chrom': chr_list, 'chromStart': cpg_pos, 'chromEnd': end_pos}
df = pd.DataFrame(df_dict)
df.to_csv('/net/snowwhite/home/beckandy/research/phasing/data/cpg_pos_mask.bed', sep="\t", index=False, header=False)
