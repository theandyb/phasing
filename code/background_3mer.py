from pyfaidx import Fasta
import pandas as pd
import itertools
ref_file = "/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_mask.fasta"

def count_string(nuc_string, result_table):
  nucs = ["A", "C", "G", "T"]
  for i in range(1, len(nuc_string)-2):
      nuc1 = nuc_string[i-1]
      nuc2 = nuc_string[i]
      nuc3 = nuc_string[i+1]
      if i % 1000000 == 0:
        print(i)
      if not(nuc1 in nucs):
        continue
      if not(nuc2 in nucs):
        continue
      if not(nuc3 in nucs):
        continue
      motif = nuc1+nuc2+nuc3
      result_table[motif] += 1
  return 0

def table_df(table):
  df = pd.DataFrame.from_dict(table, orient='index')
  df.index.name = 'Nucs'
  df.reset_index(inplace=True)
  df = df.rename({0:'N'}, axis = 1)
  return df

def main(chrom, chrom_prefix = "chr"):
  ref_file = "/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_mask.fasta"
  out_dir = "/net/snowwhite/home/beckandy/research/phasing/output/background_rates/"
  out_file = out_dir + "3mer.csv"
  fasta_obj = Fasta(ref_file)
  
  seq = fasta_obj["{}{}".format(chrom_prefix, chrom)]
  seqstr = seq[0:len(seq)].seq
  
  result_table = {''.join(key):0 for key in itertools.product(['A', 'C', 'G', 'T'], repeat = 3)}
  
  count_string(seqstr, result_table)
  df = table_df(result_table)
  df.to_csv(out_file, index=False)
  print("Done!")
  return 0

if __name__ == "__main__":
  chrom = "X"
  main(chrom)
