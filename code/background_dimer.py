from pyfaidx import Fasta
import pandas as pd
ref_file = "/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_mask.fasta"

def count_string(nuc_string, result_table):
  nucs = ["A", "C", "G", "T"]
  for i in range(len(nuc_string)-1):
      nuc1 = nuc_string[i]
      nuc2 = nuc_string[i+1]
      if i % 1000000 == 0:
        print(i)
      if not(nuc1 in nucs):
        continue
      if not(nuc2 in nucs):
        continue
      dimer = nuc1+nuc2
      result_table[dimer] += 1
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
  out_file = out_dir + "dimer.csv"
  fasta_obj = Fasta(ref_file)
  
  seq = fasta_obj["{}{}".format(chrom_prefix, chrom)]
  seqstr = seq[0:len(seq)].seq
  result_table = {"AA": 0, "AC": 0, "AG":0, "AT":0,
  "CA": 0, "CC": 0, "CG":0, "CT":0,
  "GA": 0, "GC": 0, "GG":0, "GT":0,
  "TA": 0, "TC": 0, "TG":0, "TT":0,}
  count_string(seqstr, result_table)
  df = table_df(result_table)
  df.to_csv(out_file, index=False)
  print("Done!")
  return 0

if __name__ == "__main__":
  chrom = "X"
  main(chrom)
