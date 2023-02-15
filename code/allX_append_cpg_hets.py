from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
import argparse

def get_motif(seqstr, pos, bp = 10):
    return seqstr[(pos - bp - 1):(pos + bp)]

def cpg_stat(seqstr, pos):
    motif = get_motif(seqstr, pos, 1)
    if 'CG' in motif:
        ret_val = 1
    else:
        ret_val = 0
    return ret_val
  
chrom = "X"
input_file = "/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_hets.tsv"
output_file = "/net/snowwhite/home/beckandy/research/phasing/data/ref/chrX_hets_anno.csv"

ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
print("Reading reference file...")
fasta_obj = Fasta(ref_file)
seq = fasta_obj["chr{}".format(chrom)]
seqstr = seq[0:len(seq)].seq
print("Reference read!")
output_list = []

with open(input_file) as fp:
    line = fp.readline() # header
    line = fp.readline() # first line of data
    while line:
        data = line.strip().split(" ") # CHR POS MAF
        # insert program here
        pos_start = int(data[1])
        start_cpg = cpg_stat(seqstr, pos_start)
        entry = {
            'pos_start' : pos_start,
            'cpg_start' : start_cpg
        }
        output_list.append(entry)
        line = fp.readline()

pd.DataFrame(output_list).to_csv(output_file, index = None, header=True)
