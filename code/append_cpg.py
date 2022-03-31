from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
import argparse

def get_motif(seqstr, pos, bp = 10):
    return seqstr[(pos - bp):(pos + bp + 1)]

def cpg_stat(seqstr, pos):
    motif = get_motif(seqstr, pos, 1)
    if 'CG' in motif:
        ret_val = 1
    else:
        ret_val = 0
    return ret_val
  
parser = argparse.ArgumentParser(description="Annotate genomic locations with GC status")
parser.add_argument("-c", "--chrom", help="Which chromosome?", required=True)
parser.add_argument("-s", "--switch", help="Location of switch file", required=True)
parser.add_argument("-o", "--output", help="Location of output file", required=True)
args = parser.parse_args()

chrom = args.chrom
input_file = args.switch
output_file = args.output

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
        data = line.strip().split("\t") # CHR POS_START POS_END INDV
        # insert program here
        pos_start = int(data[1])
        pos_end = int(data[2])
        start_cpg = cpg_stat(seqstr, pos_start)
        end_cpg = cpg_stat(seqstr, pos_end)
        start_motif = get_motif(seqstr, pos_start, 1)
        end_motif = get_motif(seqstr, pos_end, 1)
        entry = {
            'pos_start' : pos_start,
            'pos_end' : pos_end,
            'cpg_start' : start_cpg,
            'cpg_end' : end_cpg,
            'motif_start' : start_motif,
            'end_motif' : end_motif
        }
        output_list.append(entry)
        line = fp.readline()

pd.DataFrame(output_list).to_csv(output_file, index = None, header=True)
