import sys
from pathlib import Path
from collections import defaultdict
import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO


def compute_coverage(fn_bam, fn_references):
    with open(fn_references) as fin:
        len_lookup = {fasta_record.id: len(fasta_record.seq)
                      for fasta_record in SeqIO.parse(fin, 'fasta')}

    alignment_lookup = defaultdict(list)
    # load the alignments in a dict st they're ordered per ref
    for alignment in pysam.AlignmentFile(fn_bam, "rb"):
        alignment_lookup[alignment.reference_name].append(alignment)
    # compute coverage for each reference
    r = []
    for k, v in alignment_lookup.items():
        coverage_arr = np.zeros(len_lookup[k])
        for m in v:
            coverage_arr[m.get_reference_positions()] += 1
        breadth = (coverage_arr >= 1).sum() / len(coverage_arr)
        r.append([k, breadth, len(set(m.query_name for m in v))])
    return pd.DataFrame(r, columns=['acc', 'breadth', 'no_reads'])


if __name__ == '__main__':

    with open(snakemake.log.stderr, 'w') as fp_stderr:
        with open(snakemake.log.stdout, 'w') as fp_stdout:
            sys.stderr = fp_stderr
            sys.stdout = fp_stdout

            df = compute_coverage(snakemake.input.bam,
                                  snakemake.input.virus_fasta)
            df['dataset'] = Path(snakemake.input.bam).stem
            df['nc_accession'] = df['acc'].apply(lambda x: x.split('|')[0])
            df['tax_id'] = df['acc'].apply(lambda x: str(x.split('|')[-1]))
            df.to_csv(snakemake.output.csv, index=False)
