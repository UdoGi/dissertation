import re

from Bio import SeqIO
from collections import Counter
import pandas as pd
from pathlib import Path


def yield_nuc_counts(p):
    """
    Given a list of paths to FASTA files, return counts of nucleotides
    for each CDS, indexed by CDS and sample
    :param p:
    :return:
    """
    for fn in p:
        with open(fn) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                c = Counter(record.seq)
                # fasta-entry example Diekholzen_2012_N Consensus_ALL_D9591_S_threshold_0_quality_20_extraction
                m = re.match(r"\w*_\d*_(NSs|L|M|N)", record.id)
                assert m is not None, f'unknown record ID: {record.id}'
                yield [f"{m.group(1)}_CDS", " ".join(record.id.split('_')[:-1]),
                       c['A'], c['C'], c['G'], c['T'], c['N'], sum(c.values())]


if __name__ == '__main__':
    with open(snakemake.log.stdout, "w") as stdout_log_file, open(snakemake.log.stderr, "w") as stderr_log_file:

        paths = [Path(p) for p in snakemake.input]
        df = pd.DataFrame(yield_nuc_counts(paths),
                          columns=['CDS', 'case', 'A', 'C', 'G', 'T', 'N', 'total'])
        df['valid_nuc_count'] = df[['A', 'C', 'G', 'T']].sum(axis=1)
        df['perc_complete'] = df.valid_nuc_count/df.total * 100
        df.to_excel(snakemake.output.overall)
        df.pivot(index='case', columns='CDS',
                 values='perc_complete').sort_values('M_CDS')[::-1]\
            .reset_index().to_excel(snakemake.output.cds_completion, index=False)
