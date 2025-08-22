from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam
import pysamstats
import sys

LEN_LOOKUP = {
    b"L_segment_udo_virus": 6350,
    b"M_segment_udo_virus": 3442,
    b"S_segment_udo_virus": 1721,
}


def get_no_mapped_reads(my_bam):
    """
    :param my_bam: The input BAM file from which to count mapped reads.
    :return: A dictionary mapping reference names to the count of mapped reads.
    Given a bam file it counts per reference the number of mapped read ids. In case of paired-end
    data, only the read itself is counted. So a paired-end read only counts as one, no matter if they
    bound concordantly or not.
    """
    seg_count = defaultdict(set)
    for record in my_bam:
        seg_count[record.reference_name].add(record.query_name)
    return {ref: len(s) for ref, s in seg_count.items()}


def yield_aln_dfs(paths):
    """
    Given a list of paths to BAM files containing mapped reads against the tri-segmented bandavirus,
    this fucntion yields a dataframe listing aggregate mapping statistics per segment.
    :param paths:
    :return:
    """
    for p in paths:
        with pysam.AlignmentFile(p) as my_bam:
            df = pd.DataFrame(pysamstats.load_coverage(my_bam))[
                ["chrom", "pos", "reads_all"]
            ]
        with pysam.AlignmentFile(p) as my_bam:
            df["dataset"] = p.stem
            mapped_read_d = get_no_mapped_reads(my_bam)
            df_agg = df.groupby(["chrom", "dataset"]).agg(
                ["min", "max", "count", "mean", "sum"]
            )
            df_agg["mapped_read_count"] = [
                mapped_read_d[a.decode("utf-8")] for (a, b) in df_agg.index
            ]
            yield df_agg


def annotate_df(df):
    """
    Given a dataframe with mapping statistics per segment and dataset, this function computes further statistics
    and renames the columns, collapsing the two-level columns into one level
    :param df: dataframe with mapping statistics per segment and dataset
    :return:
    """
    df.reset_index(inplace=True)
    df["official_len"] = [LEN_LOOKUP[c] for c in df["chrom"]]
    df["my_mean_coverage"] = df[("reads_all", "sum")] / df[("official_len")]
    df["breadth"] = 100 * df[("reads_all", "count")] / df[("official_len")]
    df2 = df[
        [
            ("chrom", ""),
            ("dataset", ""),
            ("pos", "count"),
            ("reads_all", "min"),
            ("reads_all", "max"),
            ("my_mean_coverage", ""),
            ("mapped_read_count", ""),
            ("breadth", ""),
        ]
    ]
    df2.columns = [" ".join(c).strip() for c in df2.columns]
    return df2


def rename_columns(df, case_id_lookup):
    """
    Renames the columns of the dataframe and renames the samples using the passed
    `case_id_lookup` dictionary.
    :param df:
    :param case_id_lookup:
    :return:
    """
    df = df[
        [
            "chrom",
            "dataset",
            "reads_all min",
            "reads_all max",
            "mapped_read_count",
            "my_mean_coverage",
            "breadth",
        ]
    ].copy()
    df.columns = [
        "segment",
        "sample",
        "min coverage",
        "max coverage",
        "total mapped reads",
        "mean coverage",
        "breadth",
    ]
    df["segment"] = [seg.split("_")[0] for seg in df.segment.astype(str)]
    df["sample"] = [case_id_lookup[samp.split("_")[1]] for samp in df["sample"]]
    return df


if __name__ == "__main__":
    sys.stderr = open(snakemake.log.stderr, "w")
    sys.stdout = open(snakemake.log.stdout, "w")

    paths = list(Path(snakemake.params.bam_dir).glob("*.bam"))
    df = annotate_df(pd.concat(yield_aln_dfs(paths)))

    case_id_lookup = {
        case: " ".join([location, str(time)])
        for case, location, time in pd.read_excel(snakemake.input.sample_to_caseId)[
            ["caseId", "location", "time of death"]
        ].values
    }

    rename_columns(df, case_id_lookup).to_excel(snakemake.output.xlsx, index=False)
