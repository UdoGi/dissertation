import sys
sys.path.append('../scripts')
from pathlib import Path
import re
import pandas as pd

"""
based on the csv-meta data, merge the raw fastqs.
Putting all runs from all organs into one file per case.
"""


# df = pd.read_excel("../meta_data/all_paper_bat_info.xlsx", engine='openpyxl')
# fastq_dir = Path(config['fastq']['raw_fastq_dir'])
# cases = df.caseId.unique()

# rule all:
#     input:
#         excel = config['result_folder'] + config['fastq']['fastq_merged_meta_fn']
#

def get_merge_input_files(wildcards):
    """
    Based on the case-id, all associated sequencing runs are searched for in the
    provided fastq-directory and returned.
    This assumes that each sequencing run is paired-end
    :param wildcards:
    :return:
    """
    df_fn = df_raw[df_raw.caseId == wildcards.case]
    p_dir = Path(fastq_dir)
    r1s, r2s = [],[]
    for snake_name in df_fn.sample_name_snake:
        for p in p_dir.glob(f"{snake_name}*"):
            pat = f".*{snake_name}.*_(R1|R2).*\\.fastq\\.gz"
            s = str(p)
            m = re.match(pat,s)
            if m:
                r1s.append(p) if m.group(1) == 'R1' else r2s.append(p)
    assert len(r1s) == len(r2s), 'unequal number of R1 and R2 datasets!'
    assert len(r1s) > 0, f'No FASTQs found in the specified location {fastq_dir}'
    return {'R1' : r1s, 'R2' : r2s}


rule create_merged_files:
    input:
        unpack(get_merge_input_files)
    output:
        R1 = fastq_dir /  "ALL_{case}_R1_001.fastq.gz",
        R2 = fastq_dir / "ALL_{case}_R2_001.fastq.gz"
    threads:
        2
    resources:
        mem_mb='6000M',
        disk_mb='10000M',
        time="01:00:00"
    shell:
        "cat {input.R1} > {output.R1};"
        "cat {input.R2} > {output.R2};"

rule add_to_meta_data:
    input:
        meta_data_fn = "../meta_data/all_paper_bat_info.xlsx",
        R1 = [fastq_dir / f"ALL_{case}_R1_001.fastq.gz" for case in cases],
        R2 = [fastq_dir / f"ALL_{case}_R2_001.fastq.gz" for case in cases]
    output:
        excel=config['result_folder'] + config['fastq']['fastq_merged_meta_fn']
    threads:
        2
    resources:
        mem_mb='6000M',
        disk_mb='10000M',
        time="01:00:00"
    run:
        l = []
        df2 = df_raw.set_index('caseId')[['description', 'sample_id']].drop_duplicates()
        for r1 in input.R1:
            case = r1.split('ALL_')[1].split('_R1')[0]
            dft = df2.loc[[case]]
            assert len(dft) == 1, f'there is more than one entry for caseId: {case}'
            l.append(list(dft.values[0]) + [f'ALL_{case}', case ])
        df_merge = pd.DataFrame(l, columns = ['description', 'sample_id',
                                        'sample_name_snake', 'caseId'])
        df_merge.to_excel(output.excel, index=False,engine='openpyxl')



