import sys
sys.path.append('../scripts')

import pandas as pd
from Bio import SeqIO
from helpers import yield_renamed_records, validate_consensi
from pathlib import Path
import itertools


consensus_dir = config['result_folder'] + 'udo_ref_align_consensus'
sanger_dir = config['result_folder'] + 'udo_ref_align_consensus_sanger'

iqtree_binary = config['iqtree_binary']

used_cases = {case.split('_')[1] for case in cases}

df_bat_meta = pd.read_excel("../meta_data/sampleId_to_caseId.xlsx")
df_bat_meta['sample_id'] = df_bat_meta['sampleId'].astype(str)

dft = df_bat_meta[[e in used_cases for e in df_bat_meta.caseId]]
wanted_fasta_ids = [f"{c.replace(' ', '_')}_{y}_{cds}"
                    for c, y in  dft[['location', 'time of death']].values
                    for cds in ['L', 'M', 'N','NSs']]


rule validate_consensi:
    input:
        consensi = [config['result_folder'] + f'udo_ref_align_consensus/{sample}_udo_ref_align_consensus.fasta'
             for sample in cases],
        sangers = [f'flags/{sample}_consensus_sanger.done' for sample in sangered_samples]
    output:
        flag="flags/consensi_validated.done"
    run:
        if validate_consensi(consensus_dir, sanger_dir, '../CDS_of_consensus', cases):
            Path(output.flag).touch()

rule rename_cds_fasta:
    input:
        consensi_valid="flags/consensi_validated.done",
        cds_fastas=list(Path('../CDS_of_consensus').glob('*consensus.fasta'))
    output:
        renamed_cds=config['result_folder']  + config['consensi']['renamed_CDS_consensus']
    run:
        with open(output.renamed_cds, 'w') as fh:
            SeqIO.write(itertools.chain(*(yield_renamed_records(Path(fn),dft) for fn in input.cds_fastas)),
                fh, 'fasta')


rule record_CDS_completion:
    input:
        renamed_fastas=config['result_folder']  + config['consensi']['renamed_CDS_consensus']
    output:
        overall=config['result_folder']  + 'record_CDS_completion/cds_overall.xlsx',
        cds_completion=config['result_folder']  + 'record_CDS_completion/cds_completion.xlsx'
    log:
        stdout=config['result_folder'] + 'logs/record_CDS_completion/record_CDS_completion_stdout.log',
        stderr=config['result_folder'] + 'logs/record_CDS_completion/record_CDS_completion_stderr.log'
    script:
        "../scripts/record_CDS_completion.py"


rule compute_supp_tables:
    input:
        sample_to_caseId = "../meta_data/sampleId_to_caseId.xlsx",
        bams=[config['result_folder'] + f'udo_ref_align_consensus/{sample}_udo_ref_sorted.bam'
              for sample in cases]
    output:
        xlsx=config['result_folder'] + config['supp']['ngs_details']
    conda:
        '../../envs/biopython_env.yaml'
    params:
        bam_dir=config['result_folder'] + '/udo_ref_align_consensus',
    log:
        stdout=config['result_folder'] + 'logs/compute_supp_tables/stdout.log',
        stderr=config['result_folder'] + 'logs/compute_supp_tables/stderr.log'
    script:
        "../scripts/compute_supp_tables.py"


# MSA for each segment of the ancient+ dataset
rule muscle_alignment_per_segment_ancient_p:
    input:
        fasta = config['result_folder']  + config['consensi']['renamed_CDS_consensus']
    output:
        temp_fasta=temp(config['result_folder']  + 'muscle_align/ancient_p/temp_{cds}.fasta'),
        alignment=config['result_folder']  + 'muscle_align/ancient_p/ancient_p_{cds}_alignment.afa',
    conda:
        '../../envs/muscle_seqkit_env.yaml'
    log:
        stdout=config['result_folder'] + 'logs/muscle_alignment_per_segment_ancient_p/{cds}_stdout.log',
        stderr=config['result_folder'] + 'logs/muscle_alignment_per_segment_ancient_p/{cds}_stderr.log'
    params:
        fasta_ids = lambda wildcards: ','.join([f'"{e}"' for e in wanted_fasta_ids if wildcards.cds == e.split('_')[-1]])
    shell:
        'cat {input.fasta} | seqkit grep -p {params.fasta_ids} > {output.temp_fasta} 2> {log.stderr}; '
        'muscle -align {output.temp_fasta} -output {output.alignment} > {log.stdout} 2>> {log.stderr}; '


# MSA for each segment of the ancient- dataset
rule muscle_alignment_ancient_m:
    input:
        fasta = config['result_folder']  + config['consensi']['renamed_CDS_consensus']
    output:
        temp_fasta=temp(config['result_folder']  + 'muscle_align/ancient_m/temp_L.fasta'),
        alignment=config['result_folder']  + 'muscle_align/ancient_m/ancient_m_L_alignment.afa',
    conda:
        '../../envs/muscle_seqkit_env.yaml'
    log:
        stdout=config['result_folder'] + 'logs/muscle_alignment_ancient_m/stdout.log',
        stderr=config['result_folder'] + 'logs/muscle_alignment_ancient_m/stderr.log'
    params:
        fasta_ids = lambda wildcards: ','.join([f'"{e}"' for e in wanted_fasta_ids
                                                if (("L" == e.split('_')[-1]) and
                                                    ("penzlin" not in e.lower() )) ])
    shell:
        'cat {input.fasta} | seqkit grep -p {params.fasta_ids} > {output.temp_fasta} 2> {log.stderr}; '
        'muscle -align {output.temp_fasta} -output {output.alignment} > {log.stdout} 2>> {log.stderr}; '


# iqtree on each MSA
rule run_iqTree_per_segment:
    input:
        alignment=config['result_folder']  + 'muscle_align/ancient_p/ancient_p_{cds}_alignment.afa',
    output:
        tree = config['result_folder']  + 'iqtrees/ancient_p/ancient_p_{cds}.treefile',
    params:
        out_prefix = lambda wildcards: config['result_folder']  + f'iqtrees/ancient_p/ancient_p_{wildcards.cds}',
        iqtree = iqtree_binary
    log:
        stdout=config['result_folder'] + 'logs/run_iqTree_per_segment/{cds}_stdout.log',
        stderr=config['result_folder'] + 'logs/run_iqTree_per_segment/{cds}_stderr.log'
    shell:
        "{params.iqtree} -s {input.alignment} -pre {params.out_prefix} -bb 10000 -nt AUTO -redo -mtree "
        "> {log.stdout} 2> {log.stderr} "

