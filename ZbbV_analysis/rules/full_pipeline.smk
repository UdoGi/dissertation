import pandas as pd
from pathlib import Path

cases = {
                'ALL_D9595', 'ALL_D9598', 'ALL_D10643',
               'ALL_D10376', 'ALL_D10116',
                'ALL_D10118',
                'ALL_D6879',
                'ALL_D9591', 'ALL_D10758', 'ALL_D9593'
               }

df_raw = pd.read_excel("../meta_data/bat_sample_info.xlsx", engine='openpyxl')

fastq_dir = Path(config['fastq']['raw_fastq_dir'])
sangered_samples = ['ALL_D10758']

include: "merge_fastq_per_case.smk"
include: "preparations.smk"
include: "create_overview_tree.smk"
include: "create_ancient_reference.smk"
include: "create_recent_consensus.smk"
include: "add_sanger_data.smk"
include: "rename_consensi_MSA.smk"

rule all:
    input:
        # consensus sequence for all recovered genomes
        [config['result_folder'] + f'udo_ref_align_consensus/{sample}_udo_ref_align_consensus.fasta'
             for sample in cases],

        # consensus sequences, that are merged with sangered reads
        [ f'flags/{sample}_consensus_sanger.done' for sample in sangered_samples],

        # plot of IFN induction experiments
        config['result_folder'] + config['supp']['ifn_luciferase_fig4'],

        # IQTree for all bandavirus species, including my recovered genomes
        config['result_folder'] + "create_overview_tree/full_alignment_with_alrt_NAMES_mod.treefile",
        config['result_folder'] + "create_overview_tree/full_alignment_NAMES_mod.treefile",

        # table summarizing completion of the CDSs
        config['result_folder'] + 'record_CDS_completion/cds_completion.xlsx',

        # table summarizing NGS statistics
        config['result_folder'] + config['supp']['ngs_details'],

        # alignment for all the CDS
        [config['result_folder']  + f'muscle_align/ancient_p/ancient_p_{cds}_alignment.afa'
            for cds in ['L', 'M', 'N', 'NSs']],


