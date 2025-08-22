import sys
sys.path.append('../scripts')

from helpers import merge_refs_for_overview_tree

used_cases = {
                'D9595', 'D9598', 'D10643',
               'D10376', 'D10116',
               'D10118',
                'D6879',
                'D9591', 'D10758', 'D9593'
               }

acc_list = [
    "NC_021242",
    "NC_078070",
    "NC_027140",
    "JX961619",
    "NC_022630",
    "KF186494.1",
    "KF186497.1",
    "KR017835.1",
    "NC_043450.1",
    "MN509914.2",
    "PP580190.1",
    "NC_043611.1",
    "MZ440342.1",
    "NC_078064.1",
    "NC_078067.1",
    "NC_027717.1",
    "NC_014397.1", # outgroup rift valley fever
    "MN823639.1" # Zwiesel
    # PQ015116.1, Toscana virus
]

consensus_dir = config['result_folder'] + 'udo_ref_align_consensus'
sanger_dir = config['result_folder'] + 'udo_ref_align_consensus_sanger'

iqtree_binary = config['iqtree_binary']


# df = pd.read_excel(config['result_folder'] + config['fastq']['fastq_merged_meta_fn']
#     , engine='openpyxl')
# sangered_samples = ['ALL_D10758']

# rule all:
#     input:
#         config['result_folder'] + "create_overview_tree/full_alignment.treefile",
#         config['result_folder'] + "create_overview_tree/full_alignment_with_alrt.treefile",
#         config['result_folder'] + "create_overview_tree/full_alignment_NAMES_mod.treefile",
#         config['result_folder'] + "create_overview_tree/full_alignment_with_alrt_NAMES_mod.treefile"
#

rule download_refs:
    output:
        fasta=config['result_folder'] + "create_overview_tree/ext_refs/{acc}.fasta"
    log:
        stdout='logs/download_overview_references/{acc}.stdout.log',
        stderr='logs/download_overview_references/{acc}.stderr.log'
    params:
        accession = lambda wildcards: wildcards.acc
    retries: 10
    shell:
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\={params.accession}\&rettype\=fasta\&retmode\=text -O {output.fasta}"
        "> {log.stdout} 2> {log.stderr};"

rule merge_refs:
    input:
        consensi=[config['result_folder'] + f'udo_ref_align_consensus/{sample}_udo_ref_align_consensus.fasta'
                        for sample in cases],
        sangers=[config['result_folder'] +
                 f'udo_ref_align_consensus_sanger/{sample}_udo_ref_align_consensus_sanger.fasta'
                 for sample in sangered_samples],
        ext_refs=[config['result_folder'] + f"create_overview_tree/ext_refs/{acc}.fasta"
                  for acc in acc_list]
    output:
        fastas=config['result_folder'] + "create_overview_tree/temp_all.fasta"
    params:
        consensus_dir=config['result_folder'] + 'udo_ref_align_consensus',
        sanger_dir=config['result_folder'] + 'udo_ref_align_consensus_sanger'
    run:
        merge_refs_for_overview_tree(input.ext_refs, params.consensus_dir,
            params.sanger_dir, used_cases, output.fastas)

rule make_alignment:
    input:
        fastas=config['result_folder'] + "create_overview_tree/temp_all.fasta"
    output:
        alignment = config['result_folder'] + "create_overview_tree/full_alignment.afa"
    conda:
        '../envs/muscle_seqkit_env.yaml'
    log:
        stdout=config['result_folder'] + 'logs/create_overview_tree/make_alignment_stdout.log',
        stderr=config['result_folder'] + 'logs/create_overview_tree/make_alignment_stderr.log'
    shell:
        'muscle -align {input.fastas} -output {output.alignment} > {log.stdout} 2>> {log.stderr}; '



# iqtree on the MSA
rule run_iqtree:
    input:
        alignment = config['result_folder'] + "create_overview_tree/full_alignment.afa"
    output:
        tree = config['result_folder'] + "create_overview_tree/full_alignment.treefile",
        iqtree_log = config['result_folder'] + "create_overview_tree/full_alignment.iqtree"
    params:
        out_prefix = lambda wildcards: config['result_folder']  + "/create_overview_tree/full_alignment",
        iqtree = iqtree_binary
    log:
        stdout=config['result_folder'] + 'logs/create_overview_tree/iqtree_stdout.log',
        stderr=config['result_folder'] + 'logs/create_overview_tree/iqtree_stderr.log'
    shell:
        "{params.iqtree} -s {input.alignment} -pre {params.out_prefix} -bb 10000 -nt AUTO -redo -mtree "
        "> {log.stdout} 2> {log.stderr} "

rule rename_tip_labels:
    input:
        tree = config['result_folder'] + "create_overview_tree/full_alignment.treefile"
    output:
        tree = config['result_folder'] + "create_overview_tree/full_alignment_NAMES_mod.treefile"
    log:
        stdout=config['result_folder'] + 'logs/create_overview_tree/rename_tip_labels_stdout.log',
        stderr=config['result_folder'] + 'logs/create_overview_tree/rename_tip_labels_stderr.log'
    script:
        "../scripts/rename_tip_lables_overview_tree.py"

# iqtree on the MSA
rule run_iqtree_with_alrt:
    input:
        alignment = config['result_folder'] + "create_overview_tree/full_alignment.afa",
        iqtree_log= config['result_folder'] + "create_overview_tree/full_alignment.iqtree"
    output:
        tree = config['result_folder'] + "create_overview_tree/full_alignment_with_alrt.treefile"
    params:
        out_prefix = lambda wildcards: config['result_folder']  + "/create_overview_tree/full_alignment_with_alrt",
        iqtree = iqtree_binary
    log:
        stdout=config['result_folder'] + 'logs/create_overview_tree/iqtree_with_alrt_stdout.log',
        stderr=config['result_folder'] + 'logs/create_overview_tree/iqtree_with_alrt_stderr.log'
    shell:
        """
        model=$(grep "Best-fit model according to BIC" {input.iqtree_log} | awk -F': ' '{{print $2}}') > {log.stdout} 2> {log.stderr};
        {params.iqtree} -s {input.alignment} -pre {params.out_prefix} -bb 10000 -alrt 10000 -m \"$model\" -nt AUTO -redo -mtree >> {log.stdout} 2>> {log.stderr};
        """

rule rename_tip_labels_alrt:
    input:
        tree = config['result_folder'] + "create_overview_tree/full_alignment_with_alrt.treefile"
    output:
        tree = config['result_folder'] + "create_overview_tree/full_alignment_with_alrt_NAMES_mod.treefile"
    log:
        stdout=config['result_folder'] + 'logs/create_overview_tree/rename_tip_labels_alrt_stdout.log',
        stderr=config['result_folder'] + 'logs/create_overview_tree/rename_tip_labels_alrt_stderr.log'
    script:
        "../scripts/rename_tip_lables_overview_tree.py"