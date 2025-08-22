
import sys
sys.path.append('../scripts')

"""
This pipeline creates the necessary files to create the consensus sequences for each sample,
maps against oxidase-c for species confirmation and creates the figures for IFN induction experiments.

First it builds a ribosomal RNA consensus sequence for the Pipistrellus pipistrellus species.
For 18s and 28s there is only a predicted sequence available on NCBI for a different bat species.
Hence i create a consensus sequence by using one of the samples with a high rRNA content.
The rRNA reference is used to remove rRNA reads from the ancient sequence data, such that
a de-novo assembly is feasible using reasonable resources (time/memory). 

Additionally the reads of all samples are filtered with fastp for minimum length, deduplicated and overlapping
paired and reads merged into one large read. 
"""


# df = pd.read_excel(config['fastq']['fastq_merged_meta_fn'], engine='openpyxl')
#
# rule all:
#     input:
#         [f"{config['result_folder']}fastp_paired/{sample}_R1.fastp.fastq.gz"
#          for sample in df.sample_name_snake.unique() ],
#         [config['result_folder'] + f'map_vs_oxidase_bats/{sample}_oxidase_bats_ref.bam'
#                     for sample in df.sample_name_snake.unique()],
#         config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_singletons.fastp.fastq.gz",
#         fn_fig4=config['result_folder'] + config['ifn_luciferase_fig4'],


# fastp
module fastp_workflow:
    snakefile: "fastp.smk"
    config: config

module bowtie_workflow:
    snakefile: "bowtie.smk"
    config: config


use rule fastp_paired from fastp_workflow with:
    input:
        R1=config['fastq']['raw_fastq_dir'] + "{sample}_R1_001.fastq.gz",
        R2=config['fastq']['raw_fastq_dir'] + "{sample}_R2_001.fastq.gz"
    params:
        min_length=30

rule download_rRNA_mitochondria_references:
    output:
        fasta_12s_and_16s=config['result_folder'] + config['rRNA']['fasta_12s_and_16s'],
        fasta_18s=config['result_folder'] + config['rRNA']['fasta_template_18s'],
        fasta_28s=config['result_folder'] + config['rRNA']['fasta_template_28s'],
        fasta_mito=config['result_folder'] + config['anti_map']['fasta_mitochondria']
    log:
        stdout='logs/download_rRNA_references.stdout.log',
        stderr='logs/download_rRNA_references.stderr.log'
    retries: 10
    shell:
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=AF326105.1\&rettype\=fasta\&retmode\=text -O {output.fasta_12s_and_16s}"
        "> {log.stdout} 2> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=XR_002778881.1\&rettype\=fasta\&retmode\=text -O {output.fasta_18s}"
        ">> {log.stdout} 2>> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=XR_008558677.1\&rettype\=fasta\&retmode\=text -O {output.fasta_28s}"
        ">> {log.stdout} 2>> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=LR862378.1\&rettype\=fasta\&retmode\=text -O {output.fasta_mito}"
        ">> {log.stdout} 2>> {log.stderr};"



use rule run_bowtie2_paired_RNA  from bowtie_workflow as run_bowtie2_18s with:
    input:
        R1=config['result_folder'] + "fastp_paired/ALL_D10116_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/ALL_D10116_R2.fastp.fastq.gz",
        singles=config['result_folder'] + "fastp_paired/ALL_D10116_singletons.fastp.fastq.gz",
        reference_fn=config['result_folder'] + config['rRNA']['fasta_template_18s'],
    output:
        bam=config['result_folder'] + 'run_bowtie2_rRNA/ALL_D10116_18s.bam'
    log:
        stdout='logs/run_bowtie2_18s.stdout.log',
        stderr='logs/run_bowtie2_18s.stderr.log'
    params:
        single=False,
        reference_file=True,
        bowtie2Args="--sensitive --no-unal --no-discordant"


use rule run_bowtie2_paired_RNA  from bowtie_workflow as run_bowtie2_28s with:
    input:
        R1=config['result_folder'] + "fastp_paired/ALL_D10116_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/ALL_D10116_R2.fastp.fastq.gz",
        singles=config['result_folder'] + "fastp_paired/ALL_D10116_singletons.fastp.fastq.gz",
        reference_fn=config['result_folder'] + config['rRNA']['fasta_template_28s'],
    output:
        bam=config['result_folder'] + 'run_bowtie2_rRNA/ALL_D10116_28s.bam'
    log:
        stdout='logs/run_bowtie2_28s.stdout.log',
        stderr='logs/run_bowtie2_28s.stderr.log'
    params:
        single=False,
        reference_file=True,
        bowtie2Args="--sensitive --no-unal --no-discordant"

rule create_rRNA_consensus:
    input:
        bam_18s=config['result_folder'] + 'run_bowtie2_rRNA/ALL_D10116_18s.bam',
        bam_28s=config['result_folder'] + 'run_bowtie2_rRNA/ALL_D10116_28s.bam',
        fasta_12s_and_16s=config['result_folder'] + config['rRNA']['fasta_12s_and_16s'],
    output:
        cons_18s=config['result_folder'] + config['rRNA']['fasta_18s'],
        sorted_bam_18s=temp(config['result_folder'] + 'rRNA_consensus/rRNA_18s.sorted.bam'),
        cons_28s=config['result_folder'] + config['rRNA']['fasta_28s'],
        sorted_bam_28s=temp(config['result_folder'] + 'rRNA_consensus/rRNA_28s.sorted.bam'),
    log:
        stdout='logs/rRNA_consensus/stdout.log',
        stderr='logs/rRNA_consensus/stderr.log'
    conda:
        "../../envs/alignment_consensus_env.yaml"
    threads:
        16
    resources:
        mem_mb='32000M',
        disk_mb='10000M',
        time="04:00:00"
    params:
        prefix_18s=lambda wildcards: config['result_folder'] + config['rRNA']['fasta_18s'].split('.fa')[0],
        prefix_28s= lambda wildcards: config['result_folder'] + config['rRNA']['fasta_28s'].split('.fa')[0]
    shell:
        "samtools sort {input.bam_18s} -o {output.sorted_bam_18s} --threads {threads} > {log.stdout} 2> {log.stderr}; "
        "samtools mpileup -aa -A -d 0 -Q 0 {output.sorted_bam_18s} | ivar consensus -p {params.prefix_18s} >> {log.stdout} 2>> {log.stderr}; "
        "samtools sort {input.bam_28s} -o {output.sorted_bam_28s} --threads {threads} > {log.stdout} 2> {log.stderr}; "
        "samtools mpileup -aa -A -d 0 -Q 0 {output.sorted_bam_28s} | ivar consensus -p {params.prefix_28s} >> {log.stdout} 2>> {log.stderr}; "


# make bowtie2-DB to anti-map against rRNA and mitochondria to remove unwanted reads
rule make_anti_map_bowtie_db:
    input:
        fasta_12s_16s=config['result_folder'] + config['rRNA']['fasta_12s_and_16s'],
        fasta_18s=config['result_folder'] + config['rRNA']['fasta_18s'],
        fasta_28s=config['result_folder'] + config['rRNA']['fasta_28s'],
        fasta_mito=config['result_folder'] + config['anti_map']['fasta_mitochondria']
    output:
        all_anti_refs=temp('/tmp/make_rRNA_bowtie_db/all_anti_map_refs.fasta'),
        db=f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_anti_ref']}/{config['bowtie2']['bowtie2_db_name_anti']}.1.bt2"
    log:
        stdout = config['result_folder'] + 'logs/make_rRNA_bowtie_db/stdout.log',
        stderr = config['result_folder'] + 'logs/make_rRNA_bowtie_db/stderr.log',
    conda:
        "../../envs/alignment_consensus_env.yaml"
    threads:
        32
    params:
        db_dir=config['result_folder'] + config['bowtie2']["bowtie2_db_folder_anti_ref"],
        db_name=config['bowtie2']["bowtie2_db_name_anti"]
    shell:
        "rm -rf {params.db_dir} > {log.stdout} 2> {log.stderr} ; mkdir {params.db_dir} >> {log.stdout} 2>> {log.stderr};"
        "cd {params.db_dir}; "
        "cat {input} > {output.all_anti_refs};"
        "bowtie2-build {output.all_anti_refs} {params.db_name} --threads {threads} >> {log.stdout} 2>> {log.stderr}; "


use rule make_bowtie2_db_RNA from bowtie_workflow as make_bat_oxidase_bowtie_db with:
    input:
        fasta=config['result_folder'] + config['reference_seq']['oxidase_bats']
    output:
        f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_oxidase_ref']}/{config['bowtie2']['bowtie2_db_name_oxidase']}.1.bt2"
    log:
        stdout = config['result_folder'] + 'logs/make_bat_oxidase_bowtie_db/stdout.log',
        stderr= config['result_folder'] + 'logs/make_bat_oxidase_bowtie_db/stderr.log',
    params:
        db_dir=config['result_folder'] + config['bowtie2']["bowtie2_db_folder_oxidase_ref"],
        db_name=config['bowtie2']["bowtie2_db_name_oxidase"]

# anti map reads of the ancient sample to reduce the number of reads
# that are used for de-novo assembly. This is only done for the ancient sample
# as the others are not de-novo assembled.
rule run_bowtie2_anti_map_paired:
    input:
        db=f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_anti_ref']}/{config['bowtie2']['bowtie2_db_name_anti']}.1.bt2",
        R1=config['result_folder'] + "fastp_paired/ALL_D6879_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/ALL_D6879_R2.fastp.fastq.gz",
        singles=config['result_folder'] + "fastp_paired/ALL_D6879_singletons.fastp.fastq.gz"
    output:
        R1_out=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_R1.fastp.fastq.gz",
        R2_out=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_R2.fastp.fastq.gz",
        single_out=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_singletons.fastp.fastq.gz",
    log:
        stdout="logs/bowtie2_search/run_bowtie2_anti_map_paired/ALL_D6879.stdout.log",
        stderr="logs/bowtie2_search/run_bowtie2_anti_map_paired/ALL_D6879.stderr.log"
    conda:
        '../../envs/alignment_consensus_env.yaml'
    resources:
        mem_mb='64000M',
        disk_mb='10000M',
        time="04:00:00"
    params:
        db=f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_anti_ref']}/{config['bowtie2']['bowtie2_db_name_anti']}",
        paired_prefix= f"{config['result_folder']}fastp_paired/ALL_D6879_anti_map_R%.fastp.fastq.gz"
    threads:
        8
    shell:
        "bowtie2 -x {params.db} -U {input.singles} -p {threads} --un-gz {output.single_out} "
        "-S /dev/null > {log.stdout} 2> {log.stderr}; "
        "bowtie2 -x {params.db} -1 {input.R1} -2 {input.R2} -p {threads} --un-conc-gz {params.paired_prefix} "
        "-S /dev/null >> {log.stdout} 2>> {log.stderr}; "


use rule run_bowtie2_paired_RNA from bowtie_workflow as map_vs_oxidase_bats with:
    input:
        R1=config['result_folder'] + "fastp_paired/{sample}_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/{sample}_R2.fastp.fastq.gz",
        singles=config['result_folder'] + 'fastp_paired/{sample}_singletons.fastp.fastq.gz',
        index=f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_oxidase_ref']}/{config['bowtie2']['bowtie2_db_name_oxidase']}.1.bt2"
    output:
        bam = config['result_folder'] + 'map_vs_oxidase_bats/{sample}_oxidase_bats_ref.bam'
    log:
        stdout = 'logs/map_vs_oxidase_bats/{sample}.stdout.log',
        stderr = 'logs/map_vs_oxidase_bats/{sample}.stderr.log'
    params:
        single=False,
        database=True,
        bowtie2Args="--very-sensitive --no-unal"


rule create_IFN_induction_figures:
    input:
        csv="../meta_data/IFN_induction_merge.csv"
    output:
        fn_fig4=config['result_folder'] + config['supp']['ifn_luciferase_fig4'],
        fn_suppFig1=config['result_folder'] + config['supp']['ifn_luciferase_suppFig1'],
        fn_suppFig2=config['result_folder'] + config['supp']['ifn_luciferase_suppFig2'],
    log:
        stdout='logs/create_IFN_induction_figures.stdout.log',
        stderr='logs/create_IFN_induction_figures.stderr.log',
    conda:
        "../../envs/plotting_env.yaml"
    threads:
        2
    resources:
        mem_mb='3000M',
        disk_mb='10000M',
        time="01:00:00"
    script:
        "../scripts/plot_IFN_data.py"