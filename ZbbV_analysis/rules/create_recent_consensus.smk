

# df = pd.read_excel(config['fastq']['fastq_merged_meta_fn'], engine='openpyxl')
# cases = {
#                 'ALL_D9595', 'ALL_D9598', 'ALL_D10643',
#                'ALL_D10376', 'ALL_D10116',
#                 'ALL_D10118',
#                 'ALL_D6879',
#                 'ALL_D9591', 'ALL_D10758', 'ALL_D9593'
#                }

# print( df.sample_name_snake.unique())
# rule all:
#     input:
#         [config['result_folder'] + f'udo_ref_align_consensus/{sample}_udo_ref_align_consensus.fasta'
#             for sample in df.sample_name_snake.unique()]
#

module bowtie_workflow:
    snakefile: "bowtie.smk"
    config: config


use rule make_bowtie2_db_RNA from bowtie_workflow as make_udo_virus_bowtie_db with:
    input:
        fasta=config['result_folder'] + config['reference_seq']['udo_ref']
    output:
        f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_udo_ref']}/{config['bowtie2']['bowtie2_db_name_udo']}.1.bt2"
    log:
        stdout = config['result_folder'] + 'logs/make_udo_virus_bowtie_db/stdout.log',
        stderr= config['result_folder'] + 'logs/make_udo_virus_bowtie_db/stderr.log',
    params:
        db_dir=config['result_folder'] + config['bowtie2']["bowtie2_db_folder_udo_ref"],
        db_name=config['bowtie2']["bowtie2_db_name_udo"]


use rule run_bowtie2_paired_RNA from bowtie_workflow as map_vs_udo_virus with:
    input:
        R1=config['result_folder'] + "fastp_paired/{sample}_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/{sample}_R2.fastp.fastq.gz",
        singles=config['result_folder'] + 'fastp_paired/{sample}_singletons.fastp.fastq.gz',
        index=f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_udo_ref']}/{config['bowtie2']['bowtie2_db_name_udo']}.1.bt2"
    output:
        bam = config['result_folder'] + 'map_vs_udo_virus/{sample}_udo_ref.bam'
    log:
        stdout = 'logs/map_vs_udo_virus/{sample}.stdout.log',
        stderr = 'logs/map_vs_udo_virus/{sample}.stderr.log'
    params:
        single=False,
        database=True,
        bowtie2Args="--very-sensitive --no-unal"


rule udo_ref_align_consensus:
    input:
        bam=config['result_folder'] + 'map_vs_udo_virus/{sample}_udo_ref.bam'
    output:
        touch('flags/{sample}_consensus.done'),
        bam_indexed=config['result_folder'] + 'udo_ref_align_consensus/{sample}_udo_ref_sorted.bam',
        fasta_temp_L=temp(config['result_folder'] + 'udo_ref_align_consensus/{sample}_L.fa'),
        fasta_temp_M=temp(config['result_folder'] + 'udo_ref_align_consensus/{sample}_M.fa'),
        fasta_temp_S=temp(config['result_folder'] + 'udo_ref_align_consensus/{sample}_S.fa'),
        consensus_fasta=config['result_folder'] + 'udo_ref_align_consensus/{sample}_udo_ref_align_consensus.fasta'
    conda:
        '../../envs/alignment_consensus_env.yaml'
    log:
        stdout='logs/{sample}_udo_ref_align_consensus.stdout',
        stderr='logs/{sample}_udo_ref_align_consensus.stderr'
    threads:
        8
    params:
        prefix=config['result_folder'] + '/udo_ref_align_consensus',
        minimum_depth=2,
        sample_name=lambda wildcards: wildcards.sample
    resources:
        mem_mb='32000M',
        disk_mb='10000M',
        time="04:00:00"
    shell:
        "samtools sort {input.bam} -o {output.bam_indexed} --threads {threads} >> {log.stdout} 2>> {log.stderr}; "
        "samtools index {output.bam_indexed} >> {log.stdout} 2>> {log.stderr}; "
        "samtools view -h {output.bam_indexed} L_segment_udo_virus | egrep -v 'SN:M_segment_udo_virus|SN:S_segment_udo_virus' | "
        "samtools mpileup -aa -A -d 0 -Q 0 - | "
        "ivar consensus -p {params.prefix}/{params.sample_name}_L -m {params.minimum_depth} >> {log.stdout} 2>> {log.stderr}; "
        "samtools view -h {output.bam_indexed} M_segment_udo_virus | egrep -v 'SN:S_segment_udo_virus|SN:L_segment_udo_virus' | "
        "samtools mpileup -aa -A -d 0 -Q 0 - | "
        "ivar consensus -p {params.prefix}/{params.sample_name}_M -m {params.minimum_depth} >> {log.stdout} 2>> {log.stderr}; "
        "samtools view -h {output.bam_indexed} S_segment_udo_virus | egrep -v 'SN:M_segment_udo_virus|SN:L_segment_udo_virus' | "
        "samtools mpileup -aa -A -d 0 -Q 0 - | "
        "ivar consensus -p {params.prefix}/{params.sample_name}_S -m {params.minimum_depth} >> {log.stdout} 2>> {log.stderr}; "
        "cat {output.fasta_temp_L} {output.fasta_temp_M} {output.fasta_temp_S} > {output.consensus_fasta}"
