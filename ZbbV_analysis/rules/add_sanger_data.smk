
import sys
sys.path.append('../scripts')
from helpers import create_Nsplit_fasta
from pathlib import Path

# sangered_samples = ['ALL_D10758']
#
# rule all:
#     input:
#         [f'flags/{sample}_consensus_sanger.done' for sample in sangered_samples]


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


rule split_aln_based_consensus:
    input:
        consensus_fasta=config['result_folder'] + 'udo_ref_align_consensus/{sample}_udo_ref_align_consensus.fasta'
    output:
        split_consensus_fasta=config['result_folder'] + 'udo_ref_align_consensus_split/{sample}_udo_ref_align_consensus_SPLIT.fasta'
    threads:
        1
    run:
        create_Nsplit_fasta(fn_in=input.consensus_fasta, fn_out=output.split_consensus_fasta)

def get_sanger_fastas(wildcards):
    return [str(p) for p in Path(f'../sanger_data/{wildcards.sample}/').glob("*.fasta")]

rule concat_all_fastas:
    input:
        split_consensus_fasta = config['result_folder'] + 'udo_ref_align_consensus_split/{sample}_udo_ref_align_consensus_SPLIT.fasta',
        all_sanger_fastas=get_sanger_fastas
    output:
        all_fasta='../sanger_data/merged_fasta/{sample}.fasta'
    shell:
        "cat {input.split_consensus_fasta} {input.all_sanger_fastas} > {output.all_fasta}"


rule map_sanger_against_ref:
    input:
        all_fasta = '../sanger_data/merged_fasta/{sample}.fasta',
        index=f"{config['result_folder'] + config['bowtie2']['bowtie2_db_folder_udo_ref']}/{config['bowtie2']['bowtie2_db_name_udo']}.1.bt2"
    output:
        bam=config['result_folder'] + "sanger_map_vs_udo/{sample}_sanger_map.bam"
    log:
        stdout="logs/map_sagner_against_ref/{sample}.stdout",
        stderr="logs/map_sagner_against_ref/{sample}.stderr",
    conda:
        "../../envs/alignment_consensus_env.yaml"
    threads:
        2
    params:
        db=lambda wildcards, input: input.index.split('.1.bt2')[0]
    shell:
        "bowtie2 -x {params.db} -f {input.all_fasta} -p {threads} --very-sensitive "
        "-D 35 -R 12 -N 1 -L 10 --score-min L,-0.6,-1 --no-unal | "
        "samtools view -bS - > {output.bam} 2>> {log.stderr};"


rule udo_ref_align_consensus_sanger:
    input:
        bam = config['result_folder'] + "sanger_map_vs_udo/{sample}_sanger_map.bam"
    output:
        touch('flags/{sample}_consensus_sanger.done'),
        bam_indexed=config['result_folder'] + 'udo_ref_align_consensus_sanger/{sample}_udo_ref_sorted.bam',
        fasta_temp_L=temp(config['result_folder'] + 'udo_ref_align_consensus_sanger/{sample}_L.fa'),
        fasta_temp_M=temp(config['result_folder'] + 'udo_ref_align_consensus_sanger/{sample}_M.fa'),
        fasta_temp_S=temp(config['result_folder'] + 'udo_ref_align_consensus_sanger/{sample}_S.fa'),
        consensus_fasta=config['result_folder'] + 'udo_ref_align_consensus_sanger/{sample}_udo_ref_align_consensus_sanger.fasta'
    log:
        stdout='logs/{sample}_udo_ref_align_consensus_sanger.stdout',
        stderr='logs/{sample}_udo_ref_align_consensus_sanger.stderr'
    conda:
        '../../envs/alignment_consensus_env.yaml'
    params:
        prefix=config['result_folder'] + '/udo_ref_align_consensus_sanger',
        minimum_depth=1,
        sample_name=lambda wildcards: wildcards.sample.split('_ALL_')[0]
    threads:
        8
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
