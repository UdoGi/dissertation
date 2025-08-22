
rule make_bowtie2_db_RNA:
    input:
        fasta=""
    output:
        ""
    log:
        stdout="",
        stderr="",
    conda:
        "../../envs/alignment_consensus_env.yaml"
    resources:
        mem_mb='32000M',
        disk_mb='10000M',
        time="04:00:00"
    params:
        db_dir = "",
        db_name = ""
    threads:
        32
    shell:
        "rm -rf {params.db_dir} > {log.stdout} 2> {log.stderr} ; mkdir {params.db_dir} >> {log.stdout} 2>> {log.stderr};"
        "cd {params.db_dir}; "
        "bowtie2-build {input.fasta} {params.db_name} --threads {threads} >> {log.stdout} 2>> {log.stderr}; "


rule run_bowtie2_paired_RNA:
    input:
        R1="",
        R2="",
        singles="",
        index="",
        reference_fn=''
    output:
        bam=config['result_folder'] + '/bowtie2_paired/{sample}_RNA.bam'
    log:
        stdout='{sample}.stdout.log',
        stderr='{sample}.stderr.log'
    conda:
        "../../envs/alignment_consensus_env.yaml"
    threads:
        16
    resources:
        mem_mb='32000M',
        disk_mb='10000M',
        time="04:00:00"
    params:
        single=False,
        database=True,
        bowtie2Args="--very-sensitive-local --no-unal"
    script:
        "../../scripts/run_bowtie2.py"




