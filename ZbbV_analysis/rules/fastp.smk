

rule fastp_paired:
    input:
        R1="",
        R2=""
    output:
        R1=config['result_folder'] + "fastp_paired/{sample}_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/{sample}_R2.fastp.fastq.gz",
        unpaired=temp(config['result_folder'] + 'fastp_paired/{sample}_unpaired.fastp.fastq.gz'),
        merge=temp(config['result_folder'] + 'fastp_paired/{sample}_merged.fastp.fastq.gz'),
        singletons=config['result_folder'] + 'fastp_paired/{sample}_singletons.fastp.fastq.gz',
        json='fastp_paired' + '/{sample}.fastp.json',
        html='fastp_paired' + '/{sample}.html'
    conda:
        "../../envs/preprocess.yaml"
    log:
        stdout="logs/fastp_paired/{sample}.stdout.log",
        stderr="logs/fastp_paired/{sample}.stderr.log"
    threads:
        16
    resources:
        mem_mb='32000M',
        disk_mb='10000M',
        time="01:00:00"
    params:
        basename=lambda wildcards: wildcards.sample,
        min_length="30"
    shell:
         "fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} "
         "--unpaired1 {output.unpaired} --unpaired2 {output.unpaired} --merge "
         "--merged_out {output.merge} --dedup --length_required {params.min_length} -p "
         "--json {output.json} --html {output.html} --thread {threads} > {log.stdout} 2> {log.stderr};"
         "touch {output.unpaired} {output.merge};"
         "cat {output.unpaired} {output.merge} > {output.singletons}"
