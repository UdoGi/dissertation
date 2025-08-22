
#
# rule all:
#     input:
#         # config['result_folder'] + "map_vs_malsoor/malsoor.sam"
#         # config['result_folder'] + "call_ancient_consensus/ancient_consensus.fasta"
#         # config['result_folder'] + "run_spades_trusted_contig/contigs.fasta"
#         maf = config['result_folder'] + "map_vs_malsoor_with_last/malsoor.maf",
#         fasta = config['result_folder'] + "create_ancient_reference/ancient_contig_based_reference.fasta",
#         bam=config['result_folder'] + 'bowtie2_align_vs_contig/ancient_contig_based_reference_vs_6879.bam',
#         consensus_fasta = config['result_folder'] + 'call_contig_align_consensus/ancient_aln_based_reference_vs_6879.fasta',
#         final_ref=config['result_folder'] + config['reference_seq']['udo_ref']
#

# assemble with spades using the anti-map_data
rule run_spades:
    input:
        R1=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_R1.fastp.fastq.gz",
        R2=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_R2.fastp.fastq.gz",
        singles=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_singletons.fastp.fastq.gz"
    output:
        config['result_folder'] + "full_spades_assembly/contigs.fasta"
    conda:
        '../envs/spades_env.yaml'
    log:
        stdout='logs/spades.stdout.log',
        stderr='logs/spades.stderr.log'
    params:
        spades_dir=config['result_folder'] + "full_spades_assembly/"
    threads:
        60
    resources:
        mem_mb='250000M',
        disk_mb='500000M',
        time="8:00:00"
    shell:
        "spades.py -t {threads} -o {params.spades_dir} -1 {input.R1} -2 {input.R2} -s {input.singles} > {log.stdout} 2> {log.stderr}"



# download malsoor references
rule download_malsoor_references:
    output:
        fasta_L=temp('malsoor_l.fasta'),
        fasta_M=temp('malsoor_m.fasta'),
        fasta_S=temp('malsoor_s.fasta'),
        fasta_all=config['result_folder'] + config['reference_seq']['malsoor'],
        fasta_L_aa=temp('malsoor_l_aa.fasta'),
        fasta_M_aa=temp('malsoor_m_aa.fasta'),
        fasta_N_aa=temp('malsoor_N_aa.fasta'),
        fasta_NSs_aa=temp('malsoor_NSs_aa.fasta'),
        fasta_all_aa=config['result_folder'] + config['reference_seq']['malsoor_aa']
    log:
        stdout='logs/download_malsoor_references.stdout.log',
        stderr='logs/download_malsoor_references.stderr.log'
    retries: 10
    shell:
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=KF186497.1\&rettype\=fasta\&retmode\=text -O {output.fasta_L}"
        "> {log.stdout} 2> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=KF186498.1\&rettype\=fasta\&retmode\=text -O {output.fasta_M}"
        ">> {log.stdout} 2>> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=nucleotide\&id\=KF186499.1\&rettype\=fasta\&retmode\=text -O {output.fasta_S}"
        ">> {log.stdout} 2>> {log.stderr};"
        "cat {output.fasta_L} {output.fasta_M} {output.fasta_S} > {output.fasta_all};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=protein\&id\=AHF71068.1\&rettype\=fasta\&retmode\=text -O {output.fasta_L_aa}"
        "> {log.stdout} 2> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=protein\&id\=AHF71069.1\&rettype\=fasta\&retmode\=text -O {output.fasta_M_aa}"
        ">> {log.stdout} 2>> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=protein\&id\=AHF71071.1\&rettype\=fasta\&retmode\=text -O {output.fasta_N_aa}"
        ">> {log.stdout} 2>> {log.stderr};"
        "wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\?"
        "db\=protein\&id\=AHF71070.1\&rettype\=fasta\&retmode\=text -O {output.fasta_NSs_aa}"
        ">> {log.stdout} 2>> {log.stderr};"
        "cat {output.fasta_L_aa} {output.fasta_M_aa} {output.fasta_NSs_aa} {output.fasta_N_aa} "
        "> {output.fasta_all_aa};"



# filter contigs down with diamond. this results in a complete contig for the L and S segment
# and 3 contigs that together build a complete M segment
rule map_contigs_vs_malsoor:
    input:
        contigs=config['result_folder'] + "full_spades_assembly/contigs.fasta",
        fasta_malsoor=config['result_folder'] + config['reference_seq']['malsoor_aa']
    output:
        filtered_contigs = temp(config['result_folder'] + "map_contigs_vs_malsoor/filtered_contigs.fasta"),
        malsoor_db = config['result_folder'] + "map_contigs_vs_malsoor/malsoor.dmnd",
        diamond_out= config['result_folder'] + "map_contigs_vs_malsoor/diamond_out.txt",
        diamond_hit_ids=temp(config['result_folder'] + 'map_contigs_vs_malsoor/potential_hit_fasta_id.txt'),
        likely_contigs=config['result_folder'] + 'map_contigs_vs_malsoor/likely_contigs.fasta'
    conda:
        "../envs/diamond_env.yaml"
    log:
        stdout='logs/map_contigs_vs_malsoor/stdout.log',
        stderr='logs/map_contigs_vs_malsoor/stderr.log'
    shell:
        "seqkit seq -m 600 {input.contigs} > {output.filtered_contigs} 2> {log.stderr};"
        "diamond makedb --in {input.fasta_malsoor} --db {output.malsoor_db} > {log.stdout} 2>> {log.stderr};"
        "diamond blastx --db {output.malsoor_db} -q {output.filtered_contigs} > {output.diamond_out} 2>> {log.stderr};"
        "cut -f1 {output.diamond_out} | sort | uniq > {output.diamond_hit_ids} 2>> {log.stderr};"
        "seqkit grep -f {output.diamond_hit_ids} {output.filtered_contigs} > {output.likely_contigs} 2>> {log.stderr};"


# Reads are mapped with diamond vs malsoor and the matching read ids saved
# these reads together with the incomplete contigs are then used by spades to fully assemble the M-segment
# not all reads are needed here, since spades just needs to stitch together 3 overlapping contigs
rule map_reads_vs_malsoor:
    input:
        singles=config['result_folder'] + "fastp_paired/ALL_D6879_anti_map_singletons.fastp.fastq.gz",
        fasta_malsoor=config['result_folder'] + config['reference_seq']['malsoor_aa'],
        malsoor_db= config['result_folder'] + "map_contigs_vs_malsoor/malsoor.dmnd"
    output:
        diamond_out=config['result_folder'] + "maps_reads_vs_malsoor/diamond_out.txt",
        diamond_hit_ids=temp(config['result_folder'] + 'maps_reads_vs_malsoor/potential_hit_fasta_id.txt'),
        filtered_singles=config['result_folder'] + "maps_reads_vs_malsoor/filtered_singles.fastq",
    conda:
        "../envs/diamond_env.yaml"
    threads:
        64
    log:
        stdout='logs/maps_reads_vs_malsoor/stdout.log',
        stderr='logs/maps_reads_vs_malsoor/stderr.log'
    shell:
        "diamond blastx -p {threads} --db {input.malsoor_db} -q {input.singles} > {output.diamond_out} 2>> {log.stderr}; "
        "cut -f1 {output.diamond_out} | sort | uniq > {output.diamond_hit_ids} 2>> {log.stderr};"
        "seqkit grep -f {output.diamond_hit_ids} {input.singles} > {output.filtered_singles};"


# assemble with spades using the known matching singletons and passing
# the contigs that match malsoor virus, allowing spades to stitch together the M-Segment
rule run_spades_trusted_contig:
    input:
        singles=config['result_folder'] + "maps_reads_vs_malsoor/filtered_singles.fastq",
        trusted_contigs=config['result_folder'] + 'map_contigs_vs_malsoor/likely_contigs.fasta'
    output:
        config['result_folder'] + "run_spades_trusted_contig/contigs.fasta"
    conda:
        '../envs/spades_env.yaml'
    log:
        stdout='logs/run_spades_trusted_contig.stdout.log',
        stderr='logs/run_spades_trusted_contig.stderr.log'
    params:
        spades_dir=config['result_folder'] + "run_spades_trusted_contig/"
    threads:
        60
    resources:
        mem_mb='25000M',
        disk_mb='50000M',
        time="1:00:00"
    shell:
        "spades.py -t {threads} -o {params.spades_dir}  -s {input.singles} --trusted-contigs {input.trusted_contigs}"
        "> {log.stdout} 2> {log.stderr};"


# mapping the contigs against malsoor with "last" to tell me which contig matches which segment and strand
rule map_vs_malsoor_with_last:
    input:
        likely_contigs=config['result_folder'] + "run_spades_trusted_contig/contigs.fasta",
        fasta_malsoor=config['result_folder'] + config['reference_seq']['malsoor']
    output:
        lastdb=config['result_folder'] + "map_vs_malsoor_with_last/malsoor.suf",
        last_train=temp(config['result_folder'] + "map_vs_malsoor_with_last/lastdb.train"),
        maf=config['result_folder'] + "map_vs_malsoor_with_last/malsoor.maf",
    log:
        stdout='logs/map_vs_malsoor_with_last/stdout.log',
        stderr='logs/map_vs_malsoor_with_last/stderr.log'
    params:
        db_name = lambda wildcards, output: output.lastdb.split('.suf')[0]
    conda:
        "../envs/last_env.yaml"
    shell:
        "lastdb {params.db_name} {input.fasta_malsoor} > {log.stdout} 2>> {log.stderr}; "
        "last-train {params.db_name} {input.likely_contigs} > {output.last_train} 2>> {log.stderr}; "
        "lastal -p {output.last_train} {params.db_name} {input.likely_contigs} > {output.maf} 2>> {log.stderr};"


# using the alignment of "last", that tells me about the strand information and which contig matches which segment
# rename the contigs and reverse complement if needed.
# the alignment of "last" is not used to call a consensus, since it soft-clips essential bases!
rule create_ancient_reference:
    input:
        maf=config['result_folder'] + "map_vs_malsoor_with_last/malsoor.maf",
        contigs=config['result_folder'] + "run_spades_trusted_contig/contigs.fasta"
    output:
        fasta=config['result_folder'] + "create_ancient_reference/ancient_contig_based_reference.fasta"
    log:
        stdout='logs/create_ancient_reference/stdout.log',
        stderr='logs/create_ancient_reference/stderr.log'
    conda:
        "../envs/biopython_env.yaml"
    script:
        "../scripts/create_ancient_reference.py"

rule copy_ancient_virus_over:
    input:
        fasta=config['result_folder'] + "create_ancient_reference/ancient_contig_based_reference.fasta"
    output:
        fasta=config['result_folder'] + config['reference_seq']['udo_ref']
    shell:
        "cp {input.fasta} {output.fasta};"


module bowtie_workflow:
    snakefile: "bowtie.smk"
    config: config


# align all ancient reads to manually inspect the alignment in Geneious
use rule run_bowtie2_paired_RNA from bowtie_workflow as bowtie2_align_vs_contig with:
    input:
        R1 = config['result_folder'] + "fastp_paired/ALL_D6879_R1.fastp.fastq.gz",
        R2 = config['result_folder'] + "fastp_paired/ALL_D6879_R2.fastp.fastq.gz",
        singles = config['result_folder'] + "fastp_paired/ALL_D6879_singletons.fastp.fastq.gz",
        reference_fn = config['result_folder'] + "create_ancient_reference/ancient_contig_based_reference.fasta"
    output:
        bam=config['result_folder'] + 'bowtie2_align_vs_contig/ancient_contig_based_reference_vs_6879.bam'
    params:
        single=False,
        # running run_bowtie2.py using a fasta and creating a db
        reference_file=True,
        bowtie2Args="--very-sensitive --no-unal"
    log:
        stdout='logs/bowtie2_align_vs_contig/stdout.log',
        stderr='logs/bowtie2_align_vs_contig/stderr.log'



# call the consensus of this alignment
rule call_contig_align_consensus:
    input:
        bam=config['result_folder'] + 'bowtie2_align_vs_contig/ancient_contig_based_reference_vs_6879.bam'
    output:
        bam_temp=temp(config['result_folder'] + 'call_contig_align_consensus/udo_virus_temp.bam'),
        consensus_fasta=config['result_folder'] + 'call_contig_align_consensus/ancient_aln_based_reference_vs_6879.fasta'
    conda:
        '../envs/alignment_consensus_env.yaml'
    log:
        stdout='logs/call_contig_align_consensus.stdout',
        stderr='logs/call_contig_align_consensus.stderr'
    threads:
        8
    params:
        prefix=config['result_folder'] + 'call_contig_align_consensus/ancient_aln_based_reference_vs_6879',
        minimum_depth=2
    resources:
        mem_mb='32000M',
        disk_mb='10000M',
        time="04:00:00"
    shell:
        "samtools sort {input.bam} -o {output.bam_temp} --threads {threads} >> {log.stdout} 2>> {log.stderr}; "
        "samtools index {output.bam_temp} >> {log.stdout} 2>> {log.stderr}; "
        "samtools view -h {output.bam_temp} L_segment_udo_virus | egrep -v 'M_segment_udo_virus|S_segment_udo_virus' | "
        "samtools mpileup -aa -A -d 0 -Q 0 - | "
        "ivar consensus -p {params.prefix}_L -m {params.minimum_depth} >> {log.stdout} 2>> {log.stderr}; "
        "samtools view -h {output.bam_temp} M_segment_udo_virus | egrep -v 'L_segment_udo_virus|S_segment_udo_virus' | "
        "samtools mpileup -aa -A -d 0 -Q 0 - | "
        "ivar consensus -p {params.prefix}_M -m {params.minimum_depth} >> {log.stdout} 2>> {log.stderr}; "
        "samtools view -h {output.bam_temp} S_segment_udo_virus | egrep -v 'L_segment_udo_virus|M_segment_udo_virus' | "
        "samtools mpileup -aa -A -d 0 -Q 0 - | "
        "ivar consensus -p {params.prefix}_S -m {params.minimum_depth} >> {log.stdout} 2>> {log.stderr}; "
        "cat {params.prefix}_L.fa {params.prefix}_M.fa {params.prefix}_S.fa > {output.consensus_fasta}"




