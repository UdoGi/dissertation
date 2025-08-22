import sys
from Bio import AlignIO
from Bio import SeqIO


malsoor_segment_mapping = {
                            "KF186497.1": "L_segment_udo_virus",
                            "KF186498.1": "M_segment_udo_virus",
                            "KF186499.1": "S_segment_udo_virus"
                            }

def get_strand_and_name_info_from_maf(fn_maf):
    """
    Given a MAF file with alignments of Malsoor virus and Penzlin1919 virus contigs,
    determine the strand of the contig and the matching Malsoor virus segment.
    :param fn_maf:
    :return: two dictionaries: `strand_info` mapping contig-ID to strand and `name_info`
    mapping contig-ID to the matching Malsoor virus segment.
    """
    strand_info = {}
    name_info = {}
    with open(fn_maf) as fh:
        for multiple_alignment in AlignIO.parse(fh, "maf"):
            sequences = list(multiple_alignment)
            assert len(sequences) == 2, 'more than two sequences aligned. should just be contig and reference'
            contig_seq, ref_seq = identify_contig_and_reference_sequence(sequences)
            assert contig_seq.name not in strand_info, f'Contig {contig_seq.name} appears in more than one alignment'
            print(contig_seq.name, contig_seq.annotations["strand"])
            if ref_seq.name == "KF186499.1":    # This NCBI entry is in the wrong orientation
                strand_info[contig_seq.name] = contig_seq.annotations["strand"] * -1
            else:
                strand_info[contig_seq.name] = contig_seq.annotations["strand"]
            name_info[contig_seq.name] = ref_seq.name
    return strand_info, name_info


def identify_contig_and_reference_sequence(sequences):
    """
    Function to determine which sequence is the contig sequence and which one is the reference sequence given
    a list of sequences. This function assumes the contig sequence has "NODE_" in its fasta id as spades assemblies do.
    :param sequences:
    :return:
    """
    if 'NODE_' in sequences[1].name:
        return sequences[1], sequences[0]
    elif 'NODE_' in sequences[0].name:
        return sequences[0], sequences[1]
    else:
        assert False, 'NODE_ has to be in the alignment, this is only meant for contig vs malsoor analysis!'


def yield_records(fasta_file, strand_lookup, name_lookup):
    """
    Taking in a fasta file with the contigs and two dictionaries specifying
    the strand and names of the matching reference, this generator yields
    fasta records that have their sequence in the correct orientation (strand)
    and my given names in their fasta ids.
    :param fasta_file: containing the contigs from spades
    :param strand_lookup: mapping of contig name to strand
    :param name_lookup:  mapping of contig name to matching reference name
    :return:
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_id = record.id
        if strand_lookup[fasta_id] == -1:
            record.seq = record.seq.reverse_complement()

        record.id = malsoor_segment_mapping[name_lookup.get(fasta_id)]
        yield record


def rename_and_reverse_complement_contigs(input_filename, output_filename, strand_lookup, name_lookup):
    with open(input_filename, 'r') as fasta_file, open(output_filename, 'w') as outfile:
        SeqIO.write(yield_records(fasta_file, strand_lookup, name_lookup),
                    outfile, 'fasta')


if __name__ == '__main__':
    with open(snakemake.log.stdout, "w") as stdout_log_file, open(snakemake.log.stderr, "w") as stderr_log_file:
        sys.stdout, sys.stderr = stdout_log_file, stderr_log_file
        fn_maf = snakemake.input.maf
        fn_contigs = snakemake.input.contigs
        fasta_out = snakemake.output.fasta
        strand_lookup, name_lookup = get_strand_and_name_info_from_maf(fn_maf)
        rename_and_reverse_complement_contigs(fn_contigs, fasta_out, strand_lookup, name_lookup)

