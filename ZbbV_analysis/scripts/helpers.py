from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


def yield_renamed_records(fn_p, df):
    """
    Given a fasta file that contains the CDS consensus sequences,
    rename the fasta ids and keep the old ID in the fasta description
    :param fn_p: fasta file with CDS of L,M,N,NSs of udo-virus
    :param df: dataframe with metadata on bat samples
    :return: yields renamed fasta records
    """
    df = df.set_index('caseId')
    with fn_p.open() as fp: 
        cds = fn_p.stem.split('_')[0]
        for record in SeqIO.parse(fp, 'fasta'):
            rid = record.id
            assert re.match(r'Consensus_ALL_D\d+_[SML]_threshold', rid), \
                f'fasta id {record.id} not what i expected'
            sp = rid.split('_')
            case_id = sp[2]

            assert len(df.loc[[case_id]]) == 1, f'more than one sample-id for this case id {case_id}'

            row = df.loc[case_id]
            renamed_id = f"{row['location'].replace(' ', '_')}_{row['time of death']}_{cds}"
            yield SeqRecord(record.seq, id=renamed_id, description=rid)


def split_fasta_entries_generator(fn):
    """
    Given a fasta file, split each sequence by N's and return each fragment free of N's
    as a separate fasta record. The fasta-ID is extended by a suffix specifying the index of the fragment
    :param fn: fasta file to be splitted
    :return:
    """
    with open(fn) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            fragments = [e for e in record.seq.split('N') if e != '']
            for i, fragment in enumerate(fragments):
                yield SeqRecord(fragment, id=f'{record.id}_fragment_{i + 1}')


def create_Nsplit_fasta(fn_in, fn_out):
    """
    Given a fasta file, split the sequence if N's are present, such that new fragments of the passed
    sequence are created without N's
    :param fn_in: fasta file to be split
    :param fn_out: fasta file that contains multiple entries for each determined fragment
    :return:
    """
    with open(fn_out, 'wt') as fh:
        SeqIO.write(split_fasta_entries_generator(fn_in), fh, 'fasta')


def merge_refs_for_overview_tree(ext_refs, consensus_dir, sanger_dir, used_cases, fn_out):
    """
    :param ext_refs: List of bandavirus reference files in fasta format downloaded from ncbi
    :param consensus_dir: Path to directory containing consensus fasta files
    :param sanger_dir: Path to directory containing Sanger fasta files
    :param used_cases: List of cases (samples) to be used, whitelist
    :param fn_out: Output file name to store merged records of bandavirus plus my genomes
    :return: None

    Taking the bandavirus fasta files and directories of sanger and consensus sequences,
    check for each case-id if it is in the whitelist and if a sanger-sequence extended consensus exist.
    then merge all the bandavirus fasta files with either the regular consensus or the sanger-extended consensus
    into a single fasta file.
    """
    cons_d = {p.stem.split('_')[1] : p for p in Path(consensus_dir).glob("*.fasta")}
    sanger_d = {p.stem.split('_')[1] : p for p in Path(sanger_dir).glob("*.fasta")}

    # if the case is sangered the consensus from there is used
    final_list = [sanger_d.get(case) or cons_d.get(case) for case in used_cases]

    record_list = []
    for fn in final_list:
        with fn.open() as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                if "_L_" in record.id:
                    record_list.append(record)

    for fn in ext_refs:
        with open(fn) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                record_list.append(record)

    with open(fn_out, 'wt') as fo:
        SeqIO.write(record_list, fo, 'fasta')


def validate_consensi(consensus_dir, sanger_dir, ref_cds_dir, used_cases):
    """
    Given the computed consenus sequences that are either in the `consensus dir`
    or the `sanger dir`, validate that they contain the manually validated CDSs.
    :param consensus_dir: directory containing consensus fasta files
    :param sanger_dir: directory containing consensus fasta files that are extended by sangered fragments.
    they are preferred over the regular consensus, as they are more complete.
    :param ref_cds_dir: directory containing manually validated CDSs.
    :param used_cases: set of case-IDs that are to be checked
    :return:
    """
    lookup = {'L': {}, 'M': {}, 'N': {}, 'NSs': {}}
    for cds in Path(ref_cds_dir).glob("*.fasta"):
        if "_CDS_consensus" in cds.stem:
            seg = cds.stem.split('_')[0]
            with cds.open() as fh:
                for record in SeqIO.parse(fh, 'fasta'):
                    m = re.match(r"Consensus_(ALL_D\d*)_[LMS]_threshold_\d*_quality_\d*_extraction", record.id)
                    assert m, f'unknown record id: {record.id}'
                    case = m.group(1)
                    lookup[seg][case] = record.seq

    cons_d = {'_'.join(p.stem.split('_')[:2]): p for p in Path(consensus_dir).glob("*.fasta")}
    sanger_d = {'_'.join(p.stem.split('_')[:2]): p for p in Path(sanger_dir).glob("*.fasta")}
    print(cons_d)
    print(sanger_d)
    # if the case is sangered the consensus from there is used
    final_list = [sanger_d.get(case) or cons_d.get(case) for case in used_cases]
    print(final_list)
    # iterate through all created consensi and test if the manually determined CDS is present
    for fn in final_list:
        with fn.open() as fh:
            segment_counter = 0
            for record in SeqIO.parse(fh, 'fasta'):
                m = re.match(r"Consensus_(ALL_D\d*)_([LMS])_threshold_\d*_quality_\d*", record.id)
                assert m, f'unknown record id: {record.id}'
                case, seg = m.group(1), m.group(2)
                if seg == "L":
                    segment_counter+=1
                    assert lookup['L'][case] in record.seq
                elif seg == 'M':
                    segment_counter+=1
                    assert lookup['M'][case] in record.seq
                elif seg == 'S':
                    segment_counter+=1
                    assert lookup['N'][case] in record.seq
                    assert lookup['NSs'][case] in record.seq.reverse_complement(), f"NSs of case {case} not found in {record.id}"
            assert segment_counter == 3, f'Not all three segments are present in {fn.stem}'
    return True