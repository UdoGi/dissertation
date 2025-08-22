
from subprocess import run, CalledProcessError
import logging
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gzip
import sys
from multiprocessing import Process
from typing import Set, Dict


taxonomy_acc2taxId = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
genbank_and_refseq_nucleotide_sequences = "http://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/AllNucleotide.fa"
genbank_and_refseq_nucleotide_meta_data = "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/AllNuclMetadata.csv"


def download_via_http_link(http_link, out_fn=None):
    """
    Given a http link, this function downloads content from the link
    calling wget in a separate process and optionally downloads to a specified location.
    If no location is specified, it will just download to the current directory.
    @param http_link: link to download from
    @param out_fn: output file name
    @return:
    """
    params = ['wget', f"{http_link}", '-O', f"{out_fn}"] if out_fn else ['wget', f"{http_link}"]
    try:
        run(params,
            stdout=sys.stdout,
            stderr=sys.stderr,
            check=True)
    except CalledProcessError as ex:
        print(f"wget returned with exit code \
                            {ex.returncode}", file=sys.stderr)
        exit(ex.returncode)


def create_stdout_logger(level=logging.DEBUG):
    """
    Creates a logger that writes to stdout.
    Just nicer than just printing to std out, since it includes time stamps
    @return:
    """
    my_logger = logging.getLogger()
    my_logger.setLevel(logging.DEBUG)
    my_logger.addHandler(create_stdout_handler(level))
    return my_logger

def create_stdout_handler(level=logging.DEBUG):
    """
    Creates a logger handler that writes to stdout.
    Just nicer than just printing to std out, since it includes time stamps
    @return: rtype:
    """
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    return handler


def my_generator(seq_dict, tax_lookup, meta_df=None):
    """
    Generator taking a dictionary of fasta sequences and a dictionary with accession number
    as key and tax id as value.
    The generator iterates over the dictionary content, if the accession number of the sequence
    is in the tax_lookup dictionary it will add a kraken-specific signature to the fasta description
    that includes the tax-id and yield that modified fasta entry.
    @param seq_dict: dictionary containing fasta sequences {accession: SeqRecord}
    @param tax_lookup: dictionary containing tax-ids entries {accession: tax-id}
    @param meta_df: dataframe with meta information of the sequence
    @return: yields fasta-records with a kraken-specific description
    """
    for k, record in seq_dict.items():
        if k in tax_lookup:
            tax_id = tax_lookup[k]
            new_name = f"{k}|kraken:taxid|{tax_id}"
            if meta_df is not None:
                df = meta_df[meta_df["#Accession"] == k]
                assert len(df) == 1, f"meta_df contains more than one entry or none for accession: {record.id}"
                new_description = f"{new_name} {df.GenBank_Title.values[0]}"
            else:
                # extracting the description only without the id
                t = record.description.split(" ", maxsplit=1)
                new_description = new_name if len(t) == 1 else f"{new_name} {t[-1]}"
            yield SeqRecord(record.seq, id=new_name, description=new_description)


def filter_tax_id_mapping(args, my_logger, my_queue, sub_set, temp_file_list):
    # tax id mapping
    if args.acc2tax:
        p1 = Process(target=read_in_acc_to_taxid, args=(args.acc2tax,), kwargs={"filter_set": sub_set,
                                                                                "queue": my_queue})
        p1.start()
        my_logger.info("starting the acc2tax parsing ...")
    else:
        fp_acc2tax = tempfile.NamedTemporaryFile(delete=False)
        temp_file_list.append(fp_acc2tax)
        my_logger.info("downloading the acc2tax id mappings ...")
        p1 = Process(target=download_via_http_link, args=(taxonomy_acc2taxId,), kwargs={'out_fn': fp_acc2tax.name})
        p1.start(), p1.join()

        p1 = Process(target=read_in_acc_to_taxid, args=(fp_acc2tax.name,), kwargs={"filter_set": sub_set,
                                                                                   "queue": my_queue})
        p1.start()
        my_logger.info("starting the acc2tax parsing ...")
    return p1

def handle_nucleotide_fasta(args, my_logger, my_queue, sub_set, temp_file_list):
    # nucs
    if args.nuc:
        p2 = Process(target=filter_nuc_sequences, args=(args.nuc, sub_set), kwargs={"queue": my_queue})
        p2.start()
        my_logger.info("starting to parse the nuc sequences ...")
    else:
        fp_nucs = tempfile.NamedTemporaryFile(delete=False)
        temp_file_list.append(fp_nucs)
        p2 = Process(target=download_via_http_link, args=(genbank_and_refseq_nucleotide_sequences,),
                     kwargs={'out_fn': fp_nucs.name})
        p2.start(), p2.join()

        p2 = Process(target=filter_nuc_sequences, args=(fp_nucs.name, sub_set), kwargs={"queue": my_queue})
        p2.start()
    return p2


def filter_nuc_sequences(nuc_fn: str, filter_set: Set[str], queue=None) -> Dict[str, SeqRecord]:
    """
    Filters out FASTA records from a FASTA file whose Ids are not in the given set and
    returns the records in a dict[str, SeqRecord] with the fasta id as key.
    :param nuc_fn: file name of the fasta file
    :param filter_set: set with allowed ids (whitelist)
    :param queue: necessary if called within a separate process to return the dict
    :return: dict[str, SeqRecord] with the fasta id as key
    """
    ncbi_dict = dict()
    for record in SeqIO.parse(open(nuc_fn), 'fasta'):
        accession = record.id.split('|')[0]
        if accession in filter_set:
            ncbi_dict[accession] = record
    queue.put(('seq_dic', ncbi_dict)) if queue else None
    return ncbi_dict


def read_in_acc_to_taxid(in_fn, filter_set=None, queue=None):
    """
    Iterates over the content of the accession2taxid file that links
    accession ids to tax ids and returns a dictionary with accession:[taxId]
    If filter_set is given only accessions that are in the set are returned.
    @param in_fn:
    @param filter_set:
    @param queue: needs to be specified if run in a separate process to be able to return sth
    @return:
    """
    rd = dict()
    f = gzip.open(in_fn, 'rt') if is_gz_file(in_fn) else open(in_fn)
    with f as fp:
        # accession accession.version   taxid   gi
        fp.readline()  # skip header
        for line in fp:
            _, accession, tax_id, _ = line.strip().split('\t')
            if filter_set:
                if accession in filter_set:
                    assert accession not in rd, f"accession {accession} occurred twice in input file"
                    rd[accession] = tax_id
            else:
                assert accession not in rd, f"accession {accession} occurred twice in input file"
                rd[accession] = tax_id
    queue.put(('tax_id', rd)) if queue is not None else None
    return rd


def is_gz_file(filepath):
    """
    Hacky test for the magic number in the 2 first bytes of the file to determine if its gzipped.
    @param filepath:
    @return:
    """
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'