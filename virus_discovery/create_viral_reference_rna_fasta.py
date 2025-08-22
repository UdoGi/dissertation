import argparse
import os
import sys
import tempfile
from multiprocessing import Queue

import pandas as pd
from Bio import SeqIO
from general_helpers import my_generator, filter_tax_id_mapping, \
    handle_nucleotide_fasta
from general_helpers import download_via_http_link, create_stdout_logger

taxonomy_acc2taxId = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
genbank_and_refseq_nucleotide_sequences = "http://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/AllNucleotide.fa"
genbank_and_refseq_nucleotide_meta_data = "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/AllNuclMetadata.csv.gz"

meta_data_dtypes = {"#Accession": str,
                    "SRA_Accession": str,
                    "Submitters": str,
                    "Release_Date": str,
                    "Species": str,
                    "Genus": str,
                    "Family": str,
                    "Molecule_type": str,
                    "Length": int,
                    "Sequence_Type": str,
                    "Nuc_Completeness": str,
                    "Genotype": str,
                    "Segment": str,
                    "Publications": float,
                    "Geo_Location": str,
                    "USA": str,
                    "Host": str,
                    "Isolation_Source": str,
                    "Collection_Date": str,
                    "BioSample": str,
                    "GenBank_Title": str}


def handle_meta_data(in_fn):
    """
    Filter the meta data file down to RNA viruses and try to exclude phages.
    @param in_fn: csv file from NCBI giving meta info about the viral genome sequences from refseq/genbank
    @return: dataframe containing the csv file filtered down
    """
    df_meta = pd.read_csv(in_fn, dtype=meta_data_dtypes)

    def is_phage(row):
        return (("phage" in str(row.Species).lower())
                or "bacteria" in str(row.Host).lower()
                or "bacterium" in str(row.Host).lower()
                or "phage" in str(row.GenBank_Title).lower()
                or "bacterium" in str(row.GenBank_Title).lower())

    df_meta['is_phage'] = df_meta.apply(is_phage, axis=1)
    df_meta['is_rna_virus'] = df_meta.Molecule_type.apply(lambda x: "rna" in str(x).lower())

    df_meta = df_meta[(df_meta.Species != 'Severe acute respiratory syndrome-related coronavirus')
                      | (df_meta['#Accession'] == "NC_045512.2")]

    df_complete = df_meta[(df_meta.Nuc_Completeness == 'complete')
                          & (~df_meta.is_phage)
                          & df_meta.is_rna_virus]
    return df_complete


def parse_meta_data(args, my_logger):
    """
    Either downloads or gets metadata about the corresponding fasta file.
    Then it filters down to my chosen criteria.
    :param args: args from the argparser
    :param my_logger: reference to a logger
    :return: filtered data frame of the meta data
    """
    if args.meta:
        my_logger.info('meta data provided, reading and filtering it next ...')
        meta_df = handle_meta_data(args.meta)
    else:
        my_logger.info('meta data not provided, downloading it next ...')
        fp = tempfile.NamedTemporaryFile(delete=False)
        download_via_http_link(genbank_and_refseq_nucleotide_meta_data, out_fn=fp.name)
        my_logger.info('downloading done, next filtering....')
        meta_df = handle_meta_data(fp.name)
        fp.close()
        os.remove(fp.name)
    my_logger.info('finished handling of meta data')
    if args.out_fn_meta:
        my_logger.info(f"writing out meta data to: {args.out_fn_meta}")
        meta_df.to_csv(args.out_fn_meta, index=False)
    return meta_df


def create_arg_parser(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('out_fn', help='output file name for the FASTA file')
    parser.add_argument("--meta", help='file with ncbi meta data')
    parser.add_argument("--nuc", help='file with ncbi sequence records as FASTA')
    parser.add_argument("--acc2tax", help='file with accession2taxid mapping')
    parser.add_argument("--out_fn_meta", help="output the filtered meta data to this file")
    args = parser.parse_args(argv)
    return args


def main(argv):
    args = create_arg_parser(argv)
    my_logger = create_stdout_logger()

    # meta data
    meta_df = parse_meta_data(args, my_logger)
    sub_set = set(meta_df["#Accession"])

    # tax id and nucs
    my_queue = Queue()
    temp_file_list = list()  # needed as only after the spawned processes are done i can delete these files

    p1 = filter_tax_id_mapping(args, my_logger, my_queue, sub_set, temp_file_list)
    p2 = handle_nucleotide_fasta(args, my_logger, my_queue, sub_set, temp_file_list)

    qe1, qe2 = my_queue.get(), my_queue.get()
    p1.join(), p2.join()
    [(fp.close(), os.remove(fp.name)) for fp in temp_file_list]

    my_logger.info('done with parsing both, writing out the result next')
    # figure out who came out the queue first
    if qe1[0] == 'seq_dic':
        seq_dict, tax_dic = qe1[1], qe2[1]
    elif qe1[0] == 'tax_id':
        seq_dict, tax_dic = qe2[1], qe1[1]

    with open(args.out_fn, 'w') as handle:
        SeqIO.write(my_generator(seq_dict, tax_dic, meta_df), handle, "fasta")
    my_logger.info('done with everything')


if __name__ == "__main__":
    main(sys.argv[1:])
