import pandas as pd
import psycopg2
import os
from ConfigParser import SafeConfigParser


class DBConnection(object):
    """
    Opens connection with database based on the development.ini credentials
    """
    def _get_db_url(self):
        p = SafeConfigParser()
        p.read(os.path.join(os.getenv("HOME"), "development.ini"))

        return p.get('Database', 'sqlalchemy.url')

    def __enter__(self):
        self.conn = psycopg2.connect(self._get_db_url())
        return self.conn

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()


def query_by_chain(chain):

    query = """select name as gene_name,
                       allele as allele,
                       sigpend as signal_peptide_end,
                       cdna as cdna_seq,
                       cds as cds_seq,
                       prot as pep_seq
                       from binf_germline_seq                   
                where release_id = '383364251'
                and name like '{}%';""".format(chain)

    return query


def create_fasta_files(chain, query_results, path, is_gapped):

    if len(query_results) > 0:
        with open(path + "/" + chain + "_AA_new.fasta", "w") as file_aa:
            with open(path + "/" + chain + "_NT_new.fasta", "w") as file_nt:
                for index in query_results.index:
                    result = query_results.iloc[index]

                    parsed_results = parse_germline_sequences(result["gene_name"],
                                                              result["allele"],
                                                              result["cds_seq"],
                                                              result["pep_seq"],
                                                              result["signal_peptide_end"],
                                                              is_gapped)

                    file_aa.write(">" + str(index) + "|" + parsed_results[0] + "|" + " |" * 13 + "\n" + parsed_results[1] + "\n")
                    file_nt.write(">" + str(index) + "|" + parsed_results[0] + "|" + " |" * 13 + "\n" + parsed_results[2] + "\n")

                file_nt.close()
            file_aa.close()


def parse_germline_sequences(gene_name, allele, nt_seq, aa_seq, sig_pep_end, gapped=False):

    # allele_name = get_proper_allele_format(allele)
    full_name = gene_name + "*" + allele

    if sig_pep_end == -1 or pd.isnull(sig_pep_end):
        sequence_aa = aa_seq
        sequence_nt = nt_seq
    elif sig_pep_end != -1 and not pd.isnull(sig_pep_end):
        sequence_aa = aa_seq[int(sig_pep_end):]
        sequence_nt = nt_seq[int(sig_pep_end)*3:]

    if gapped is False:
        sequence_aa = sequence_aa.replace(".", "")
        sequence_nt = sequence_nt.replace(".", "")

    return (full_name, sequence_aa, sequence_nt)


def get_germline_seqs_from_db(path):
    with DBConnection() as conn:
        query_result = pd.read_sql_query(query_by_chain("IGH"), conn)
        create_fasta_files("IGH", query_result, path, True)

        query_result = pd.read_sql_query(query_by_chain("IGL"), conn)
        create_fasta_files("IGL", query_result, path, True)

        query_result = pd.read_sql_query(query_by_chain("IGK"), conn)
        create_fasta_files("IGK", query_result, path, True)


if __name__ == "__main__":
    get_germline_seqs_from_db("../data")