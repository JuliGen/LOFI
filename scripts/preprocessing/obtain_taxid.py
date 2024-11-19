import argparse
import json
import pandas as pd

from Bio import Entrez


def accession_to_taxid(acc: str, db="assembly") -> int:
    """
    Converts an accession number to a taxid from NCBI database.
    Additional function for obtain_taxid function.

    :param acc: accession number
    :param db: database to search
    :return: taxid
    """
    handle = Entrez.esearch(db=db, term=acc)
    record = Entrez.read(handle)
    gb_id = record["IdList"][0]

    handle = Entrez.esummary(db=db, id=gb_id, retmode="json")
    result = json.load(handle)["result"]

    taxid = result[gb_id]["taxid"]

    return taxid


def obtain_taxid(path_dist_df: str) -> (int, pd.DataFrame):
    """
    Obtains taxid from NCBI database based on an input dataframe distances best hit.

    :param df: dataframe with distances calculated with mash
    :return: taxid
    """
    dist_df = pd.read_csv(path_dist_df,
                          sep="\t",
                          names=["Reference-ID", "Query-ID", "Mash-distance", "P-value", "Matching-hashes"])

    dist_df.sort_values(by=["Mash-distance", "P-value"], inplace=True)
    dist_df["Query-ID"] = dist_df["Query-ID"].str.replace("genomes/", "")
    dist_df["Reference-ID"] = dist_df["Reference-ID"].str.replace(".gz", "")

    accession = "_".join(dist_df["Reference-ID"].iloc[0].split("_")[:2])
    taxid = accession_to_taxid(accession)

    return taxid, dist_df


def parse_args():
    parser = argparse.ArgumentParser(
        usage="obtain_taxid.py --genome GENOME.FNA --dist GENOME_DISTANCES.TAB --email EMAIL",
        description="""
        Obtains taxid from NCBI database based on Mash distances best hit.
        Email is required due to the use of the Entrez library.
        """,
    )
    parser.add_argument("--genome", nargs="?", help="genome.fna")
    parser.add_argument("--dist", nargs="?", help="mash_distances.tab")
    parser.add_argument("--email", nargs="?", help="email")

    return parser.parse_args()


if __name__ == "__main__":
    path_genome = parse_args().genome
    path_dist_df = parse_args().dist
    email = parse_args().email

    Entrez.email = email

    taxid, sorted_dist_df = obtain_taxid(path_dist_df)

    path_dist_df_sorted = path_dist_df.replace(".tab", "_sorted.tab")
    sorted_dist_df.to_csv(
        path_dist_df_sorted,
        sep="\t",
    )

    print(f"TAXON ID for provided genome is {taxid}")
