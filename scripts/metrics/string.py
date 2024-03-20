import gzip
import shutil
import wget
import pandas as pd


def get_protein_links(taxon_id: str) -> pd.DataFrame:
    """
    Downloads and unzips a protein score table from the STRING database.
    :param taxon_id: taxon_id for species
    :return: pd.DataFrame with protein combined scores
    """
    suffix = "protein.links.v12.0"
    url = f"https://stringdb-downloads.org/download/protein.links.v12.0/{taxon_id}.{suffix}.txt.gz"
    wget.download(url, bar=None)

    with gzip.open(f"{taxon_id}.{suffix}.txt.gz", "rb") as in_file:
        with open(f"{taxon_id}.{suffix}.txt", "wb") as out_file:
            shutil.copyfileobj(in_file, out_file)

    protein_links = pd.read_csv(f"{taxon_id}.{suffix}.txt", sep=" ")
    protein_links.protein1 = protein_links.protein1.str.replace(taxon_id + ".", "")
    protein_links.protein2 = protein_links.protein2.str.replace(taxon_id + ".", "")

    return protein_links
