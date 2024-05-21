import argparse
import gzip
import os
import shutil

import wget


def download_string_files(taxon_id: int) -> None:
    """
    Downloads and unzips a protein score table and protein sequences
    for provided taxon_id from the STRING database.

    :param taxon_id: taxon_id for analyzing species
    """
    db_version = "v12.0"
    keywords = [["sequences", "fa"], ["links", "txt"]]
    output_dir = f"results/{taxon_id}/string"
    urls = []

    for keyword in keywords:
        urls.append(
            f"https://stringdb-downloads.org/download/protein.{keyword[0]}.{db_version}/"
            f"{taxon_id}.protein.{keyword[0]}.{db_version}.{keyword[1]}.gz"
        )

    os.makedirs(
        output_dir, exist_ok=True
    )  # `if not os.exists(output_dir)` doesn't work with snakemake

    for url in urls:
        wget.download(url, out=output_dir, bar=None)

    for keyword in keywords:
        with gzip.open(
            f"{output_dir}/{taxon_id}.protein.{keyword[0]}.{db_version}.{keyword[1]}.gz",
            "rb",
        ) as in_file:
            with open(
                f"{output_dir}/{taxon_id}.protein.{keyword[0]}.{db_version}.{keyword[1]}",
                "wb",
            ) as out_file:
                shutil.copyfileobj(in_file, out_file)


def parse_args():
    parser = argparse.ArgumentParser(
        usage="download_string_files.py --taxid TAXON_ID",
        description="""Downloads files for subsequent analysis from the STRING database.""",
    )
    parser.add_argument("--taxid", nargs="?", help="taxon_id")

    return parser.parse_args()


if __name__ == "__main__":
    taxon_id = parse_args().taxid
    download_string_files(taxon_id)
    print("Files from STRING database downloaded and unzipped successfully")
