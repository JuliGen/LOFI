import argparse

import pandas as pd
from BCBio import GFF


def parse_gff(path: str) -> pd.DataFrame:
    """
    Parses file.gff from given path.
    :param path: path to file
    :return: pd.DataFrame with the columns necessary for subsequent analysis
    """
    info = []
    limit_info = dict(gff_type=["CDS"])

    handle = open(path)
    for record in GFF.parse(handle, limit_info=limit_info):
        for feature in record.features:

            try:
                gene_name = "".join(feature.qualifiers["gene"])
            except KeyError:
                gene_name = ""

            info.append(
                [
                    record.id,
                    int(feature.location.start) + 1,
                    int(feature.location.end),
                    1 if feature.location.strand == 1 else 0,
                    gene_name,
                    "".join(feature.qualifiers["locus_tag"]),
                ]
            )
    handle.close()

    columns = ["contig", "start", "end", "strand", "gene_name", "locus_name"]
    parsed_gff = pd.DataFrame(data=info, columns=columns)

    return parsed_gff


def parse_args():
    parser = argparse.ArgumentParser(
        usage="parse_gff.py --input-gff FILE.GFF3 --output-gff PARSED_FILE.TSV",
        description="""Parses .gff3 file with annotation obtained by bakta.""",
    )
    parser.add_argument("-i", "--input-gff", nargs="?", help="gff_file_to_parse.gff3")
    parser.add_argument("-o", "--output-gff", nargs="?", help="parsed_gff_file.tsv")

    return parser.parse_args()


if __name__ == "__main__":
    input_gff = parse_args().input_gff
    output_gff = parse_args().output_gff

    parsed_gff = parse_gff(input_gff)
    parsed_gff.to_csv(output_gff, sep="\t")

    print("Bakta annotation file parsed successfully")
