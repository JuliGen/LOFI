import argparse

import pandas as pd


def filter_diamond_results(path: str) -> pd.DataFrame:
    """
    Filters results after Diamond alignment.

    :param path: file path to diamond alignment
    :return: filtered table
    """
    column_names = [
        "query_accession",
        "target_accession",
        "sequence_identity",
        "length",
        "mismatches",
        "gap_openings",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "e_value",
        "bit_score",
    ]
    diamond_result = pd.read_csv(path, sep="\t", header=None, names=column_names)
    diamond_result_filtered = diamond_result.loc[
        diamond_result.groupby("query_accession")["sequence_identity"].idxmax()
    ]

    return diamond_result_filtered


def parse_args():
    parser = argparse.ArgumentParser(
        usage="filter_diamond_results.py --input DIAMOND_RESULT.TSV --output FILTERED_DIAMOND_RESULT.TSV",
        description="""Filters results after Diamond alignment.""",
    )
    parser.add_argument("-i", "--input", nargs="?", help="diamond result file.tsv")
    parser.add_argument(
        "-o", "--output", nargs="?", help="filtered diamond result file.tsv"
    )

    return parser.parse_args()


if __name__ == "__main__":
    diamond_result = parse_args().input
    output_filename = parse_args().output

    diamond_result_filtered = filter_diamond_results(diamond_result)
    diamond_result_filtered.to_csv(output_filename, sep="\t")

    print("Diamond results were filtered successfully")
