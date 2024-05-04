import argparse

import pandas as pd


def combine_predictions(genome, taxon_id, string_scores):  # TODO add other metrics

    final_table = pd.read_csv(f"results/{taxon_id}/bakta/{genome}_parsed.tsv", sep="\t")
    final_table.drop(columns=final_table.columns[0], axis=1, inplace=True)
    string_scores.drop(columns=string_scores.columns[0], axis=1, inplace=True)

    final_table["string_scores"] = string_scores

    return final_table


def parse_args():
    parser = argparse.ArgumentParser(
        usage="main.py --genome GENOME.FNA --taxid TAXON_ID",
        description="""TODO""",
    )
    parser.add_argument("--taxid", nargs="?", help="taxon_id")
    parser.add_argument("--genome", nargs="?", help="genome.fna")

    return parser.parse_args()


if __name__ == "__main__":
    genome = parse_args().genome
    taxon_id = parse_args().taxid

    string_scores = pd.read_csv(
        f"results/{taxon_id}/string/{taxon_id}_{genome}_string_scores.tsv", sep="\t"
    )
    #
    # other metrics
    #
    predictions = combine_predictions(genome, taxon_id, string_scores)
    predictions.to_csv(
        f"results/{taxon_id}/predictions/{genome}_predictions.tsv",
        sep="\t",
        index=False,
    )

    print("Job is done!")
