import argparse

import numpy as np
import pandas as pd


def format_output(genome, taxon_id, predictions):  # TODO add other metrics

    final_table = pd.read_csv(f"results/{taxon_id}/bakta/{genome}_parsed.tsv", sep="\t")
    final_table.drop(columns=final_table.columns[0], axis=1, inplace=True)
    predictions.drop(columns=predictions.columns[0], axis=1, inplace=True)

    final_table["prediction"] = predictions
    final_table.prediction = final_table.prediction.replace(
        [1, 0], ["operon", "non_operon"]
    )

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

    predictions = pd.read_csv(
        f"results/{taxon_id}/predictions/temp_dir/{taxon_id}_{genome}_predictions.tsv",
        sep="\t",
    )
    final_table = format_output(genome, taxon_id, predictions)
    output_filename = f"results/{taxon_id}/predictions/{genome}_predictions.tsv"
    final_table.to_csv(
        output_filename,
        sep="\t",
    )

    print(f"Job is done! Prediction results are in the folder {output_filename}")
