import argparse

import pandas as pd


def make_predictions():
    data = [1, 2, 3]
    predictions = pd.DataFrame(data=data)

    return predictions


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

    protein_links = pd.read_csv(
        f"results/{taxon_id}/string/{taxon_id}.protein.links.v12.0.txt"
    )
    parsed_gff = pd.read_csv(f"results/{taxon_id}/bakta/{genome}_parsed.tsv")
    diamond_result_filtered = pd.read_csv(
        f"results/{taxon_id}/diamond/{taxon_id}_{genome}_filtered.tsv"
    )

    predictions = make_predictions()
    predictions.to_csv(
        f"results/{taxon_id}/predictions/{genome}_predictions.tsv", sep="\t"
    )

    print("Job is done!")
