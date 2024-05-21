import argparse

import numpy as np
import pandas as pd


def format_output(
    genome: str, taxon_id: int, predictions: pd.DataFrame
) -> pd.DataFrame:
    """
    Creates final table with results of operons predictions.

    :param genome: input genome fasta file
    :param taxon_id: NCBI taxon id
    :param predictions: file with predictions
    :return: table with final results
    """
    parsed_gff = pd.read_csv(f"results/{taxon_id}/bakta/{genome}_parsed.tsv", sep="\t")
    parsed_gff.drop(columns=parsed_gff.columns[0], axis=1, inplace=True)
    predictions.drop(columns=predictions.columns[0], axis=1, inplace=True)

    final_table = parsed_gff.copy()
    kegg_annotation = pd.read_csv(
        f"results/{taxon_id}/predictions/temp_dir/{taxon_id}_{genome}_kegg.tsv",
        sep="\t",
    )
    string_result = pd.read_csv(
        f"results/{taxon_id}/predictions/temp_dir/{taxon_id}_{genome}_string_scores.tsv",
        sep="\t",
    )

    final_table["kegg_orthology"] = kegg_annotation["kegg_orthology"]
    final_table["metabolic_pathway_kegg"] = kegg_annotation[
        "metabolic_pathway_kegg"
    ].replace(np.nan, "")
    final_table["description_kegg"] = kegg_annotation["description_kegg"]
    final_table["prediction"] = predictions

    current_operon = 1
    dict_operons = {1: []}
    for i in range(len(final_table) - 1):
        if final_table["prediction"].iloc[i] == 1:
            dict_operons[current_operon].append(i)
            map_current_gene = final_table["metabolic_pathway_kegg"].iloc[i]
            map_next_gene = final_table["metabolic_pathway_kegg"].iloc[i + 1]
            string_next_gene = string_result["next_scores"].iloc[i]
            if (
                string_next_gene < 810
                and len(
                    set(map_current_gene.split(",")).intersection(
                        map_next_gene.split(",")
                    )
                )
                < 1
            ) or final_table["prediction"].iloc[i + 1] == 0:
                current_operon += 1
                dict_operons[current_operon] = []

    final_table.insert(loc=0, column="operon_number", value="non_operon")

    for operon_number, index_cds in dict_operons.items():
        final_table.loc[index_cds, "operon_number"] = f"operon_{operon_number}"

    final_table = final_table.drop(["prediction", "strand"], axis=1)

    return final_table


def parse_args():
    parser = argparse.ArgumentParser(
        usage="main.py --genome GENOME.FNA --taxid TAXON_ID",
        description="""Starts operon prediction pipeline.""",
    )
    parser.add_argument("--genome", nargs="?", help="genome.fna")
    parser.add_argument("--taxid", nargs="?", help="ncbi taxon id")

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
