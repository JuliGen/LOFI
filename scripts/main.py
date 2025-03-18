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


def result_to_gff(path_gff: str, path_tsv: str, output_filename: str) -> None:
    """
    Adds additional qualifier with operon prediction to .gff3 file annotated by bakta.

    :param path_gff: path to file with annotation
    :param path_tsv: path to file with operon prediction results
    :param output_filename: output filename
    """
    df_tsv = pd.read_csv(path_tsv, sep="\t")

    rows_to_keep = list(range(8, df_tsv.shape[0]+9))  #FIXME borders (left=comment, right=comment_length + 1 region)
    df_gff = pd.read_csv(path_gff, sep="\t", header=None, skiprows=lambda x: x not in rows_to_keep)

    for i in range(df_gff.shape[0]):
        if df_gff.iloc[i][2] == "CDS":
            locus_tag = df_gff.iloc[i][8].split(";")[0][3:]
            operon_info = df_tsv.query("locus_name == @locus_tag").operon_number.iloc[0]
            df_gff.iloc[i, 8] = f"{df_gff.iloc[i][8]};operon_info={operon_info}"

    df_gff.to_csv(output_filename, sep="\t", header=False, index=False)


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

    raw_filename_gff = f"results/{taxon_id}/bakta/{genome}.gff3"
    output_filename_tsv = f"results/{taxon_id}/predictions/{genome}_predictions.tsv"
    output_filename_gff = f"results/{taxon_id}/predictions/{genome}_predictions.gff3"

    final_table.to_csv(
        output_filename_tsv,
        sep="\t",
    )

    result_to_gff(raw_filename_gff, output_filename_tsv, output_filename_gff)

    print(f"Job is done! Prediction results are in the folder results/{taxon_id}/predictions/")
