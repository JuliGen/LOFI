import argparse

import numpy as np
import pandas as pd


def get_score_from_df_subset(df_subset: pd.DataFrame, id: str) -> int:
    """
    Additional function for predict_string function to match
    the id of the second protein in the data subset for the first protein.

    :param df_subset: subset of data from the protein_links table for a specific protein id
    :param id: the second id for the protein that needs to be matched with the first one from df_subset
    :return:
    """
    try:
        score = df_subset.query("protein2 == @id").combined_score.iloc[0]
    except IndexError:
        score = np.nan
    return score


def get_string_scores(
        parsed_gff: pd.DataFrame,
        diamond_result_filtered: pd.DataFrame,
        protein_links: pd.DataFrame,
) -> pd.DataFrame:
    """
    Predicts the presence of a gene in an operon based on a score from the STRING database at a given threshold.

    :param parsed_gff: parsed file.gff with annotation: the result of the parse_gff function
    :param diamond_result_filtered:
    :param protein_links: pd.DataFrame with protein combined scores: the result of get_protein_links function
    :return: pd.Series with predictions for each gene
    """
    # FIXME (for now contigs must be in the correct order)
    final_scores = []
    next_scores = []
    prev_scores = []

    df_length = parsed_gff.shape[0]

    for n_id in range(df_length):
        id_cur = parsed_gff.locus_name.iloc[n_id]
        try:
            id_cur_string = diamond_result_filtered.query(
                "query_accession == @id_cur"
            ).target_accession.iloc[0]
        except IndexError:
            final_scores.append(500)  # mean threshold for prediction
            next_scores.append(500)
            prev_scores.append(500)
            continue

        df_subset = protein_links.query("protein1 == @id_cur_string")

        id_prev = parsed_gff.locus_name.iloc[n_id - 1]
        try:
            id_prev_string = diamond_result_filtered.query(
                "query_accession == @id_prev"
            ).target_accession.iloc[0]
        except IndexError:
            score_prev = 500
        else:
            score_prev = get_score_from_df_subset(df_subset, id_prev_string)

        # Take into account the ring structure of the bacterial chromosome
        if n_id == parsed_gff.shape[0] - 1:
            id_next = parsed_gff.locus_name.iloc[0]
        else:
            id_next = parsed_gff.locus_name.iloc[n_id + 1]

        try:
            id_next_string = diamond_result_filtered.query(
                "query_accession == @id_next"
            ).target_accession.iloc[0]
        except IndexError:
            score_next = 500
        else:
            score_next = get_score_from_df_subset(df_subset, id_next_string)

        final_score = max(score_prev, score_next)
        final_scores.append(final_score)
        next_scores.append(score_next)
        prev_scores.append(score_prev)

    scores = pd.DataFrame(
        data={
            "prev_scores": prev_scores,
            "next_scores": next_scores,
            "final_scores": final_scores
        }
    )
    scores = scores.fillna(500)

    return scores


def parse_args():
    parser = argparse.ArgumentParser(
        usage="get_string_scores.py \
        --parsed-gff PARSED_FILE.TSV \
        --filtered-diamond-result FILTERED_DIAMOND_RESULT.TSV \
        --protein-links PROTEIN_LINKS.TXT \
        --output STRING_SCORES.TSV",
        description="""Parses .gff3 file with annotation obtained by bakta.""",
    )
    parser.add_argument("--parsed-gff", nargs="?", help="parsed gff file.tsv")
    parser.add_argument(
        "--filtered-diamond-result", nargs="?", help="filtered diamond result.tsv"
    )
    parser.add_argument(
        "--protein-links", nargs="?", help="protein links from STRING db.txt"
    )
    parser.add_argument(
        "-o", "--output", nargs="?", help="result file with obtained STRING scores.txt"
    )

    return parser.parse_args()


if __name__ == "__main__":
    path_parsed_gff = parse_args().parsed_gff
    path_filtered_diamond_result = parse_args().filtered_diamond_result
    path_protein_links = parse_args().protein_links
    output_filename = parse_args().output

    parsed_gff = pd.read_csv(path_parsed_gff, sep="\t")
    filtered_diamond_result = pd.read_csv(path_filtered_diamond_result, sep="\t")
    protein_links = pd.read_csv(path_protein_links, sep=" ")

    string_scores = get_string_scores(
        parsed_gff, filtered_diamond_result, protein_links
    )

    string_scores.to_csv(output_filename, sep="\t")
    print("STRING scores obtained")
