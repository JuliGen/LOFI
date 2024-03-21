import os
import gzip
import shutil
import wget

import numpy as np
import pandas as pd


def get_protein_links(taxon_id: str) -> pd.DataFrame:
    """
    Downloads and unzips a protein score table from the STRING database.

    :param taxon_id: taxon_id for species
    :return: pd.DataFrame with protein combined scores
    """
    suffix = "protein.links.v12.0"
    output_dir = f"results/string_protein_links/{taxon_id}"
    url = f"https://stringdb-downloads.org/download/protein.links.v12.0/{taxon_id}.{suffix}.txt.gz"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        wget.download(url, out=output_dir, bar=None)

        with gzip.open(f"{output_dir}/{taxon_id}.{suffix}.txt.gz", "rb") as in_file:
            with open(f"{output_dir}/{taxon_id}.{suffix}.txt", "wb") as out_file:
                shutil.copyfileobj(in_file, out_file)

    protein_links = pd.read_csv(f"{output_dir}/{taxon_id}.{suffix}.txt", sep=" ")
    protein_links.protein1 = protein_links.protein1.str.replace(taxon_id + ".", "")
    protein_links.protein1 = protein_links.protein1.str.replace("_", "")
    protein_links.protein2 = protein_links.protein2.str.replace(taxon_id + ".", "")
    protein_links.protein2 = protein_links.protein2.str.replace("_", "")

    return protein_links


def get_score_from_df_subset(df_subset: pd.DataFrame, id: str, threshold: int) -> bool:
    """
    Additional function for predict_string function to match
    the id of the second protein in the data subset for the first protein.

    :param df_subset: subset of data from the protein_links table for a specific protein id
    :param id: the second id for the protein that needs to be matched with the first one from df_subset
    :param threshold: score threshold by which a prediction is made
    :return:
    """
    try:
        score = df_subset.query("protein2 == @id").combined_score.iloc[0]
    except IndexError:
        score = np.nan
    return score >= threshold


def predict_string(parsed_gff: pd.DataFrame, protein_links: pd.DataFrame, threshold: int = 810) -> pd.Series:
    """
    Predicts the presence of a gene in an operon based on a score from the STRING database at a given threshold.

    :param parsed_gff: parsed file.gff with annotation: the result of the parse_gff function
    :param protein_links: pd.DataFrame with protein combined scores: the result of get_protein_links function
    :param threshold: threshold for predicting operonicity based on a score from the STRING database
    :return: pd.Series with predictions for each gene
    """
    predictions = []
    df_length = parsed_gff.shape[0]
    threshold = threshold

    for n_id in range(df_length):
        id_cur = parsed_gff.locus_name.iloc[n_id]
        id_prev = parsed_gff.locus_name.iloc[n_id - 1]
        df_subset = protein_links.query("protein1 == @id_cur")

        # 0 - operon, 1 - non_operon
        score_prev = get_score_from_df_subset(df_subset, id_prev, threshold)
        if score_prev:
            operon_predict = 0
        else:
            # Take into account the ring structure of the bacterial chromosome
            if n_id == parsed_gff.shape[0] - 1:
                id_next = parsed_gff.locus_name.iloc[0]
            else:
                id_next = parsed_gff.locus_name.iloc[n_id + 1]

            score_next = get_score_from_df_subset(df_subset, id_next, threshold)
            operon_predict = 0 if score_next else 1

        predictions.append(operon_predict)

    return pd.Series(predictions)
