import os
import gzip
import shutil
import wget

import numpy as np
import pandas as pd


def get_protein_seqs_links(taxon_id: int) -> pd.DataFrame:
    """
    Downloads and unzips a protein score table and protein sequences
    for provided taxon_id from the STRING database.

    :param taxon_id: taxon_id for species
    :return: pd.DataFrame with protein combined scores
    """
    taxon_id = taxon_id

    db_version = "v12.0"
    keywords = [["sequences", "fa"], ["links", "txt"]]
    output_dir = f"results/{taxon_id}/string"

    urls = []

    for keyword in keywords:
        urls.append(
            f"https://stringdb-downloads.org/download/protein.{keyword[0]}.{db_version}/"
            f"{taxon_id}.protein.{keyword[0]}.{db_version}.{keyword[1]}.gz"
        )

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

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

    protein_links = pd.read_csv(
        f"{output_dir}/{taxon_id}.protein.links.{db_version}.txt", sep=" "
    )

    return protein_links


def filter_diamond_results(path: str):
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


def predict_string(
    parsed_gff: pd.DataFrame,
    diamond_result_filtered: pd.DataFrame,
    protein_links: pd.DataFrame,
    threshold: int = 810,
) -> pd.Series:
    """
    Predicts the presence of a gene in an operon based on a score from the STRING database at a given threshold.

    :param parsed_gff: parsed file.gff with annotation: the result of the parse_gff function
    :param diamond_result_filtered:
    :param protein_links: pd.DataFrame with protein combined scores: the result of get_protein_links function
    :param threshold: threshold for predicting operonicity based on a score from the STRING database
    :return: pd.Series with predictions for each gene
    """
    # TODO пока что предполагается, что контиги в геномной последовательности расположены в правильном порядке
    predictions = []
    df_length = parsed_gff.shape[0]

    for n_id in range(df_length):
        id_cur = parsed_gff.locus_name.iloc[n_id]
        try:
            id_cur_string = diamond_result_filtered.query(
                "query_accession == @id_cur"
            ).target_accession.iloc[0]
        except IndexError:
            predictions.append(0)
            continue

        df_subset = protein_links.query("protein1 == @id_cur_string")

        id_prev = parsed_gff.locus_name.iloc[n_id - 1]
        try:
            id_prev_string = diamond_result_filtered.query(
                "query_accession == @id_prev"
            ).target_accession.iloc[0]
        except IndexError:
            score_prev = False  #TODO
        else:
            score_prev = get_score_from_df_subset(df_subset, id_prev_string, threshold)

        # 1 - operon, 0 - non_operon
        if score_prev:
            operon_predict = 1
        else:
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
                score_next = False
            else:
                score_next = get_score_from_df_subset(
                    df_subset, id_next_string, threshold
                )
            operon_predict = 1 if score_next else 0

        predictions.append(operon_predict)

    return pd.Series(predictions)
