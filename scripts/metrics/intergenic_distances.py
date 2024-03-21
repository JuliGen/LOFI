import numpy as np
import pandas as pd
from hmmlearn import hmm


def calculate_intergenic_dist(parsed_gff: pd.DataFrame) -> pd.DataFrame:
    """
    Counts the distance between genes.
    :param parsed_gff: parsed gff file
    :return: table with additional columns (intergenic_distance_next - the distance to the next gene,
    intergenic_distance_prev - the distance to the previous gene)
    """
    df_inter_dist = pd.DataFrame(columns=["intergenic_distance_next", "intergenic_distance_prev"])

    # Add columns with intergenic distance
    df_inter_dist["intergenic_distance_next"] = (parsed_gff.start - parsed_gff.end.shift(1)).shift(-1)
    df_inter_dist["intergenic_distance_prev"] = df_inter_dist["intergenic_distance_next"].shift(1)

    # Find gene indexes at the beginning and end of the contigs
    ind_start_contig = parsed_gff.loc[parsed_gff["contig"] != parsed_gff["contig"].shift(1)].index
    ind_end_contig = parsed_gff.loc[parsed_gff["contig"] != parsed_gff["contig"].shift(-1)].index

    # FIXME
    df_inter_dist.loc[ind_start_contig, "intergenic_distance_prev"] = 1
    df_inter_dist.loc[ind_end_contig, "intergenic_distance_next"] = 1

    # FIXME
    # Replace the negative distance
    df_inter_dist["intergenic_distance_next"] = np.where(
        df_inter_dist.intergenic_distance_next < 0,
        0,
        df_inter_dist.intergenic_distance_next,
    )

    # Dividing the observed states into categories
    df_inter_dist["cat_dist"] = pd.cut(
        df_inter_dist["intergenic_distance_next"],
        [i for i in range(-1, 800, 15)] + [10000],
        labels=[i for i in range(54)],
    )

    return df_inter_dist


def predict_operon_inter_dist(df_inter_dist: pd.DataFrame) -> pd.Series:
    """
    Calculates the оperon score for genes based on the intergenic distance.
    :param df_inter_dist: pd.DataFrame
    :return: operonic predictions for each gene
    """

    # Predict operon
    model = hmm.CategoricalHMM(n_components=2, algorithm="map")  # build model (lib hmm learn), forward-backward
    model.startprob_ = np.array([0.5, 0.5])  # initialize start probability
    model.transmat_ = np.array([[0.72143634, 0.27856366], [0.19284369, 0.80715631]])  # initialize transition matrix
    model.emissionprob_ = np.load("data/matrix_emission_15.npy")  # initialize emission matrix

    obs_states = np.array(df_inter_dist.cat_dist.values).reshape(-1, 1)
    hid_states = pd.Series(model.predict(obs_states))

    return hid_states
