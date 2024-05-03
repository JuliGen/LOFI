import pandas as pd

from sys import argv
from BCBio import GFF
from typing import Any

from metrics.string import get_protein_seqs_links, predict_string
from metrics.intergenic_distances import calculate_intergenic_dist, predict_operon_inter_dist


def parse_gff(path: str) -> tuple[pd.DataFrame, str]:
    """
    Parses file.gff from given path.
    :param path: path to file
    :return: pd.DataFrame with the columns necessary for subsequent analysis; taxon_id for species
    """
    taxon_id = ""
    info = []
    limit_info = dict(gff_type=["gene"])

    handle = open(path)
    for record in GFF.parse(handle, limit_info=limit_info):
        taxon_id = record.annotations["species"][0].split("=")[-1]  # FIXME
        for feature in record.features:
            info.append(
                [
                    record.id,
                    int(feature.location.start) + 1,
                    int(feature.location.end),
                    "+" if feature.location.strand == 1 else "-",
                    "".join(feature.qualifiers["Name"]),
                    "".join(feature.qualifiers["locus_tag"]).replace("_", "")
                ]
            )
    handle.close()

    columns = ["contig", "start", "end", "strand", "gene_name", "locus_name"]
    parsed_gff = pd.DataFrame(data=info, columns=columns)

    return parsed_gff, taxon_id


def final_prediction(*data: Any) -> pd.DataFrame:
    """
    Takes parsed gff table and results from various metrics to predict operonicity for each gene.
    :param data: parsed gff table and metrics
    :return: table with predictions
    """
    df, inter_dist, string = data

    prediction = []
    for n_id in range(df.shape[0]):
        inter_dist_value = inter_dist[n_id]
        string_value = string[n_id]

        if inter_dist_value == string_value:
            if inter_dist_value == 0:
                prediction.append("operon")
            else:
                prediction.append("non_operon")
        elif string_value == 0:  # string > inter_dist_value
            prediction.append("operon")
        else:
            prediction.append("non_operon")

    final_predictions = pd.DataFrame(
        {
            "contig": df.contig,
            "locus_name": df.locus_name,
            "gene_name": df.gene_name,
            "inter_dist": inter_dist,
            "string": string,
            "prediction": prediction,
        }
    )

    return final_predictions


if __name__ == "__main__":
    _, path_to_gff, output_filename = argv

    print("Data processing...")

    # Parsing gff file
    parsed_gff_file, species_id = parse_gff(path_to_gff)

    # Predictions for intergenic distance metric
    inter_dist_df = calculate_intergenic_dist(parsed_gff_file)
    inter_dist_predictions = predict_operon_inter_dist(inter_dist_df)

    # Predictions for STRING metric
    protein_links = get_protein_seqs_links(species_id)
    string_predictions = predict_string(parsed_gff_file, protein_links)

    # Combining all metrics into one table
    df_prediction = final_prediction(parsed_gff_file, inter_dist_predictions, string_predictions)

    df_prediction.to_csv(
        output_filename,
        sep="\t",
        encoding="utf-8",
    )

    print("Job done!")
