import argparse
import pickle

import pandas as pd


def predict_operon(model_file: str, data_for_predict: pd.DataFrame) -> pd.Series:
    """
    Makes predictions using a trained model.

    :param data_for_predict: table with data required for prediction
    :param model_file: path to model
    :return: series with predicted operons
    """
    with open(model_file, "rb") as f:
        model = pickle.load(f)

    predictions = model.predict(data_for_predict)
    return pd.Series(predictions)


def parse_args():
    parser = argparse.ArgumentParser(
        usage="predict_operon.py"
        "--parsed-gff PARSED_GFF.TSV"
        "--string STRING_SCORES.TSV"
        "--inter-dist INTER_DIST.TSV"
        "--kegg KEGG.TSV"
        "--model MODEL.PKL"
        "--output RESULT.TSV",
        description="""Makes predictions using a trained RandomForest model.""",
    )
    parser.add_argument("--parsed-gff", nargs="1", help="path to parsed gff file.tsv")
    parser.add_argument("--string", nargs="1", help="path to STRING scores.tsv")
    parser.add_argument(
        "--inter-dist", nargs="1", help="path to intergenic distances result file.tsv"
    )
    parser.add_argument("--kegg", nargs="1", help="path to kegg result file.npy")
    parser.add_argument("--model", nargs="1", help="model for predictions.pkl")
    parser.add_argument("-o", "--output", nargs="1", help="predictions results.tsv")
    return parser.parse_args()


if __name__ == "__main__":
    path_parsed_gff = parse_args().parsed_gff
    path_string = parse_args().string
    path_inter_dist = parse_args().inter_dist
    path_kegg = parse_args().kegg
    model = parse_args().model
    output_filename = parse_args().output

    parsed_gff = pd.read_csv(path_parsed_gff, sep="\t")
    string = pd.read_csv(path_string, sep="\t")
    inter_dist = pd.read_csv(path_inter_dist, sep="\t")
    kegg = pd.read_csv(path_kegg, sep="\t")

    data_for_predict = pd.DataFrame(
        data={
            "strand": parsed_gff.strand,
            "prob_operon": inter_dist.iloc[:, 1],
            "pred_string": string.final_scores,
            "intersection_map_count": kegg.intersection_map_count,
        }
    )

    result = predict_operon(model, data_for_predict)
    result.to_csv(output_filename, sep="\t")

    print("Operons predicted successfully")
