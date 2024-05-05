import argparse
import pickle

import pandas as pd


def predict_operon(model_file: str, data_for_predict: pd.DataFrame) -> pd.Series:
    """
    :param data_for_predict:
    :param model_file: path to model
    :return:pd.Series with predicted operons
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
        description="""TODO""",
    )
    parser.add_argument("--parsed-gff", nargs="?", help="path to parsed gff file.tsv")
    parser.add_argument("--string", nargs="?", help="path to STRING scores.tsv")
    parser.add_argument(
        "--inter-dist", nargs="?", help="path to intergenic distances result file.tsv"
    )
    parser.add_argument("--kegg", nargs="?", help="path to kegg result file.npy")
    parser.add_argument("--model", nargs="?", help="model for predictions.pkl")
    parser.add_argument("-o", "--output", nargs="?", help="predictions results.tsv")
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
