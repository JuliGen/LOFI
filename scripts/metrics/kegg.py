import argparse
import json
from typing import Tuple, Dict, List

import pandas as pd


def get_ko_map(
    path_hmm_result: str, path_to_ko_map: str = "data/ko_map.json"
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Gets information about KO and metabolic pathways for each protein_id
    :param path_to_ko_map:
    :param path_hmm_result: path to hmm result file (after kofam scan)
    :return:
    dict_ko: dict with KO for each protein_id
    dict_map: dict with metabolic pathways for each protein_id
    """

    with open(path_to_ko_map, "r", encoding="utf-8") as fh:
        ko_map = json.load(fh)
    dict_ko = {}
    dict_map = {}
    with open(path_hmm_result, "r") as file:
        for line in file:
            map = []
            protein_id = line.strip().split("\t")[0]
            kos = line.strip().split("\t")[1:]
            if not kos:
                map.extend(["NaN"])
            else:
                for ko in kos:
                    if not ko_map[ko]:
                        map.extend(["NaN"])
                    map.extend(ko_map[ko])
            dict_ko[protein_id] = kos
            dict_map[protein_id] = list(set(map))
    return dict_ko, dict_map


def calc_intersection_map(parsed_gff: pd.DataFrame, path_hmm_result: str) -> pd.Series:
    """
    Counts the number of matches between metabolic pathways
    :param parsed_gff: parsed file.gff with annotation: the result of the parse_gff function
    :return: pd.Series with the number of intersections between metabolic pathways
    """

    dict_ko, dict_map = get_ko_map(path_hmm_result)
    parsed_gff["ko"] = parsed_gff["locus_name"].map(dict_ko)
    parsed_gff["map"] = parsed_gff["locus_name"].map(dict_map)

    map_list = parsed_gff.dropna()["map"].to_list()
    map_list_sh1 = map_list[1:]
    map_list_sh1.append(["NaN"])
    map_list_sh2 = map_list[:-1]
    map_list_sh2.insert(0, ["NaN"])

    inter_map = []
    for list1, list2 in zip(map_list, map_list_sh1):
        inter_map.append(len(set(list1).intersection(list2)))
    inter_map2 = []
    for list1, list2 in zip(map_list, map_list_sh2):
        inter_map2.append(len(set(list1).intersection(list2)))

    parsed_gff.loc[parsed_gff.dropna()["map"].index, "intersection_map_count"] = [
        max(itm) for itm in zip(inter_map, inter_map2)
    ]

    return parsed_gff.intersection_map_count


def parse_args():
    parser = argparse.ArgumentParser(
        usage="kegg.py --input_gff PATH_TO_PARSED_GFF_FILE.TSV --input_hmm PATH_TO_HMM_RESULT.txt --output INTERSECTION_MAP_COUNT.tsv",  # TODO extension
        description="""TODO""",
    )
    parser.add_argument("-igff", "--input-gff", nargs="?", help="parsed gff file.tsv")
    parser.add_argument(
        "-ihmm", "--input-hmm", nargs="?", help="hmm results after kofam scan.tsv"
    )
    parser.add_argument("-o", "--output", nargs="?", help="path to hmm result file.txt")

    return parser.parse_args()


if __name__ == "__main__":
    path_parsed_gff = parse_args().input_gff
    path_hmm_result = parse_args().input_hmm
    output_filename = parse_args().output

    parsed_gff = pd.read_csv(path_parsed_gff, sep="\t")
    result = calc_intersection_map(parsed_gff, path_hmm_result)
    result = result.fillna(0)

    result.to_csv(output_filename, sep="\t")

    print("Kegg analysis completed")
