import argparse
import json
from typing import Tuple, Dict, List

import pandas as pd


def get_ko_map(
    path_hmm_result: str, path_to_ko_map: str = "data/ko_map.json"
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Gets information about KO and metabolic pathways for each protein_id.

    :param path_to_ko_map: .json KO map file
    :param path_hmm_result: path to hmm result file (after kofam scan)
    :return: dict with KO for each protein_id, dict with metabolic pathways for each protein_id
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
                map.extend([])
            else:
                for ko in kos:
                    if not ko_map[ko]:
                        map.extend([])
                    map.extend(ko_map[ko])
            dict_ko[protein_id] = kos
            dict_map[protein_id] = list(set(map))
    return dict_ko, dict_map


def calc_intersection_map(
    parsed_gff: pd.DataFrame,
    path_hmm_result: str,
    path_to_ko_desc: str = "data/ko_descriptions.json",
) -> pd.DataFrame:
    """
    Counts the number of matches between metabolic pathways.

    :param parsed_gff: parsed file.gff with annotation: the result of the parse_gff function
    :param path_hmm_result: path to hmm result file (after kofam scan)
    :param path_to_ko_desc: .json file with KO descriptions
    :return: table with the number of intersections between metabolic pathways
    """
    dict_ko, dict_map = get_ko_map(path_hmm_result)
    parsed_gff["ko"] = parsed_gff["locus_name"].map(dict_ko)
    parsed_gff["map"] = parsed_gff["locus_name"].map(dict_map)

    map_list = parsed_gff.dropna()["map"].to_list()
    map_list_sh1 = map_list[1:]
    map_list_sh1.append([])
    map_list_sh2 = map_list[:-1]
    map_list_sh2.insert(0, [])

    inter_map = []
    for list1, list2 in zip(map_list, map_list_sh1):
        inter_map.append(len(set(list1).intersection(list2)))
    inter_map2 = []
    for list1, list2 in zip(map_list, map_list_sh2):
        inter_map2.append(len(set(list1).intersection(list2)))

    parsed_gff.loc[parsed_gff.dropna()["map"].index, "intersection_map_count"] = [
        max(itm) for itm in zip(inter_map, inter_map2)
    ]

    with open(path_to_ko_desc, "r", encoding="utf-8") as fh:
        ko_desc = json.load(fh)

    list_desc_ko = []
    ko_description = []
    for kos in parsed_gff["ko"]:
        for ko in kos:
            ko_description.append("".join(ko_desc[ko]))
        list_desc_ko.append("; ".join(ko_description))
        ko_description = []

    kegg_annotation = pd.DataFrame(
        data={
            "kegg_orthology": [",".join(map(str, l)) for l in parsed_gff["ko"]],
            "metabolic_pathway_kegg": [
                ",".join(map(str, l)) for l in parsed_gff["map"]
            ],
            "description_kegg": list_desc_ko,
            "intersection_map_count": parsed_gff.intersection_map_count,
        }
    )

    kegg_annotation = kegg_annotation.fillna("")

    return kegg_annotation


def parse_args():
    parser = argparse.ArgumentParser(
        usage="kegg.py \
        --input_gff PARSED_GFF.TSV \
        --input_hmm PATH_TO_HMM_RESULT.TXT \
        --output INTERSECTION_MAP_COUNT.TSV",
        description="""Gets intersection of metabolic pathways for each protein.""",
    )
    parser.add_argument("-igff", "--input-gff", nargs="1", help="parsed gff file.tsv")
    parser.add_argument(
        "-ihmm", "--input-hmm", nargs="1", help="hmm results after kofam scan.tsv"
    )
    parser.add_argument("-o", "--output", nargs="1", help="path to hmm result file.txt")

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
