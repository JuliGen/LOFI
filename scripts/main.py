import pandas as pd

from BCBio import GFF


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
                    "".join(feature.qualifiers["gene"]),
                    feature.id.split("-")[1],
                ]
            )
    handle.close()

    columns = ["contig", "start", "end", "strand", "gene_name", "locus_name"]
    parsed_gff = pd.DataFrame(data=info, columns=columns)

    return parsed_gff, taxon_id
