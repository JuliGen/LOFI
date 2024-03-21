# 👩🏻‍💻OperonMapper <img src ='https://papik.pro/uploads/posts/2021-12/1639240390_33-papik-pro-p-dinozavr-klipart-33.png' width =400 align="right">
> *This repository is for a project to create a tool for searching for operons*

## Warning

Unfortunately, the tool <u>does not currently work</u> on the **macOS** operating system.

All performance tests were carried out on **Ubuntu 22.04 LTS**.

## Description

This tool is designed to predict the location of bacterial genes in an operon based on 2 metrics: 
**intergenic distance** and a **score from the STRING database** (the number of metrics will increase soon).

At this stage of development, <wonderful name of the tool> is working with fully sequenced and annotated bacterial
species that are in the [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/) database, and at the same time 
are in the [STRING](https://string-db.org/) database.

## Installation:

Clone this repo and go to created folder: 

```shell
git clone git@github.com:JuliGen/OperonMapper.git && \
cd OperonMapper
```

Create mamba environment from `environment.yaml` and activate it:

```shell
mamba env create -f environment.yaml && \
mamba activate operonmapper
```

## Usage

```shell
python3 scripts/main.py \
path_to_annotation.gff \
output_filename.csv
```

Supported filetypes for annotation: `.gff`, `.gff3`

### Example usage

```shell
python3 scripts/main.py \
data/example/E_coli_K_12_from_NCBI.gff \
results/E_coli_predicted_operons.csv
```

## Troubleshooting

If you encounter any issues or have questions, ~~try to use _E. coli K-12_~~ directly contact the authors for support.
