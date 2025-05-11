# 👩🏻‍💻LOFI (Locally run Operon Finder based on Integrated metrics)  <a href=""><img src ='https://drive.google.com/uc?export=view&id=1c0OLbJXtGd3ZQvNOSGUGb6lvS2ngdevB' width =300 align="right"></a>
> *This repository is for a project to create a tool for searching for operons*

\
All performance tests were carried out on **Ubuntu 22.04 LTS**, 16 GB RAM, 16 threads.

The tool will take up about **12 GB** of disk space.
<br>
<br>

## Description

This tool is designed to predict the location of bacterial genes in an operon based on 3 metrics: 
**intergenic distance**, **score from the [STRING](https://string-db.org) database** and information about metabolic pathways 
from **[KEGG](https://www.kegg.jp) database**.

## Installation:

Clone this repo and go to created folder: 

```shell
git clone https://github.com/JuliGen/LOFI.git && \
cd LOFI
```

Create mamba environment from `environment.yaml` and activate it:

**Linux:**
```shell
mamba env create -f environment_linux.yaml && \
mamba activate lofi
```
**MacOS:**
```shell
mamba env create -f environment_macos.yaml && \
mamba activate lofi

#Install bakta (Required. You need to install BAKTA separately, because at the moment it is not possible to install it using mamba or conda.
#If you already have BAKTA, you can skip this step.)
git clone https://github.com/oschwengers/bakta.git
cd bakta
python setup.py install
cd ../
```

## Usage

The first run will take longer than subsequent runs due to the loading of the databases required for analysis.

Bash scripts are used to perform all stages of the analysis (downloading the required databases, searching for Taxon ID and predicting operons). However, all the basic commands are implemented using snakemake.

First of all, you need to download all the necessary databases. Before reading further, it is recommended to run this command:

```shell
./download_dbs_tools.sh
```

- In case of errors, simply run the command again. Unfortunately, sometimes there is a problem with the servers, which means you will have to wait until they are fixed at the source.

As input, you must provide a `FASTA file` with the **nucleotide sequence** of the bacterial genome in the `.fna` format (if you have a different extension, please <u>change it manually</u>). Contigs within a FASTA file are expected to be in the **correct order**.

**File with sequence must be in the `genomes` folder**.

As an example of using the tool, we will use a file with the *Escherichia coli K-12* genome, which is already in the genomes folder.

It is also necessary to provide the `taxid` of the species being researched. For proper operation of the tool, this ID must be in the **STRING database**.

If you don't know the `taxid`, you can use the following command, where the input is the same `genome` file as in the main analysis:

```shell
./get_taxid.sh -g <genome>
```

To start the prediction, use the command below:

```shell
./start_prediction.sh -t <taxid> -g <genome> -f
```

The `-f` (force) flag is used to overwrite the prediction.

### Example usage

The [_E. coli_ K-12](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/511145/) genome, which is already uploaded to the `genomes` folder, is used as an example.

Downloading DBs:

```shell
./download_dbs_tools.sh
```

Obtaining Taxon ID:

```shell
./get_taxid.sh -g GCF_000005845.2_ASM584v2_genomic.fna
```

Predicting operons:

```shell
./start_prediction.sh -t 511145 -g GCF_000005845.2_ASM584v2_genomic.fna -f
```

- If you want to launch **kofam_scan** in parallel while doing pipeline, go to the appropriate rule in `Snakefile` and change the `params.threads` variable to the desired one. Default value is 8.
- If you want to change the email to your own (currently it is set to a generic one), go to the `Snakefile` and change the global variable at the beginning.

## Troubleshooting

Make sure **taxid** is in the file `data/species.v12.0.txt`, that is, it is in the STRING database.

Make sure that the **genome** of the species is located in the `genomes` folder.

If you encounter any issues or have questions, ~~try to use _E. coli_ K-12~~ directly contact the authors for support.

## Authors

* Yulia Nechaeva (Perm State University)
* Artem Vasiliev (Saint Petersburg State University)

*Supervisors*: Oksana Kotovskaya (Skoltech), Nikita Vaulin (Skoltech, vaulin@ro.ru), Nadezhda Pavlova (Lomonosov Moscow State University)

Initial idea was proposed by Polina Kuchur and Alexey Komissarov (ITMO) and researhed by Anna Churkina and Anna Rybina in [[O-antigens](https://github.com/rybinaanya/O-antigens)]
