#!/bin/bash

print_usage() {
    printf "
Welcome to the additional script for LOFI tool.
This script gets Taxon ID for provided genome FASTA-file needed for next predictions.
Run it if you don't know Taxon ID for your species.

Usage %s:
    -g  <name_of_file>
        genome file name in the LOFI/genomes directory in FASTA format
        only '.fna' extension is supported, if you have '.fa' or another, please change in manually (required).
    -h
        help. Shows this message.

Example usage:
./get_taxid.sh -g GCF_000005845.2_ASM584v2_genomic.fna
" "$0"
}

while getopts ':g:h?' flag; do
    case "${flag}" in
        g)
            snake=true
            genome=$OPTARG
            echo "Using genome file from LOFI/genomes: $genome" >&2
            ;;
        \? | h | *)
            print_usage
            exit 1
            ;;
        esac
    done

if [ "$snake" ] # avoiding positional arguments
then
  genome=${genome%.*} # removes extension to work with snakemake rule
#  snakemake --configfile=config.yaml --cores=all -p genomes/"${genome}"_distances.tab --force # force to recalculate distances
#  snakemake --configfile=config.yaml --cores=all -p genomes/"${genome}"_distances_sorted.tab --force # force to obtain taxid with used genome
  snakemake --cores=all -p genomes/"${genome}"_distances.tab --force # force to recalculate distances
  snakemake --cores=all -p genomes/"${genome}"_distances_sorted.tab --force # force to obtain taxid with used genome
else
  print_usage
  exit 1
fi
