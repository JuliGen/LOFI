#!/bin/bash

print_usage() {
    printf "
Welcome to the main script for LOFI tool.
This script starts pipeline to predict operons in the specified genome file.

Usage %s:
    -t  <taxon_id>
        taxon ID for specified genome FASTA file (required).
    -g  <name_of_file>
        genome file name in the LOFI/genomes directory in FASTA format
        only '.fna' extension is supported, if you have '.fa' or another, please change in manually (required).
    -f  <force>
        use if you want to rewrite results (optional).
    -h
        help. Shows this message.

Example usage:
./start_prediction.sh -t 511145 -g GCF_000005845.2_ASM584v2_genomic.fna -f
" "$0"
}

while getopts ':t:g:fh?' flag; do
    case "${flag}" in
        t)
            snake_taxid=true
            taxid=$OPTARG
            echo "Using this Taxon ID: $taxid" >&2
            ;;
        g)
            snake_genome=true
            genome=$OPTARG
            echo "Using genome file from LOFI/genomes: $genome" >&2
            ;;
        f)
            force_run=true
            ;;
        \? | h | *)
            print_usage
            exit 1
            ;;
        esac
    done

if [ "$snake_taxid" ] && [ "$snake_genome" ] # avoiding positional arguments
then
  genome=${genome%.*} # removes extension to work with snakemake rule
  if [ "$force_run" ]
    then
      snakemake --cores=all -p results/"${taxid}"/predictions/"${genome}"_predictions.tsv --force # TODO all pipeline?
    else
      snakemake --cores=all -p results/"${taxid}"/predictions/"${genome}"_predictions.tsv
  fi
else
  print_usage
  exit 1
fi
