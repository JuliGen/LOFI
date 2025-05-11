#!/bin/bash

print_usage() {
    printf "
Welcome to the additional script for LOFI tool.
This script downloads all DBs and tools required for future operon prediction.
All the steps are included in the main pipeline, but you can use this script to preload everything.

Usage %s:
    This script doesn't require any arguments to work.

    -h
        help. Shows this message.

Example usage:
./download_dbs_tools.sh
" "$0"
}

while getopts 'h?' flag; do
    case "${flag}" in
        \? | h | *)
            print_usage
            exit 1
            ;;
        esac
    done

snakemake --cores=all \
databases/db-light \
databases/ko_list \
databases/profiles/prokaryote.hal \
databases/mash/refseq.genomes.k21s1000.msh \
kofam_scan-1.3.0/exec_annotation

#snakemake --cores=all databases/db-light # download_db_light
#
#snakemake --cores=all mash-Linux64-v2.3 # download_mash
#snakemake --cores=all databases/mash/refseq.genomes.k21s1000.msh # download_mash_db
#
#snakemake --cores=all databases/ko_list # download_db_for_kofam_scan (get_ko_list)
#snakemake --cores=all databases/profiles/prokaryote.hal # download_db_for_kofam_scan (get_profiles)
#snakemake --cores=all kofam_scan-1.3.0/exec_annotation # download_kofam_scan
