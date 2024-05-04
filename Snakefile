# snakemake --cores=all -p databases/db-light
rule download_db_light:  # DONE
    output:
        directory("databases/db-light")
    shell:
        "bakta_db download --output databases/ --type light"


# snakemake --cores=all -p results/511145/string/511145.protein.links.v12.0.txt
rule download_string_files:  # DONE
    output:  # 4 files in total
        "results/{taxid}/string/{taxid}.protein.links.v12.0.txt",
        "results/{taxid}/string/{taxid}.protein.sequences.v12.0.fa"
    shell:
        "python3 scripts/preprocessing/download_string_files.py --taxid {wildcards.taxid}"


# snakemake --cores=all -p results/511145/diamond/511145.dmnd
rule diamond_makedb:  # DONE
    input:
        "results/{taxid}/string/{taxid}.protein.sequences.v12.0.fa"
    output:
        "results/{taxid}/diamond/{taxid}.dmnd"
    shell:
        "diamond makedb --in {input} --db {output}"


# snakemake --cores=all -p results/511145/bakta/GCF_000005845.2_ASM584v2_genomic.faa
rule bakta_annotation:  # DONE
    input:
        db="databases/db-light",
        genome="data/{genome}.fna"  # files with genome must be in data folder without subfolders
    output:  # 15 files in total
        gff3="results/{taxid}/bakta/{genome}.gff3",
        faa="results/{taxid}/bakta/{genome}.faa"
    params:
        folder="results/{taxid}/bakta/"
    shell:
        "bakta --db {input.db} --output {params.folder} {input.genome} --force"


# snakemake --cores=all -p results/511145/diamond/511145.tsv
rule diamond_blastp:  # DONE
    input:
        faa=rules.bakta_annotation.output.faa,
        db=rules.diamond_makedb.output
    output:
        "results/{taxid}/diamond/{taxid}_{genome}.tsv"  # FIXME
    shell:
        "diamond blastp -q {input.faa} -d {input.db} -o {output} --fast --quiet"


# snakemake --cores=all -p results/511145/bakta/GCF_000005845.2_ASM584v2_genomic_parsed.tsv
rule parse_gff:  # DONE
    input:
        rules.bakta_annotation.output.gff3
    output:
        "results/{taxid}/bakta/{genome}_parsed.tsv"
    shell:
        "python3 scripts/preprocessing/parse_gff.py --input-gff {input} --output-gff {output}"


# snakemake --cores=all -p results/511145/diamond/511145_GCF_000005845.2_ASM584v2_genomic_filtered.tsv
rule filter_diamond_results: # DONE
    input:
        rules.diamond_blastp.output
    output:
        "results/{taxid}/diamond/{taxid}_{genome}_filtered.tsv"
    shell:
        "python3 scripts/preprocessing/filter_diamond_results.py --input {input} --output {output}"


# snakemake --cores=all -p results/511145/predictions/GCF_000005845.2_ASM584v2_genomic_predictions.tsv
rule make_predictions:  # DONE
    input:
        "results/{taxid}/string/{taxid}.protein.links.v12.0.txt",
        "results/{taxid}/bakta/{genome}_parsed.tsv",
        "results/{taxid}/diamond/{taxid}_{genome}_filtered.tsv"
    output:
        "results/{taxid}/predictions/{genome}_predictions.tsv"
    shell:
        "python3 scripts/main.py --genome {wildcards.genome} --taxid {wildcards.taxid}"
