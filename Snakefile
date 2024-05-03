rule download_db_light:
    output:
        "databases/db-light"
    shell:
        "bakta_db download --output {output} --type light"


# Как вариант разбить эту часть на 2: скачивание 2-х нужных файлов через терминал и просто чтение в string.py, и реализация main.py уже здесь
rule get_protein_seqs_links:
    output:
        "results/{taxid}/string/{taxid}.protein.links.v12.0.txt",
        "results/{taxid}/string/{taxid}.protein.sequences.v12.0.fa"
    shell:
        "python3 scripts/main.py --genome {wildcards.genome}.fna --taxid {wildcards.taxid}"


rule diamond_makedb:
    output:
        "results/{taxid}/diamond/{taxid}.dmnd"
    input:
        "results/{taxid}/string/{taxid}.protein.sequences.v12.0.fa"
    shell:
        "diamond makedb --in {input} --db {output}"


rule bakta_annotation:
    input:
        "{genome}.fna"
    output:
        gff3="results/{taxid}/bakta/{genome}.gff3",
        faa="results/{taxid}/bakta/{genome}.faa"
    params:
        db="databases/db-light",
        folder="results/{taxid}/bakta/"
    shell:
        "bakta --db {params.db} --output {params.folder} {input}.fna --force"


rule diamond_blastp:
    input:
        rules.bakta_annotation.output.faa
    output:
        "results/{taxid}/diamond/{taxid}.tsv"
    params:
        db="results/{taxid}/diamond/{taxid}.dmnd"
    shell:
        "-q {input} -d {params.db} -o {output} --fast --quiet"
