# Comments represent single launch of rule

# configfile: "config.yaml"
# print("Config is: ", config)

EMAIL = "****@example.com"


# snakemake --cores=all -p databases/db-light
rule download_db_light:
    output:
        directory("databases/db-light")
    shell:
        "bakta_db download --output databases/ --type light"


rule download_db_for_kofam_scan:
    output:
        "databases/ko_list.gz",
        "databases/profiles.tar.gz"
    run:
        shell("gdown --fuzzy \
        https://drive.google.com/file/d/1csxtWD2UBbs6XAwTCVEqzGVVDzJKQhOl/view\?usp\=drive_link \
        -O databases/ko_list.gz")
        shell("gdown --fuzzy \
        https://drive.google.com/file/d/1m_907sx0XhL5UC8pTYg24b4caA4OZ_Cx/view?usp=drive_link \
        -O databases/profiles.tar.gz")


# snakemake --cores=all -p databases/ko_list
rule get_ko_list:
    input:
        "databases/ko_list.gz"
    output:
        "databases/ko_list"
    shell:
        "gunzip {input} -c > {output} && rm {input}"


# snakemake --cores=all -p databases/profiles
rule get_profiles:
    input:
        "databases/profiles.tar.gz"
    output:
        "databases/profiles/prokaryote.hal"
    shell:
        "tar -xzf {input} -C databases/ && rm {input}"


# snakemake --cores=all -p results/511145/string/511145.protein.links.v12.0.txt
rule download_string_files:
    output:  # 4 files in total
        "results/{taxid}/string/{taxid}.protein.links.v12.0.txt",
        "results/{taxid}/string/{taxid}.protein.sequences.v12.0.fa"
    shell:
        "python3 scripts/preprocessing/download_string_files.py --taxid {wildcards.taxid}"


# # snakemake --cores=all -p mash-Linux64-v2.3
# rule download_mash:
#     output:
#         directory("mash-Linux64-v2.3")
#     run:
#         shell("wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar")
#         shell("tar -xf mash-Linux64-v2.3.tar && rm mash-Linux64-v2.3.tar")


# snakemake --cores=all -p databases/mash/refseq.genomes.k21s1000.msh
rule download_mash_db:
    output:
        "databases/mash/refseq.genomes.k21s1000.msh"
    shell:
        """
        wget -P databases/mash/ \
        https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh \
        --tries=10 \
        --retry-connrefused
        """


# snakemake --cores=all -p genomes/GCF_000005845.2_ASM584v2_genomic_distances.tab
rule run_mash_dist:
    input:
        genome="genomes/{genome}.fna", #FIXME
        ref=rules.download_mash_db.output
    output:
        "genomes/{genome}_distances.tab"
    shell:
        "mash dist {input.ref} {input.genome} > {output}"


# snakemake --cores=all -p genomes/GCF_000005845.2_ASM584v2_genomic_distances_sorted.tab
rule obtain_taxid:
    input:
        genome="genomes/{genome}.fna",
        dist="genomes/{genome}_distances.tab"
    output:
        "genomes/{genome}_distances_sorted.tab"
    shell:
        """
        python3 scripts/preprocessing/obtain_taxid.py \
        --genome {input.genome} \
        --dist {input.dist} \
        --email {EMAIL}
        """


# snakemake --cores=all -p results/511145/diamond/511145.dmnd
rule diamond_makedb:
    input:
        "results/{taxid}/string/{taxid}.protein.sequences.v12.0.fa"
    output:
        "results/{taxid}/diamond/{taxid}.dmnd"
    shell:
        "diamond makedb --in {input} --db {output}"


# snakemake --cores=all -p results/511145/bakta/GCF_000005845.2_ASM584v2_genomic.faa
rule bakta_annotation:
    input:
        db="databases/db-light",
        genome="genomes/{genome}.fna"  # files with genome must be in genomes folder without subfolders
    output:  # 15 files in total
        gff3="results/{taxid}/bakta/{genome}.gff3",
        faa="results/{taxid}/bakta/{genome}.faa"
    params:
        folder="results/{taxid}/bakta/"
    shell:
        """
        bakta --db {input.db} --output {params.folder} {input.genome} --force \
        --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region \
        --skip-crispr --skip-pseudo --skip-sorf --skip-gap --skip-ori --skip-plot
        """


# snakemake --cores=all -p results/511145/diamond/511145.tsv
rule diamond_blastp:
    input:
        faa=rules.bakta_annotation.output.faa,
        db=rules.diamond_makedb.output
    output:
        "results/{taxid}/diamond/{taxid}_{genome}.tsv"  # FIXME
    shell:
        "diamond blastp -q {input.faa} -d {input.db} -o {output} --fast --quiet"


# snakemake --cores=all -p results/511145/bakta/GCF_000005845.2_ASM584v2_genomic_parsed.tsv
rule parse_gff:
    input:
        rules.bakta_annotation.output.gff3
    output:
        "results/{taxid}/bakta/{genome}_parsed.tsv"
    shell:
        "python3 scripts/preprocessing/parse_gff.py --input-gff {input} --output-gff {output}"


# snakemake --cores=all -p results/511145/diamond/511145_GCF_000005845.2_ASM584v2_genomic_filtered.tsv
rule filter_diamond_results:
    input:
        rules.diamond_blastp.output
    output:
        "results/{taxid}/diamond/{taxid}_{genome}_filtered.tsv"
    shell:
        "python3 scripts/preprocessing/filter_diamond_results.py --input {input} --output {output}"


# snakemake --cores=all -p results/511145/predictions/temp_dir/511145_GCF_000005845.2_ASM584v2_genomic_string_scores.tsv
rule get_string_scores:
    input:
        parsed_gff="results/{taxid}/bakta/{genome}_parsed.tsv",
        filtered_diamond_result="results/{taxid}/diamond/{taxid}_{genome}_filtered.tsv",
        protein_links="results/{taxid}/string/{taxid}.protein.links.v12.0.txt"
    output:
        "results/{taxid}/predictions/temp_dir/{taxid}_{genome}_string_scores.tsv"
    shell:
        """
        python3 scripts/metrics/get_string_scores.py \
        --parsed-gff {input.parsed_gff} \
        --filtered-diamond-result {input.filtered_diamond_result} \
        --protein-links {input.protein_links} \
        --output {output}
        """


# snakemake --cores=all -p kofam_scan-1.3.0/exec_annotation
rule download_kofam_scan:
    output:
        "kofam_scan-1.3.0/exec_annotation"
    run:
        shell(
            """
            wget https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz \
            --tries=10 \
            --retry-connrefused
            """
        )
        shell("tar -xzf kofam_scan-1.3.0.tar.gz && rm kofam_scan-1.3.0.tar.gz")


# snakemake --cores=all -p results/511145/hmm/511145_GCF_000005845.2_ASM584v2_genomic_hmm.txt
rule kofam_scan:
    input:
        kofam=rules.download_kofam_scan.output,
        faa=rules.bakta_annotation.output.faa,
        ko_list="databases/ko_list",
        profile="databases/profiles/prokaryote.hal"
    output:
        "results/{taxid}/hmm/{taxid}_{genome}_hmm.txt"
    params:
        format="mapper-one-line",
        threads=8
    shell:
        """
        {input.kofam} \
        --cpu={params.threads} \
        -p {input.profile} \
        -k {input.ko_list} \
        -f {params.format} \
        -o {output} \
        {input.faa}
        """


# snakemake --cores=all -p results/511145/predictions/temp_dir/511145_GCF_000005845.2_ASM584v2_genomic_kegg.tsv
rule kegg:
    input:
        gff=rules.parse_gff.output,
        hmm=rules.kofam_scan.output
    output:
        "results/{taxid}/predictions/temp_dir/{taxid}_{genome}_kegg.tsv"
    shell:
        "python3 scripts/metrics/kegg.py --input-gff {input.gff} --input-hmm {input.hmm} --output {output}"


# snakemake --cores=all -p results/511145/predictions/temp_dir/511145_GCF_000005845.2_ASM584v2_inter_dist.tsv
rule intergenic_distances:
    input:
        rules.parse_gff.output
    output:
        "results/{taxid}/predictions/temp_dir/{taxid}_{genome}_inter_dist.tsv"
    params:
        matrix="data/matrix_emission_15.npy"
    shell:
        """
        python3 scripts/metrics/intergenic_distances.py \
        --input {input} \
        --emission-matrix {params.matrix} \
        --output {output}
        """


# snakemake --cores=all -p results/511145/predictions/GCF_000005845.2_ASM584v2_genomic_predictions.tsv
rule predict_operons:
    input:
        gff=rules.parse_gff.output,
        string=rules.get_string_scores.output,
        inter_dist=rules.intergenic_distances.output,
        kegg=rules.kegg.output
    output:
        "results/{taxid}/predictions/temp_dir/{taxid}_{genome}_predictions.tsv"
    params:
        model="models/model.pkl"
    shell:
        """
        python3 scripts/metrics/predict_operon.py \
        --parsed-gff {input.gff} \
        --string {input.string} \
        --inter-dist {input.inter_dist} \
        --kegg {input.kegg} \
        --model {params.model} \
        --output {output}
        """


# snakemake --cores=all -p results/511145/predictions/GCF_000005845.2_ASM584v2_genomic_predictions.tsv
rule main:
    input:
        rules.parse_gff.output,
        rules.predict_operons.output
    output:
        "results/{taxid}/predictions/{genome}_predictions.tsv"
    shell:
        "python3 scripts/main.py --genome {wildcards.genome} --taxid {wildcards.taxid}"
