def get_shotgun_unifrac_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_shotgun_unifrac":
            return Path(fp)
    raise Error(
        "Filepath for sbx_shotgun_unifrac not found, are you sure it's installed under extensions/sbx_shotgun_unifrac?"
    )


SBX_SHOTGUN_UNIFRAC_VERSION = (
    open(get_shotgun_unifrac_path() / "VERSION").read().strip()
)
SBX_SHOTGUN_UNIFRAC_GG_VERSION = Cfg["sbx_shotgun_unifrac"]["green_genes_version"]

try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


localrules:
    all_shotgun_unifrac,


rule all_shotgun_unifrac:
    input:
        UNIFRAC_FP / "faith.tsv",
        UNIFRAC_FP / "unweighted",
        UNIFRAC_FP / "weighted",


rule su_download_green_genes:
    """Download greengenes db"""
    output:
        dir_fp=dir(Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]),
        seqs_fp=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"bwa.{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.seqs",
    log:
        LOG_FP / "su_download_green_genes.log",
    benchmark:
        BENCHMARK_FP / "su_download_green_genes.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        echo "GG2 db already downloaded" > {log}
        """


rule su_import_green_genes_objects_to_qiime:
    """Probably necessary to create the qza versions of the database objects"""
    input:
        Cfg["sbx_shotgun_unifrac"]["green_genes_fp"],
    output:
        tax=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.taxonomy.id.nwk.qza",
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
    log:
        LOG_FP / "su_import_green_genes_objects_to_qiime.log",
    benchmark:
        BENCHMARK_FP / "su_import_green_genes_objects_to_qiime.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        echo "GG2 db already imported" > {log}
        """


rule su_align_to_green_genes:
    """Align reads to greengenes db"""
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        green_genes_bwa_seqs_fp=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"bwa.{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.seqs",
    output:
        UNIFRAC_FP / "aligned" / "{sample}.sam",
    log:
        LOG_FP / "su_align_to_green_genes.log",
    benchmark:
        BENCHMARK_FP / "su_align_to_green_genes.tsv"
    threads: Cfg["sbx_shotgun_unifrac"]["threads"]
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        "bwa mem -t {threads} {input.green_genes_bwa_seqs_fp} {input.reads} > {output} 2> {log}"


rule su_woltka_classify:
    """Classify reads using woltka"""
    input:
        expand(UNIFRAC_FP / "aligned" / "{sample}.sam", sample=Samples),
    output:
        biom=UNIFRAC_FP / "classified" / "ogu.biom",
        qza=UNIFRAC_FP / "classified" / "ogu.table.qza",
    log:
        LOG_FP / "su_woltka_classify.log",
    benchmark:
        BENCHMARK_FP / "su_woltka_classify.tsv"
    params:
        aligned_fp=UNIFRAC_FP / "aligned",
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        woltka classify -i {params.aligned_fp} -f sam -o {output.biom} 2> {log}
        
        qiime tools import --type FeatureTable[Frequency] \
        --input-path {output.biom} \
        --output-path {output.qza}
        """


rule su_taxonomy_from_table:
    """Get taxonomy from table"""
    input:
        ogu=UNIFRAC_FP / "classified" / "ogu.table.qza",
        tax=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.taxonomy.id.nwk.qza",
    output:
        UNIFRAC_FP / "classified" / "ogu.taxonomy.qza",
    log:
        LOG_FP / "su_taxonomy_from_table.log",
    benchmark:
        BENCHMARK_FP / "su_taxonomy_from_table.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        qiime greengenes2 taxonomy-from-table \
        --i-reference-taxonomy {input.tax} \
        --i-table {input.ogu} \
        --o-classification {output} 2> {log}
        """


rule su_filter_table:
    """Filter table"""
    input:
        ogu=UNIFRAC_FP / "classified" / "ogu.table.qza",
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
    output:
        UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    log:
        LOG_FP / "su_filter_table.log",
    benchmark:
        BENCHMARK_FP / "su_filter_table.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        qiime greengenes2 filter-features \
        --i-feature-table {input.ogu} \
        --i-reference {input.phy} \
        --o-filtered-feature-table {output} 2> {log}
        """


rule su_alpha_phylogenetic_diversity:
    """Calculate alpha phylogenetic diversity using Faith's PD"""
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "faith.qza"),
        tsv=UNIFRAC_FP / "faith.tsv",
    log:
        LOG_FP / "su_alpha_phylogenetic_diversity.log",
    benchmark:
        BENCHMARK_FP / "su_alpha_phylogenetic_diversity.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        qiime diversity alpha-phylogenetic \
        --i-phylogeny {input.phy} \
        --i-table {input.ogu} \
        --p-metric faith_pd \
        --o-alpha-diversity {output.qza}
        
        qiime tools export \
        --input-path {output.qza} \
        --output-path {output.tsv}
        """


rule su_weighted_unifrac_distance:
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "weighted.qza"),
        wu=UNIFRAC_FP / "weighted",
    log:
        LOG_FP / "su_weighted_unifrac_distance.log",
    benchmark:
        BENCHMARK_FP / "su_weighted_unifrac_distance.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        qiime diversity beta-phylogenetic \
        --i-phylogeny {input.phy} \
        --i-table {input.ogu} \
        --p-metric weighted_unifrac \
        --o-distance-matrix {output.qza}
        
        qiime tools export \
        --input-path {output.qza} \
        --output-path {output.tsv}
        """


rule su_unweighted_unifrac_distance:
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "unweighted.qza"),
        uu=UNIFRAC_FP / "unweighted",
    log:
        LOG_FP / "su_unweighted_unifrac_distance.log",
    benchmark:
        BENCHMARK_FP / "su_unweighted_unifrac_distance.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """
        qiime diversity beta-phylogenetic \
        --i-phylogeny {input.phy} \
        --i-table {input.ogu} \
        --p-metric unweighted_unifrac \
        --o-distance-matrix {output.qza}
        
        qiime tools export \
        --input-path {output.qza} \
        --output-path {output.tsv}
        """
