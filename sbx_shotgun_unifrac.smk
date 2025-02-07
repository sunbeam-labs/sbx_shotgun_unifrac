def get_shotgun_unifrac_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_shotgun_unifrac":
            return Path(fp)
    raise Error(
        "Filepath for sbx_shotgun_unifrac not found, are you sure it's installed under extensions/sbx_shotgun_unifrac?"
    )


UNIFRAC_FP = Cfg["all"]["output_fp"] / "shotgun_unifrac"
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
        UNIFRAC_FP / "faith",
        UNIFRAC_FP / "unweighted",
        UNIFRAC_FP / "weighted",


# rule su_download_green_genes:
#     """Download greengenes db"""
#     output:
#         phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
#         / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk",
#         seqs_fp=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
#         / f"bwa.{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.seqs",
#     log:
#         LOG_FP / "su_download_green_genes.log",
#     benchmark:
#         BENCHMARK_FP / "su_download_green_genes.tsv"
#     conda:
#         "envs/sbx_shotgun_unifrac_env.yml"
#     container:
#         f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
#     shell:
#         """
#         echo "RULE NOT IMPLEMENTED, ASSUMING PREEXISTING DB" > {log}
#         """


# rule su_import_green_genes_objects_to_qiime:
#     """Probably necessary to create the qza versions of the database objects"""
#     input:
#         Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
#         / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk",
#     output:
#         phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
#         / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
#     log:
#         LOG_FP / "su_import_green_genes_objects_to_qiime.log",
#     benchmark:
#         BENCHMARK_FP / "su_import_green_genes_objects_to_qiime.tsv"
#     conda:
#         "envs/sbx_shotgun_unifrac_env.yml"
#     container:
#         f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
#     shell:
#         """
#         echo "RULE NOT IMPLEMENTED, ASSUMING PREEXISTING IMPORTS" > {log}
#         """


rule su_temp_install_pip:
    """TEMPORARY: install pip packages because the conda file can't handle it"""
    output:
        temp(UNIFRAC_FP / ".pip_installed"),
    log:
        LOG_FP / "su_temp_install_pip.log",
    benchmark:
        BENCHMARK_FP / "su_temp_install_pip.tsv"
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        ${{CONDA_PREFIX}}/bin/python -m pip install numpy > {log}
        ${{CONDA_PREFIX}}/bin/python -m pip install cython >> {log}
        ${{CONDA_PREFIX}}/bin/python -m pip install q2-greengenes2 >> {log}
        touch {output}
        """


rule su_align_to_green_genes:
    """Align reads to greengenes db"""
    input:
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        green_genes_bwa_seqs_fp=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"bwa.{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.seqs",
        green_genes_bwa_seqs_indexes_fp=[
            Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
            / f"bwa.{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.seqs.{ext}"
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
        ],
        pip=UNIFRAC_FP / ".pip_installed",
    output:
        temp(UNIFRAC_FP / "aligned" / "{sample}.sam"),
    log:
        LOG_FP / "su_align_to_green_genes_{sample}.log",
    benchmark:
        BENCHMARK_FP / "su_align_to_green_genes_{sample}.tsv"
    threads: Cfg["sbx_shotgun_unifrac"]["threads"]
    resources:
        runtime=240,
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        "bwa mem -t {threads} {input.green_genes_bwa_seqs_fp} {input.reads} > {output} 2> {log}"


rule su_woltka_classify:
    """Classify reads using woltka"""
    input:
        expand(UNIFRAC_FP / "aligned" / "{sample}.sam", sample=Samples),
    output:
        biom=UNIFRAC_FP / "classified" / "ogu.biom",
    log:
        LOG_FP / "su_woltka_classify.log",
    benchmark:
        BENCHMARK_FP / "su_woltka_classify.tsv"
    params:
        aligned_fp=UNIFRAC_FP / "aligned",
    resources:
        runtime=240,
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        woltka classify -i {params.aligned_fp} -f sam -o {output.biom} > {log} 2>&1
        """


rule su_convert_biom_to_qza:
    input:
        biom=UNIFRAC_FP / "classified" / "ogu.biom",
    output:
        qza=UNIFRAC_FP / "classified" / "ogu.table.qza",
    log:
        LOG_FP / "su_convert_biom_to_qza.log",
    benchmark:
        BENCHMARK_FP / "su_convert_biom_to_qza.tsv"
    resources:
        runtime=240,
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime tools import --type FeatureTable[Frequency] \
        --input-path {input.biom} \
        --output-path {output.qza} > {log} 2>&1
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
    resources:
        runtime=240,
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime greengenes2 taxonomy-from-table \
        --i-reference-taxonomy {input.tax} \
        --i-table {input.ogu} \
        --o-classification {output} > {log} 2>&1
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
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime greengenes2 filter-features \
        --i-feature-table {input.ogu} \
        --i-reference {input.phy} \
        --o-filtered-feature-table {output} > {log} 2>&1
        """


rule su_alpha_phylogenetic_diversity:
    """Calculate alpha phylogenetic diversity using Faith's PD"""
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "faith.qza"),
    log:
        LOG_FP / "su_alpha_phylogenetic_diversity.log",
    benchmark:
        BENCHMARK_FP / "su_alpha_phylogenetic_diversity.tsv"
    conda:
        "envs/sbx_q2_diversity_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime diversity alpha-phylogenetic \
        --i-phylogeny {input.phy} \
        --i-table {input.ogu} \
        --p-metric faith_pd \
        --o-alpha-diversity {output.qza} > {log} 2>&1
        """


rule su_weighted_unifrac_distance:
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "weighted.qza"),
    log:
        LOG_FP / "su_weighted_unifrac_distance.log",
    benchmark:
        BENCHMARK_FP / "su_weighted_unifrac_distance.tsv"
    conda:
        "envs/sbx_q2_diversity_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime diversity beta-phylogenetic \
        --i-phylogeny {input.phy} \
        --i-table {input.ogu} \
        --p-metric weighted_unifrac \
        --o-distance-matrix {output.qza} > {log} 2>&1
        """


rule su_unweighted_unifrac_distance:
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]
        / f"{SBX_SHOTGUN_UNIFRAC_GG_VERSION}.phylogeny.id.nwk.qza",
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "unweighted.qza"),
    log:
        LOG_FP / "su_unweighted_unifrac_distance.log",
    benchmark:
        BENCHMARK_FP / "su_unweighted_unifrac_distance.tsv"
    conda:
        "envs/sbx_q2_diversity_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime diversity beta-phylogenetic \
        --i-phylogeny {input.phy} \
        --i-table {input.ogu} \
        --p-metric unweighted_unifrac \
        --o-distance-matrix {output.qza} > {log} 2>&1
        """


rule su_export_qzas:
    input:
        faith=UNIFRAC_FP / "faith.qza",
        weighted=UNIFRAC_FP / "weighted.qza",
        unweighted=UNIFRAC_FP / "unweighted.qza",
    output:
        faith=UNIFRAC_FP / "faith",
        weighted=UNIFRAC_FP / "weighted",
        unweighted=UNIFRAC_FP / "unweighted",
    log:
        LOG_FP / "su_export_qzas.log",
    benchmark:
        BENCHMARK_FP / "su_export_qzas.tsv"
    conda:
        "envs/sbx_q2_diversity_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime tools export \
        --input-path {input.faith} \
        --output-path {output.faith} > {log} 2>&1

        qiime tools export \
        --input-path {input.weighted} \
        --output-path {output.weighted} >> {log} 2>&1

        qiime tools export \
        --input-path {input.unweighted} \
        --output-path {output.unweighted} >> {log} 2>&1
        """
