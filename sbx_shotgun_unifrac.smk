def get_shotgun_unifrac_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_shotgun_unifrac":
            return Path(fp)
    raise Error(
        "Filepath for sbx_shotgun_unifrac not found, are you sure it's installed under extensions/sbx_shotgun_unifrac?"
    )


SBX_SHOTGUN_UNIFRAC_VERSION = open(get_shotgun_unifrac_path() / "VERSION").read().strip()

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
        UNIFRAC_FP / "faith" / "{sample}",
        UNIFRAC_FP / "unweighted" / "{sample}",
        UNIFRAC_FP / "weighted" / "{sample}",


rule su_download_green_genes:
    """Download greengenes db"""
    output:
        dir(Cfg["sbx_shotgun_unifrac"]["green_genes_fp"]),
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
        
        """


rule su_align_to_green_genes:
    """Align reads to greengenes db"""
    input:
        reads=expand(QC_FP / "cleaned" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        green_genes_fp=Cfg["sbx_shotgun_unifrac"]["green_genes_fp"],
    output:
        UNIFRAC_FP / "aligned" / "{sample}.sam",
    log:
        LOG_FP / "su_align_to_green_genes.log",
    benchmark:
        BENCHMARK_FP / "su_align_to_green_genes.tsv"
    params:
        green_genes_version=Cfg["sbx_shotgun_unifrac"]["green_genes_version"],
    threads: Cfg["sbx_shotgun_unifrac"]["threads"],
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        "bwa mem -t {threads} {input.green_genes_fp}/bwa.{params.green_genes_version}.seqs {input.reads} > {output} 2> {log}"


rule su_alpha_phylogenetic_diversity:
    """Calculate alpha phylogenetic diversity using Faith's PD"""
    input:

    output:
        qza=temp(UNIFRAC_FP / "faith" / "{sample}.qza"),
        tsv=UNIFRAC_FP / "faith" / "{sample}.tsv",
    log:
        LOG_FP / "su_alpha_phylogenetic_diversity.log",
    benchmark:
        BENCHMARK_FP / "su_alpha_phylogenetic_diversity.tsv"
    params:
        green_genes_version=Cfg["sbx_shotgun_unifrac"]["green_genes_version"],
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """qiime diversity alpha-phylogenetic \
        --i-phylogeny "{input.green_genes_fp}/{params.green_genes_version}.phylogeny.id.nwk.qza" \
        --i-table "${DIR_OUT}/ogu.filtered.table.qza" \
        --p-metric faith_pd \
        --o-alpha-diversity {output.qza}
        
        qiime tools export \
        --input-path {output.qza} \
        --output-path {output.tsv}
        """

rule su_weighted_unifrac_distance:
    input:

    output:
        qza=temp(UNIFRAC_FP / "weighted" / "{sample}.qza"),
        tsv=UNIFRAC_FP / "weighted" / "{sample}.tsv",
    log:
        LOG_FP / "su_weighted_unifrac_distance.log",
    benchmark:
        BENCHMARK_FP / "su_weighted_unifrac_distance.tsv"
    params:
        green_genes_version=Cfg["sbx_shotgun_unifrac"]["green_genes_version"],
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """qiime diversity beta-phylogenetic \
        --i-phylogeny "{input.green_genes_fp}/{params.green_genes_version}.phylogeny.id.nwk.qza" \
        --i-table "${DIR_OUT}/ogu.filtered.table.qza" \
        --p-metric weighted_unifrac \
        --o-distance-matrix {output.qza}
        
        qiime tools export \
        --input-path {output.qza} \
        --output-path {output.tsv}
        """


rule su_unweighted_unifrac_distance:
    input:

    output:
        qza=temp(UNIFRAC_FP / "unweighted" / "{sample}.qza"),
        tsv=UNIFRAC_FP / "unweighted" / "{sample}.tsv",
    log:
        LOG_FP / "su_unweighted_unifrac_distance.log",
    benchmark:
        BENCHMARK_FP / "su_unweighted_unifrac_distance.tsv"
    params:
        green_genes_version=Cfg["sbx_shotgun_unifrac"]["green_genes_version"],
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_TEMPLATE_VERSION}"
    shell:
        """qiime diversity beta-phylogenetic \
        --i-phylogeny "{input.green_genes_fp}/{params.green_genes_version}.phylogeny.id.nwk.qza" \
        --i-table "${DIR_OUT}/ogu.filtered.table.qza" \
        --p-metric unweighted_unifrac \
        --o-distance-matrix {output.qza}
        
        qiime tools export \
        --input-path {output.qza} \
        --output-path {output.tsv}
        """