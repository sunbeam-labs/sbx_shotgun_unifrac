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
    su_temp_install_pip,
    su_extract_outputs,


rule all_shotgun_unifrac:
    input:
        faith=UNIFRAC_FP / "faith_pd_unrarefied.tsv",
        weighted=UNIFRAC_FP / "wu_unrarefied.tsv",
        unweighted=UNIFRAC_FP / "uu_unrarefied.tsv",


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


rule su_align_to_wolr:
    """Align reads to WoLr db"""
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
        wolr=expand(
            (Path(Cfg["sbx_shotgun_unifrac"]["wolr_fp"]) / "WoLr2").with_suffix(
                "{ext}"
            ),
            ext=[
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            ],
        ),
        pip=UNIFRAC_FP / ".pip_installed",
    output:
        sam=temp(UNIFRAC_FP / "aligned" / "{sample}.sam"),
    log:
        LOG_FP / "su_align_to_green_genes_{sample}.log",
    benchmark:
        BENCHMARK_FP / "su_align_to_green_genes_{sample}.tsv"
    params:
        wolr=Cfg["sbx_shotgun_unifrac"]["wolr_fp"],
    threads: Cfg["sbx_shotgun_unifrac"]["threads"]
    resources:
        runtime=240,
        mem_mb=100000,
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        bowtie2 -p 8 -x {params.wolr} \
        -1 {input.r1} \
        -2 {input.r2} \
        --very-sensitive --no-head \
        --no-unal | cut -f1-9 | sed 's/$/\t*\t*/' > {output.sam} 2> {log}
        """


rule su_filter_on_sequence_number:
    input:
        expand(UNIFRAC_FP / "aligned" / "{sample}.sam", sample=Samples),
    output:
        temp(UNIFRAC_FP / "aligned" / "filtered" / ".done"),
    log:
        LOG_FP / "su_filter_on_sequence_number.log",
    benchmark:
        BENCHMARK_FP / "su_filter_on_sequence_number.tsv"
    params:
        fp=UNIFRAC_FP / "aligned" / "filtered",
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        echo "Files with 100 or fewer sequences:" > {log}
        mkdir -p {params.fp}
        for f in {input} ; do
            if [ $(samtools view -c -F 4 $f) -le 100 ]; then
                echo $f >> {log}
            else
                cp $f {params.fp}
            fi
        done
        touch {output}
        """


rule su_woltka_classify:
    """Classify reads using woltka"""
    input:
        aligned_fp=UNIFRAC_FP / "aligned" / "filtered" / ".done",
    output:
        biom=UNIFRAC_FP / "classified" / "ogu.biom",
    log:
        LOG_FP / "su_woltka_classify.log",
    benchmark:
        BENCHMARK_FP / "su_woltka_classify.tsv"
    resources:
        runtime=240,
    params:
        fp=UNIFRAC_FP / "aligned" / "filtered",
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        woltka classify -i {params.fp} -f sam -o {output.biom} > {log} 2>&1
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


rule su_filter_table:
    """Filter table"""
    input:
        ogu=UNIFRAC_FP / "classified" / "ogu.table.qza",
        phy=Cfg["sbx_shotgun_unifrac"]["tree_fp"],
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
        phy=Cfg["sbx_shotgun_unifrac"]["tree_fp"],
        ogu=UNIFRAC_FP / "classified" / "ogu.filtered.table.qza",
    output:
        qza=temp(UNIFRAC_FP / "faith_pd_vector.qza"),
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
        phy=Cfg["sbx_shotgun_unifrac"]["tree_fp"],
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
        --p-metric weighted_normalized_unifrac \
        --o-distance-matrix {output.qza} > {log} 2>&1
        """


rule su_unweighted_unifrac_distance:
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["tree_fp"],
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
        faith=UNIFRAC_FP / "faith_pd_vector.qza",
        weighted=UNIFRAC_FP / "weighted.qza",
        unweighted=UNIFRAC_FP / "unweighted.qza",
    output:
        faith=temp(UNIFRAC_FP / "faith"),
        weighted=temp(UNIFRAC_FP / "weighted"),
        unweighted=temp(UNIFRAC_FP / "unweighted"),
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


rule su_extract_outputs:
    input:
        faith=UNIFRAC_FP / "faith",
        weighted=UNIFRAC_FP / "weighted",
        unweighted=UNIFRAC_FP / "unweighted",
    output:
        faith=UNIFRAC_FP / "faith_pd_unrarefied.tsv",
        weighted=UNIFRAC_FP / "wu_unrarefied.tsv",
        unweighted=UNIFRAC_FP / "uu_unrarefied.tsv",
    shell:
        """
        mv {input.faith} / alpha-diversity.tsv {output.faith}
        mv {input.weighted} / distance-matrix.tsv {output.weighted}
        mv {input.unweighted} / distance-matrix.tsv {output.unweighted}
        """
