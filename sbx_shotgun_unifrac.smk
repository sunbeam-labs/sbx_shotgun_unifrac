try:
    SBX_SHOTGUN_UNIFRAC_VERSION = get_ext_version("sbx_shotgun_unifrac")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_SHOTGUN_UNIFRAC_VERSION = "0.0.0"


UNIFRAC_FP = Cfg["all"]["output_fp"] / "shotgun_unifrac"


localrules:
    all_shotgun_unifrac,
    su_extract_outputs,


rule all_shotgun_unifrac:
    input:
        genus=UNIFRAC_FP / "classified" / "genus.tsv",
        phylum=UNIFRAC_FP / "classified" / "phylum.tsv",
        species=UNIFRAC_FP / "classified" / "species.tsv",
        faith=UNIFRAC_FP / "faith_pd_unrarefied.tsv",
        weighted=UNIFRAC_FP / "wu_unrarefied.tsv",
        unweighted=UNIFRAC_FP / "uu_unrarefied.tsv",


rule su_align_to_wolr:
    """Align reads to WoLr db"""
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        sam=temp(UNIFRAC_FP / "aligned" / "{sample}.sam"),
    log:
        LOG_FP / "su_align_to_wolr_{sample}.log",
    benchmark:
        BENCHMARK_FP / "su_align_to_wolr_{sample}.tsv"
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
        (
            bowtie2 \
            -p 8 \
            -x {params.wolr}/WoLr2 \
            -1 {input.r1} \
            -2 {input.r2} \
            --very-sensitive --no-head \
            --no-unal \
            | cut -f1-9 \
            | sed 's/$/\t*\t*/' \
            > {output.sam}
        ) > {log} 2>&1
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
        (
            woltka classify \
            --input {params.fp} \
            -f sam \
            --output {output.biom} \
            --to-biom
        ) > {log} 2>&1
        """


rule su_woltka_classify_map:
    """Classify reads using woltka"""
    input:
        aligned_fp=UNIFRAC_FP / "aligned" / "filtered" / ".done",
    output:
        genus=UNIFRAC_FP / "classified" / "genus.biom",
        phylum=UNIFRAC_FP / "classified" / "phylum.biom",
        species=UNIFRAC_FP / "classified" / "species.biom",
    log:
        LOG_FP / "su_woltka_classify.log",
    benchmark:
        BENCHMARK_FP / "su_woltka_classify.tsv"
    resources:
        runtime=240,
    params:
        fp=UNIFRAC_FP / "aligned" / "filtered",
        map_fp=str(Path(Cfg["sbx_shotgun_unifrac"]["woltka_map_fp"]) / "taxid.map.txt"),
        nodes_fp=str(Path(Cfg["sbx_shotgun_unifrac"]["woltka_map_fp"]) / "nodes.dmp"),
        names_fp=str(Path(Cfg["sbx_shotgun_unifrac"]["woltka_map_fp"]) / "names.dmp"),
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        (
            woltka \
            classify \
            --input {params.fp} \
            --map {params.map_fp} \
            --nodes {params.nodes_fp} \
            --names {params.names_fp} \
            --rank phylum,genus,species \
            --output $(dirname {output.genus}) \
            --to-biom
        ) > {log} 2>&1
        """


rule su_convert_biom_to_tsv:
    input:
        genus=UNIFRAC_FP / "classified" / "genus.biom",
        phylum=UNIFRAC_FP / "classified" / "phylum.biom",
        species=UNIFRAC_FP / "classified" / "species.biom",
    output:
        genus=UNIFRAC_FP / "classified" / "genus.tsv",
        phylum=UNIFRAC_FP / "classified" / "phylum.tsv",
        species=UNIFRAC_FP / "classified" / "species.tsv",
    log:
        LOG_FP / "su_convert_biom_to_tsv.log",
    benchmark:
        BENCHMARK_FP / "su_convert_biom_to_tsv.tsv"
    resources:
        runtime=240,
    conda:
        "envs/sbx_shotgun_unifrac_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        (
            biom convert -i {input.genus} -o {output.genus} --to-tsv
            biom convert -i {input.phylum} -o {output.phylum} --to-tsv
            biom convert -i {input.species} -o {output.species} --to-tsv
        ) > {log} 2>&1
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
        "envs/sbx_q2_diversity_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_shotgun_unifrac:{SBX_SHOTGUN_UNIFRAC_VERSION}"
    shell:
        """
        qiime tools import --type FeatureTable[Frequency] \
        --input-path {input.biom} \
        --output-path {output.qza} > {log} 2>&1
        """


rule su_alpha_phylogenetic_diversity:
    """Calculate alpha phylogenetic diversity using Faith's PD"""
    input:
        phy=Cfg["sbx_shotgun_unifrac"]["tree_fp"],
        ogu=UNIFRAC_FP / "classified" / "ogu.table.qza",
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
        ogu=UNIFRAC_FP / "classified" / "ogu.table.qza",
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
        ogu=UNIFRAC_FP / "classified" / "ogu.table.qza",
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
        faith=temp(UNIFRAC_FP / "faith" / "alpha-diversity.tsv"),
        weighted=temp(UNIFRAC_FP / "weighted" / "distance-matrix.tsv"),
        unweighted=temp(UNIFRAC_FP / "unweighted" / "distance-matrix.tsv"),
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
        (
            qiime tools export \
            --input-path {input.faith} \
            --output-path $(dirname {output.faith})

            qiime tools export \
            --input-path {input.weighted} \
            --output-path $(dirname {output.weighted})

            qiime tools export \
            --input-path {input.unweighted} \
            --output-path $(dirname {output.unweighted})
        ) > {log} 2>&1
        """


rule su_extract_outputs:
    input:
        faith=UNIFRAC_FP / "faith" / "alpha-diversity.tsv",
        weighted=UNIFRAC_FP / "weighted" / "distance-matrix.tsv",
        unweighted=UNIFRAC_FP / "unweighted" / "distance-matrix.tsv",
    output:
        faith=UNIFRAC_FP / "faith_pd_unrarefied.tsv",
        weighted=UNIFRAC_FP / "wu_unrarefied.tsv",
        unweighted=UNIFRAC_FP / "uu_unrarefied.tsv",
    shell:
        """
        mv {input.faith} {output.faith}
        mv {input.weighted} {output.weighted}
        mv {input.unweighted} {output.unweighted}
        """
