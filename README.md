<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_shotgun_unifrac

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_shotgun_unifrac/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_shotgun_unifrac/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_shotgun_unifrac)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_shotgun_unifrac/)
<!-- Badges end -->

## Introduction

sbx_shotgun_unifrac is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for assigning shotgun reads to a bacterial tree of life and performing population statistics like UniFrac. This pipeline uses [bwa](https://bio-bwa.sourceforge.net/) for initial alignment steps and [QIIME2](https://qiime2.org/) with the [GreenGenes plugin](https://forum.qiime2.org/t/introducing-greengenes2-2022-10/25291) for the reference phylogeny and performing the core methods.

## Config

  - threads: Number of threads to use with `bwa`
  - wolr_fp: Directory with WoLr2 bowtie2 database
  - tree_fp: Qiime2 phylogeny
  - woltka_map_fp: Path to the directory containing the Woltka map file

## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).