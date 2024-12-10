<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_shotgun_unifrac

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_shotgun_unifrac/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_shotgun_unifrac/actions/workflows/tests.yml)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_shotgun_unifrac)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_shotgun_unifrac/)
<!-- Badges end -->

## Introduction

sbx_shotgun_unifrac is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for assigning shotgun reads to a bacterial tree of life and performing population statistics like UniFrac. This pipeline uses [bwa](https://bio-bwa.sourceforge.net/) for initial alignment steps and [QIIME2](https://qiime2.org/) with the [GreenGenes plugin](https://forum.qiime2.org/t/introducing-greengenes2-2022-10/25291) for the reference phylogeny and performing the core methods.

## Installation

Extension install is as simple as passing the extension's URL on GitHub to `sunbeam extend`:

    sunbeam extend https://github.com/sunbeam-labs/sbx_shotgun_unifrac

Any user-modifiable parameters specified in `config.yml` are automatically added on `sunbeam init`. If you're installing an extension in a project where you already have a config file, run the following to add the options for your newly added extension to your config (the `-i` flag means in-place config file modification; remove the `-i` flag to see the new config in stdout):

    sunbeam config update -i /path/to/project/sunbeam_config.yml

Installation instructions for older versions of Sunbeam are included at the end of this README.

## Running

To run an extension, simply run Sunbeam as usual with your extension's target rule specified:

    sunbeam run --profile /path/to/project/ all_shotgun_unifrac

### Options for config.yml

  - threads: Number of threads to use with `bwa`
  - green_genes_fp: Directory with all GreenGenes reference files
  - green_genes_version: Version of the GreenGenes db (E.g. "2024.09")
    
## Installing an extension (legacy instructions for sunbeam <3.0)

Installing an extension is as simple as cloning (or moving) your extension directory into the sunbeam/extensions/ folder, installing requirements through Conda, and adding the new options to your existing configuration file: 

    git clone https://github.com/sunbeam-labs/sbx_shotgun_unifrac/ sunbeam/extensions/sbx_shotgun_unifrac
    cat sunbeam/extensions/sbx_shotgun_unifrac/config.yml >> sunbeam_config.yml

## Issues with pipeline

Please post any issues with this extension [here](https://github.com/sunbeam-labs/sbx_shotgun_unifrac/issues).
