FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_shotgun_unifrac_env

COPY envs/sbx_shotgun_unifrac_env.yml ./

# Install environment
RUN conda env create --file sbx_shotgun_unifrac_env.yml --name sbx_shotgun_unifrac

ENV PATH="/opt/conda/envs/sbx_shotgun_unifrac/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_shotgun_unifrac", "/bin/bash", "-c"]

# Run
CMD "bash"