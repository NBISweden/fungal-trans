FROM continuumio/miniconda3:4.9.2

LABEL maintainer="John Sundh" email=john.sundh@nbis.se
LABEL description="Docker image for a snakemake workflow for fungal metatranscriptomics"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /analysis

# Set tmpdir
ENV TMPDIR="/scratch"

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && apt-get clean

# Add environment file
COPY environment.yml .

# Install environment into base
RUN conda env update -n base -f environment.yml && conda clean -a

# Add workflow
COPY config config
COPY config.yml config.yml
COPY source source
COPY Snakefile Snakefile
COPY envs envs
COPY sample_list.tsv sample_list.tsv

# Run workflow
ENTRYPOINT ["snakemake", "--use-conda", "--conda-frontend mamba", "-j", "1"]
