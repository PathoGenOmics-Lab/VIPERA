FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7071b22b1161190c06be3ac061ab0019a1a8d3038532c9134e070a3875414ef5"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/afwdist.yaml
#   prefix: /conda-envs/9c24a867826615972cc288081976e7fc
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - afwdist==1.0.0
RUN mkdir -p /conda-envs/9c24a867826615972cc288081976e7fc
COPY workflow/envs/afwdist.yaml /conda-envs/9c24a867826615972cc288081976e7fc/environment.yaml

# Conda environment:
#   source: workflow/envs/biopython.yaml
#   prefix: /conda-envs/162796cecea22d99c8702138f0c48e2f
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python==3.10
#     - biopython==1.81
#     - pandas==2.0.3
RUN mkdir -p /conda-envs/162796cecea22d99c8702138f0c48e2f
COPY workflow/envs/biopython.yaml /conda-envs/162796cecea22d99c8702138f0c48e2f/environment.yaml

# Conda environment:
#   source: workflow/envs/fetch.yaml
#   prefix: /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - entrez-direct==16.2
#     - curl
RUN mkdir -p /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a
COPY workflow/envs/fetch.yaml /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a/environment.yaml

# Conda environment:
#   source: workflow/envs/freyja.yaml
#   prefix: /conda-envs/bb4c5f3a509433cc08861582fab4a705
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - freyja==2.0.1
RUN mkdir -p /conda-envs/bb4c5f3a509433cc08861582fab4a705
COPY workflow/envs/freyja.yaml /conda-envs/bb4c5f3a509433cc08861582fab4a705/environment.yaml

# Conda environment:
#   source: workflow/envs/gisaidr.yaml
#   prefix: /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f
#   channels:
#     - conda-forge
#   dependencies:
#     - r-base=4.1.3
#     - r-tidyverse==2.0.0
#     - r-devtools==2.4.5
#     - r-logger==0.2.2
RUN mkdir -p /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f
COPY workflow/envs/gisaidr.yaml /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f/environment.yaml
COPY workflow/envs/gisaidr.post-deploy.sh /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f/post-deploy.sh

# Conda environment:
#   source: workflow/envs/iqtree.yaml
#   prefix: /conda-envs/0a608afb24723cb6fa8aef748f5efbc8
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - iqtree==2.2.2.3
RUN mkdir -p /conda-envs/0a608afb24723cb6fa8aef748f5efbc8
COPY workflow/envs/iqtree.yaml /conda-envs/0a608afb24723cb6fa8aef748f5efbc8/environment.yaml

# Conda environment:
#   source: workflow/envs/nextalign.yaml
#   prefix: /conda-envs/04a3347f94ddf7e21c34bc49e5246076
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - nextalign==2.13
RUN mkdir -p /conda-envs/04a3347f94ddf7e21c34bc49e5246076
COPY workflow/envs/nextalign.yaml /conda-envs/04a3347f94ddf7e21c34bc49e5246076/environment.yaml

# Conda environment:
#   source: workflow/envs/pangolin.yaml
#   prefix: /conda-envs/fd645c541ee7a3d43fb9167441b77888
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - pangolin==4.3.1
RUN mkdir -p /conda-envs/fd645c541ee7a3d43fb9167441b77888
COPY workflow/envs/pangolin.yaml /conda-envs/fd645c541ee7a3d43fb9167441b77888/environment.yaml

# Conda environment:
#   source: workflow/envs/quarto_render.yaml
#   prefix: /conda-envs/96f3c1cec4b3ce5d72f708992272e9c1
#   channels:
#     - conda-forge
#   dependencies:
#     - r-base==4.5.2
#     - r-gt==1.1.0
#     - quarto==1.8.25
#     - deno==2.3.1
#     - r-tidyverse==2.0.0
#     - r-heatmaply==1.6.0
RUN mkdir -p /conda-envs/96f3c1cec4b3ce5d72f708992272e9c1
COPY workflow/envs/quarto_render.yaml /conda-envs/96f3c1cec4b3ce5d72f708992272e9c1/environment.yaml

# Conda environment:
#   source: workflow/envs/renv.yaml
#   prefix: /conda-envs/8ad6cdcf265d30289788da99d5bf9fff
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base=4.3.2
#     - r-tidyverse==2.0.0
#     - r-ggrepel==0.9.3
#     - r-ggpubr==0.6.0
#     - bioconductor-ggtree==3.10.0
#     - r-ape==5.8
#     - r-adephylo==1.1_13
#     - r-pegas==1.2
#     - r-data.table==1.14.8
#     - r-future.apply==1.11.0
#     - r-scales==1.3.0
#     - r-showtext==0.9_6
#     - r-logger==0.2.2
RUN mkdir -p /conda-envs/8ad6cdcf265d30289788da99d5bf9fff
COPY workflow/envs/renv.yaml /conda-envs/8ad6cdcf265d30289788da99d5bf9fff/environment.yaml

# Conda environment:
#   source: workflow/envs/snpeff.yaml
#   prefix: /conda-envs/0adafb79cb1bec58ef4c77bf4cca4f95
#   channels:
#     - bioconda
#   dependencies:
#     - snpeff==5.1d
#     - snpsift==5.1d
RUN mkdir -p /conda-envs/0adafb79cb1bec58ef4c77bf4cca4f95
COPY workflow/envs/snpeff.yaml /conda-envs/0adafb79cb1bec58ef4c77bf4cca4f95/environment.yaml

# Conda environment:
#   source: workflow/envs/var_calling.yaml
#   prefix: /conda-envs/5150d0f0a91d7f7a789a06f453d63479
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ivar==1.4.2
#     - samtools==1.17
RUN mkdir -p /conda-envs/5150d0f0a91d7f7a789a06f453d63479
COPY workflow/envs/var_calling.yaml /conda-envs/5150d0f0a91d7f7a789a06f453d63479/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/9c24a867826615972cc288081976e7fc --file /conda-envs/9c24a867826615972cc288081976e7fc/environment.yaml && \
    conda env create --prefix /conda-envs/162796cecea22d99c8702138f0c48e2f --file /conda-envs/162796cecea22d99c8702138f0c48e2f/environment.yaml && \
    conda env create --prefix /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a --file /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a/environment.yaml && \
    conda env create --prefix /conda-envs/bb4c5f3a509433cc08861582fab4a705 --file /conda-envs/bb4c5f3a509433cc08861582fab4a705/environment.yaml && \
    conda env create --prefix /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f --file /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f/environment.yaml && \
    conda env create --prefix /conda-envs/0a608afb24723cb6fa8aef748f5efbc8 --file /conda-envs/0a608afb24723cb6fa8aef748f5efbc8/environment.yaml && \
    conda env create --prefix /conda-envs/04a3347f94ddf7e21c34bc49e5246076 --file /conda-envs/04a3347f94ddf7e21c34bc49e5246076/environment.yaml && \
    conda env create --prefix /conda-envs/fd645c541ee7a3d43fb9167441b77888 --file /conda-envs/fd645c541ee7a3d43fb9167441b77888/environment.yaml && \
    conda env create --prefix /conda-envs/96f3c1cec4b3ce5d72f708992272e9c1 --file /conda-envs/96f3c1cec4b3ce5d72f708992272e9c1/environment.yaml && \
    conda env create --prefix /conda-envs/8ad6cdcf265d30289788da99d5bf9fff --file /conda-envs/8ad6cdcf265d30289788da99d5bf9fff/environment.yaml && \
    conda env create --prefix /conda-envs/0adafb79cb1bec58ef4c77bf4cca4f95 --file /conda-envs/0adafb79cb1bec58ef4c77bf4cca4f95/environment.yaml && \
    conda env create --prefix /conda-envs/5150d0f0a91d7f7a789a06f453d63479 --file /conda-envs/5150d0f0a91d7f7a789a06f453d63479/environment.yaml && \
    conda clean --all -y

# Step 4: Run post-deploy scripts

RUN conda init && . /root/.bashrc && \
    conda activate /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f && \
    bash /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f/post-deploy.sh && \
    conda deactivate && \
    conda clean --all -y
