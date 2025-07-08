FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bd9b69b913e7b4e17d7d5ca7169a5815c3145ea775af87163f18c0b92abc1bf8"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/biopython.yaml
#   prefix: /conda-envs/bd81c49fcb540d7706807c1683ba7200
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python==3.10
#     - biopython==1.81
#     - pandas==2.0.3
#     - pip==23.2.1
#     - mafft==7.525
#     - pip:
#       - gb2seq==0.2.20
RUN mkdir -p /conda-envs/bd81c49fcb540d7706807c1683ba7200
COPY workflow/envs/biopython.yaml /conda-envs/bd81c49fcb540d7706807c1683ba7200/environment.yaml

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
#   prefix: /conda-envs/ee7a2e1b4ec9a7a9999f34dddaea0605
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - freyja==1.4.2
RUN mkdir -p /conda-envs/ee7a2e1b4ec9a7a9999f34dddaea0605
COPY workflow/envs/freyja.yaml /conda-envs/ee7a2e1b4ec9a7a9999f34dddaea0605/environment.yaml

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
#   prefix: /conda-envs/f2a098519cf1f8c4cecb3c13f8c92883
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base==4.3.1
#     - r-gt==0.9.0
#     - quarto==1.3.450
#     - r-jsonlite==1.8.5
#     - r-tidyverse==2.0.0
#     - r-quarto==1.2
#     - r-heatmaply==1.4.2
#     - r-readr==2.1.4
RUN mkdir -p /conda-envs/f2a098519cf1f8c4cecb3c13f8c92883
COPY workflow/envs/quarto_render.yaml /conda-envs/f2a098519cf1f8c4cecb3c13f8c92883/environment.yaml

# Conda environment:
#   source: workflow/envs/renv.yaml
#   prefix: /conda-envs/4b57bfc237ddc217c1f0b04d34dc06ef
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base=4.1.3
#     - r-tidyverse==2.0.0
#     - r-ggrepel==0.9.3
#     - r-stringi==1.7.12
#     - r-ggpubr==0.6.0
#     - bioconductor-ggtree==3.2.0
#     - r-ape==5.7
#     - r-adephylo==1.1_13
#     - r-pegas==1.2
#     - r-data.table==1.14.8
#     - r-future.apply==1.11.0
#     - r-scales==1.2.1
#     - r-showtext==0.9_6
#     - r-jsonlite==1.8.5
#     - r-logger==0.2.2
RUN mkdir -p /conda-envs/4b57bfc237ddc217c1f0b04d34dc06ef
COPY workflow/envs/renv.yaml /conda-envs/4b57bfc237ddc217c1f0b04d34dc06ef/environment.yaml

# Conda environment:
#   source: workflow/envs/snpeff.yaml
#   prefix: /conda-envs/1934df0e4df02a7ee33c52f53f9e3c30
#   channels:
#     - bioconda
#   dependencies:
#     - snpeff==5.1d
RUN mkdir -p /conda-envs/1934df0e4df02a7ee33c52f53f9e3c30
COPY workflow/envs/snpeff.yaml /conda-envs/1934df0e4df02a7ee33c52f53f9e3c30/environment.yaml

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

RUN conda env create --prefix /conda-envs/bd81c49fcb540d7706807c1683ba7200 --file /conda-envs/bd81c49fcb540d7706807c1683ba7200/environment.yaml && \
    conda env create --prefix /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a --file /conda-envs/9439457f932a4fbca3665c9ea1ac2f0a/environment.yaml && \
    conda env create --prefix /conda-envs/ee7a2e1b4ec9a7a9999f34dddaea0605 --file /conda-envs/ee7a2e1b4ec9a7a9999f34dddaea0605/environment.yaml && \
    conda env create --prefix /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f --file /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f/environment.yaml && \
    conda env create --prefix /conda-envs/0a608afb24723cb6fa8aef748f5efbc8 --file /conda-envs/0a608afb24723cb6fa8aef748f5efbc8/environment.yaml && \
    conda env create --prefix /conda-envs/04a3347f94ddf7e21c34bc49e5246076 --file /conda-envs/04a3347f94ddf7e21c34bc49e5246076/environment.yaml && \
    conda env create --prefix /conda-envs/fd645c541ee7a3d43fb9167441b77888 --file /conda-envs/fd645c541ee7a3d43fb9167441b77888/environment.yaml && \
    conda env create --prefix /conda-envs/f2a098519cf1f8c4cecb3c13f8c92883 --file /conda-envs/f2a098519cf1f8c4cecb3c13f8c92883/environment.yaml && \
    conda env create --prefix /conda-envs/4b57bfc237ddc217c1f0b04d34dc06ef --file /conda-envs/4b57bfc237ddc217c1f0b04d34dc06ef/environment.yaml && \
    conda env create --prefix /conda-envs/1934df0e4df02a7ee33c52f53f9e3c30 --file /conda-envs/1934df0e4df02a7ee33c52f53f9e3c30/environment.yaml && \
    conda env create --prefix /conda-envs/5150d0f0a91d7f7a789a06f453d63479 --file /conda-envs/5150d0f0a91d7f7a789a06f453d63479/environment.yaml && \
    conda clean --all -y

# Step 4: Run post-deploy scripts

RUN conda init && . /root/.bashrc && \
    conda activate /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f && \
    bash /conda-envs/3fad3c9cdfa40bee9404f6a2e8fda69f/post-deploy.sh && \
    conda deactivate && \
    conda clean --all -y
