#!env bash

# Install GISAIDR
# https://github.com/Wytamma/GISAIDR
Rscript -e 'devtools::install_github("Wytamma/GISAIDR")'

# Install curl
Rscript -e 'install.packages("curl", repos = "https://cloud.r-project.org/")'
