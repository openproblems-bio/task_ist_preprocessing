# Use a bioconductor base image
FROM bioconductor/bioconductor:3.17

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    python3-pip \
    python3-venv  # Install python3-venv package

# Set environment variables for Python
ENV PYTHON_VERSION=3.8

# Install R packages
RUN R -e "install.packages(c('data.table','Seurat', 'yaml', 'optparse', 'devtools','purrr'))"
RUN R -e "install.packages('ggdendro', repos='https://cran.rstudio.com')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('GO.db', 'WGCNA'))"

# Install Python packages using pip
RUN pip3 install anndata

RUN R -e "install.packages(c('anndata'))"
# Install other R packages from GitHub
RUN R -e "devtools::install_github('AllenInstitute/scrattch.vis')"
RUN R -e "devtools::install_github('AllenInstitute/mfishtools')"

# Set the working directory
WORKDIR /usr/src/app

# Set entry point
CMD ["R"]
