
FROM njjai/scrattch_mapping:0.55.4

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    python3-pip \
    python3-venv  # Install python3-venv package

# Set environment variables for Python
ENV PYTHON_VERSION=3.8


# Install Python packages using pip
RUN pip3 install anndata
RUN R -e "install.packages(c('anndata'))"




# Set the working directory
WORKDIR /usr/src/app

# Set entry point
CMD ["R"]
