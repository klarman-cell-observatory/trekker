FROM continuumio/miniconda3:latest

# Set shell to bash for predictable script sourcing
SHELL ["/bin/bash", "-c"]

# Set the working directory inside the container
WORKDIR /mnt/

# 1. Installation: Download and extract the resource archive
RUN apt-get update && apt-get install -y wget gcc g++ && \
    wget https://takarabio.com -O - | tar -xzf -

# 2. Change the output directory programmatically
RUN mkdir -p /home/jupyter/trekker_out && \
    cd trekker-* && \
    sed -i 's|#OUTPUT_DIR=.*|OUTPUT_DIR=/home/jupyter/trekker_out|g' nuclei_locater_toplevel.sh

# 3. Create the environment & Install baseline workflow dependencies
RUN cd trekker-* && \
    conda env create -f environment.yml && \
    conda run -n trekker conda install -y cytoolz toolz && \
    conda run -n trekker conda install -y -c conda-forge importlib_metadata

# 4. Install Pandas, Scanpy, and Cellbender into the 'trekker' environment
# Using standard bioconda channels for scanpy/pandas, and pip for cellbender to safely track PyTorch
RUN conda run -n trekker conda install -y -c conda-forge -c bioconda \
        pandas \
        scanpy \
        python-igraph \
        leidenalg && \
    conda run -n trekker pip install --no-cache-dir cellbender

# 5. Change the environment path in the shell script
RUN cd trekker-* && \
    sed -i 's|PROFILE_CONDA_PATH=.*|PROFILE_CONDA_PATH=/opt/conda/envs/trekker|g' nuclei_locater_conda.sh

# Set the trekker environment to be permanently active when entering the container
RUN echo "conda activate trekker" >> ~/.bashrc

# Default command when the container starts
CMD ["/bin/bash"]