FROM tbusaspatial/trekker:v1.4.11

WORKDIR /opt

ARG TREKKER_REPO=https://github.com/klarman-cell-observatory/trekker.git
ARG TREKKER_BRANCH=main

# Install system dependencies and Google Cloud CLI
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        curl \
        ca-certificates \
        apt-transport-https \
        gnupg && \
    mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /etc/apt/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        google-cloud-cli && \
    rm -rf /var/lib/apt/lists/*

# Pull the current Trekker repository
RUN git clone \
        --depth 1 \
        --branch "${TREKKER_BRANCH}" \
        "${TREKKER_REPO}" \
        /tmp/trekker-repo

# Install the Trekker launcher
RUN mkdir -p /opt/trekker-v1.4.11/terra && \
    cp \
        /tmp/trekker-repo/trekker-v1.4.11/nuclei_locater_toplevel.sh \
        /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh

# Install CellBender -> 10x conversion utility
RUN cp \
        /tmp/trekker-repo/terra/cellbender_to_10x.py \
        /opt/trekker-v1.4.11/terra/cellbender_to_10x.py

RUN chmod +x \
        /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh \
        /opt/trekker-v1.4.11/terra/cellbender_to_10x.py

# Clean up source repository
RUN rm -rf /tmp/trekker-repo

WORKDIR /opt/trekker-v1.4.11