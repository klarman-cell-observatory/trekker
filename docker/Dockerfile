FROM tbusaspatial/trekker:v1.4.11

WORKDIR /opt

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

ARG TREKKER_REPO=https://github.com/chenebroad/trekker.git
ARG TREKKER_BRANCH=main

# ----------------------------------------------------------------------
# Install git if it is not already available in the Trekker image.
# ----------------------------------------------------------------------

RUN if ! command -v git >/dev/null 2>&1; then \
        apt-get update && \
        apt-get install -y --no-install-recommends git && \
        rm -rf /var/lib/apt/lists/*; \
    fi

# ----------------------------------------------------------------------
# Clone the current version of the Terra-specific Trekker code.
#
# The repository contains:
#
#   trekker-v1.4.11/
#       nuclei_locater_toplevel.sh
#
#   terra/
#       cellbender_to_10x.py
#
# A fresh clone is performed every Docker build.
# ----------------------------------------------------------------------

RUN git clone \
        --depth 1 \
        --branch "${TREKKER_BRANCH}" \
        "${TREKKER_REPO}" \
        /tmp/trekker-repo

# ----------------------------------------------------------------------
# Install the modified Trekker top-level script.
#
# The Trekker base image supplies the Trekker environment/dependencies.
# The GitHub repository supplies our Terra-specific version of the
# pipeline entrypoint.
# ----------------------------------------------------------------------

RUN mkdir -p /opt/trekker-v1.4.11/terra && \
    cp \
        /tmp/trekker-repo/trekker-v1.4.11/nuclei_locater_toplevel.sh \
        /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh

# ----------------------------------------------------------------------
# Install the CellBender -> 10x conversion utility.
# ----------------------------------------------------------------------

RUN cp \
        /tmp/trekker-repo/terra/cellbender_to_10x.py \
        /opt/trekker-v1.4.11/terra/cellbender_to_10x.py

# ----------------------------------------------------------------------
# Make the scripts executable.
# ----------------------------------------------------------------------

RUN chmod +x \
        /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh \
        /opt/trekker-v1.4.11/terra/cellbender_to_10x.py

# ----------------------------------------------------------------------
# Remove the temporary GitHub checkout.
# ----------------------------------------------------------------------

RUN rm -rf /tmp/trekker-repo

# ----------------------------------------------------------------------
# Trekker working directory.
# ----------------------------------------------------------------------

WORKDIR /opt/trekker-v1.4.11