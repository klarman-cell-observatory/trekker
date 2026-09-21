FROM tbusaspatial/trekker:v1.4.11

WORKDIR /opt

# ----------------------------------------------------------------------
# Add the Trekker v1.4.11 pipeline files.
#
# This is the local copy of the Takara Trekker distribution, including
# the Terra-specific edit to nuclei_locater_toplevel.sh:
#
#     OUT_DIR="/mnt/disks/cromwell_root/out"
#
# SCRIPT_DIR remains dynamically determined by the Trekker script.
# ----------------------------------------------------------------------

COPY trekker-v1.4.11 /opt/trekker-v1.4.11

# ----------------------------------------------------------------------
# Add the CellBender -> 10x conversion utility.
# ----------------------------------------------------------------------

COPY terra/cellbender_to_10x.py \
    /opt/trekker-v1.4.11/terra/cellbender_to_10x.py

# ----------------------------------------------------------------------
# Make pipeline scripts executable.
# ----------------------------------------------------------------------

RUN chmod +x \
    /opt/trekker-v1.4.11/*.sh \
    /opt/trekker-v1.4.11/terra/cellbender_to_10x.py

# ----------------------------------------------------------------------
# Use Trekker's environment as the default environment.
# ----------------------------------------------------------------------

WORKDIR /opt/trekker-v1.4.11