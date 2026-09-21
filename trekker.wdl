version 1.0

workflow trekker_sc {

    input {

        # Original Trekker samplesheet stored in GCS.
        File input_samplesheet

        # Retained for compatibility with the existing workflow configuration.
        String output_directory

        # Full Docker image URI.
        String docker_registry

        Int num_cpu = 16
        String memory = "256G"
        Int disk_space = 300
    }

    # Read the complete samplesheet.
    Array[String] all_rows = read_lines(input_samplesheet)

    # Remove the first row, which is the samplesheet header.
    #
    # sub(array, start, length)
    #
    # Starting at index 1 means the header at index 0 is excluded.
    Array[String] sample_rows = sub(
        all_rows,
        1,
        length(all_rows) - 1
    )

    # Process one sample per task.
    scatter (sample_line in sample_rows) {

        call process_sample {
            input:
                line = sample_line,
                input_samplesheet = input_samplesheet,
                docker_registry = docker_registry,
                num_cpu = num_cpu,
                memory = memory,
                disk_space = disk_space
        }
    }

    output {

        Array[Directory] trekker_outputs = process_sample.trekker_output
    }
}


task process_sample {

    input {

        # One data row from the Trekker samplesheet.
        String line

        # Original samplesheet.
        File input_samplesheet

        # Docker image containing Trekker and the CellBender converter.
        String docker_registry

        Int num_cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail

        # ------------------------------------------------------------------
        # Terra working directories
        # ------------------------------------------------------------------

        MNT_PATH="/mnt/disks/cromwell_root"
        RAW_DIR="${MNT_PATH}/raw"
        OUT_DIR="${MNT_PATH}/out"

        mkdir -p "${RAW_DIR}"
        mkdir -p "${OUT_DIR}"

        chmod -R 777 "${MNT_PATH}"

        # ------------------------------------------------------------------
        # Parse samplesheet row
        #
        # Expected columns:
        #
        #   1 = sample
        #   2 = ...
        #   3 = ...
        #   4 = barcode_file
        #   5 = fastq_1
        #   6 = fastq_2
        #   7 = sc_outdir
        #
        # ------------------------------------------------------------------

        SAMPLE_NAME=$(echo "~{line}" | awk -F',' '{print $1}')
        BARCODE_PATH=$(echo "~{line}" | awk -F',' '{print $4}')
        FASTQR1_PATH=$(echo "~{line}" | awk -F',' '{print $5}')
        FASTQR2_PATH=$(echo "~{line}" | awk -F',' '{print $6}')
        SC_OUTDIR=$(echo "~{line}" | awk -F',' '{print $7}')

        echo "============================================================"
        echo "Processing sample: ${SAMPLE_NAME}"
        echo "============================================================"

        echo "Barcode/tile path: ${BARCODE_PATH}"
        echo "FASTQ R1:          ${FASTQR1_PATH}"
        echo "FASTQ R2:          ${FASTQR2_PATH}"
        echo "sc_outdir:         ${SC_OUTDIR}"
        echo ""

        # ------------------------------------------------------------------
        # Localize the original samplesheet.
        #
        # Cromwell/Terra already localizes WDL File inputs, so no GCS
        # download is needed here.
        # ------------------------------------------------------------------

        LOCAL_SAMPLE_SHEET="${RAW_DIR}/$(basename "~{input_samplesheet}")"

        cp "~{input_samplesheet}" "${LOCAL_SAMPLE_SHEET}"

        # ------------------------------------------------------------------
        # Download FASTQs and barcode/tile input.
        #
        # These paths are stored as GCS URIs in the samplesheet.
        # ------------------------------------------------------------------

        echo "Downloading FASTQ R1..."

        gcloud storage cp \
            "${FASTQR1_PATH}" \
            "${RAW_DIR}/"

        echo "Downloading FASTQ R2..."

        gcloud storage cp \
            "${FASTQR2_PATH}" \
            "${RAW_DIR}/"

        echo "Downloading barcode/tile file..."

        gcloud storage cp \
            "${BARCODE_PATH}" \
            "${RAW_DIR}/"

        LOCAL_R1="${RAW_DIR}/$(basename "${FASTQR1_PATH}")"
        LOCAL_R2="${RAW_DIR}/$(basename "${FASTQR2_PATH}")"
        LOCAL_BARCODE="${RAW_DIR}/$(basename "${BARCODE_PATH}")"

        echo ""
        echo "Localized files:"
        echo "  R1:      ${LOCAL_R1}"
        echo "  R2:      ${LOCAL_R2}"
        echo "  Barcode: ${LOCAL_BARCODE}"
        echo ""

        # ------------------------------------------------------------------
        # Handle CellBender / single-cell input.
        #
        # If sc_outdir is an H5 file, treat it as CellBender output and
        # convert it into a Trekker-compatible 10x Matrix Market directory.
        #
        # Otherwise, download the supplied sc_outdir directly.
        # ------------------------------------------------------------------

        if [[ "${SC_OUTDIR}" == *.h5 ]]; then

            echo "CellBender H5 input detected."

            echo "Downloading CellBender output..."

            gcloud storage cp \
                "${SC_OUTDIR}" \
                "${RAW_DIR}/"

            LOCAL_CELLBENDER_H5="${RAW_DIR}/$(basename "${SC_OUTDIR}")"

            CELLBENDER_MATRIX_DIR="${RAW_DIR}/cellbender_matrix"

            echo ""
            echo "CellBender H5:"
            echo "  ${LOCAL_CELLBENDER_H5}"
            echo ""

            echo "Converting CellBender H5 to 10x Matrix Market..."

            python3 \
                /opt/trekker-v1.4.11/terra/cellbender_to_10x.py \
                "${LOCAL_CELLBENDER_H5}" \
                "${CELLBENDER_MATRIX_DIR}"

            echo ""
            echo "CellBender conversion complete."

            echo "Converted matrix contents:"

            ls -lh "${CELLBENDER_MATRIX_DIR}"

            LOCAL_SC_OUTDIR="${CELLBENDER_MATRIX_DIR}"

        else

            echo "Non-H5 sc_outdir detected."

            echo "Downloading supplied sc_outdir..."

            gcloud storage cp \
                -r \
                "${SC_OUTDIR}" \
                "${RAW_DIR}/"

            LOCAL_SC_OUTDIR="${RAW_DIR}/$(basename "${SC_OUTDIR}")"

        fi

        # ------------------------------------------------------------------
        # Create a local Trekker samplesheet.
        #
        # Trekker receives local filesystem paths instead of GCS URIs.
        #
        # We preserve the original samplesheet and replace the paths for
        # the sample being processed.
        # ------------------------------------------------------------------

        LOCAL_TREKKER_SAMPLESHEET="${RAW_DIR}/trekker_samplesheet.csv"

        python3 - \
            "${LOCAL_SAMPLE_SHEET}" \
            "${LOCAL_TREKKER_SAMPLESHEET}" \
            "${SAMPLE_NAME}" \
            "${LOCAL_BARCODE}" \
            "${LOCAL_R1}" \
            "${LOCAL_R2}" \
            "${LOCAL_SC_OUTDIR}" \
            <<'PY'
import csv
import sys


(
    input_csv,
    output_csv,
    sample_name,
    barcode_path,
    fastq1_path,
    fastq2_path,
    sc_outdir,
) = sys.argv[1:]


# ----------------------------------------------------------------------
# Read samplesheet
# ----------------------------------------------------------------------

with open(input_csv, "r", newline="") as infile:
    reader = csv.reader(infile)
    rows = list(reader)


if not rows:
    raise RuntimeError("Samplesheet is empty.")


header = rows[0]


# ----------------------------------------------------------------------
# Find required columns
# ----------------------------------------------------------------------

header_lookup = {
    name.strip(): index
    for index, name in enumerate(header)
}


required_columns = [
    "sample",
    "barcode_file",
    "fastq_1",
    "fastq_2",
    "sc_outdir",
]


missing = [
    column
    for column in required_columns
    if column not in header_lookup
]


if missing:
    raise RuntimeError(
        "Samplesheet is missing required columns: "
        + ", ".join(missing)
    )


sample_idx = header_lookup["sample"]
barcode_idx = header_lookup["barcode_file"]
fastq1_idx = header_lookup["fastq_1"]
fastq2_idx = header_lookup["fastq_2"]
sc_outdir_idx = header_lookup["sc_outdir"]


# ----------------------------------------------------------------------
# Replace paths for the requested sample
# ----------------------------------------------------------------------

updated_rows = [header]

found_sample = False


for row in rows[1:]:

    if not row:
        continue

    if len(row) < len(header):
        raise RuntimeError(
            f"Malformed samplesheet row with {len(row)} columns; "
            f"expected at least {len(header)}."
        )

    if row[sample_idx] == sample_name:

        row[barcode_idx] = barcode_path
        row[fastq1_idx] = fastq1_path
        row[fastq2_idx] = fastq2_path
        row[sc_outdir_idx] = sc_outdir

        found_sample = True

    updated_rows.append(row)


if not found_sample:
    raise RuntimeError(
        f"Could not find sample '{sample_name}' in samplesheet."
    )


# ----------------------------------------------------------------------
# Write localized samplesheet
# ----------------------------------------------------------------------

with open(output_csv, "w", newline="") as outfile:

    writer = csv.writer(outfile)

    writer.writerows(updated_rows)


print(
    f"Wrote localized samplesheet: {output_csv}"
)
PY

        echo ""
        echo "============================================================"
        echo "Localized Trekker samplesheet"
        echo "============================================================"

        cat "${LOCAL_TREKKER_SAMPLESHEET}"

        echo ""

        # ------------------------------------------------------------------
        # Run Trekker
        #
        # The Trekker script is installed by the Docker image at:
        #
        #   /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh
        #
        # The modified script should contain:
        #
        #   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        #
        # and:
        #
        #   OUT_DIR="/mnt/disks/cromwell_root/out"
        # ------------------------------------------------------------------

        echo "============================================================"
        echo "Starting Trekker"
        echo "============================================================"

        bash \
            /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh \
            "${LOCAL_TREKKER_SAMPLESHEET}"

        echo ""
        echo "============================================================"
        echo "Trekker complete"
        echo "============================================================"

        echo "Trekker output:"

        find \
            "${OUT_DIR}" \
            -maxdepth 3 \
            -type f \
            -print

    >>>

    output {

        Directory trekker_output = "/mnt/disks/cromwell_root/out"
    }

    runtime {

        docker: docker_registry

        cpu: num_cpu

        memory: memory

        bootDiskSizeGB: 25

        disks: "local-disk ${disk_space} HDD"
    }
}