version 1.0

workflow trekker_sc {

    input {
        File input_samplesheet
        String output_directory
        String docker_registry

        Int num_cpu = 16
        String memory = "256G"
        Int disk_space = 300
    }

    Array[String] all_rows = read_lines(input_samplesheet)

    scatter (row_index in range(length(all_rows) - 1)) {

        Int row_index_plus_one = row_index + 1

        String sample_line = all_rows[row_index_plus_one]

        call process_sample {
            input:
                line = sample_line,
                input_samplesheet = input_samplesheet,
                output_directory = output_directory,
                docker_registry = docker_registry,
                num_cpu = num_cpu,
                memory = memory,
                disk_space = disk_space
        }
    }

    output {
        Array[File] trekker_outputs = process_sample.trekker_output_files
    }
}


task process_sample {

    input {
        String line
        File input_samplesheet
        String output_directory
        String docker_registry

        Int num_cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail

        MNT_PATH="/mnt/disks/cromwell_root"
        RAW_DIR="${MNT_PATH}/raw"
        OUT_DIR="${MNT_PATH}/out"

        mkdir -p "${RAW_DIR}"
        mkdir -p "${OUT_DIR}"

        chmod -R 777 "${MNT_PATH}"

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
        echo "Output directory:  ~{output_directory}"
        echo ""

        echo "============================================================"
        echo "Container diagnostics"
        echo "============================================================"

        echo "PATH=${PATH}"

        echo "gcloud location:"
        which gcloud || true

        echo "gcloud version:"
        gcloud --version

        echo ""

        LOCAL_SAMPLE_SHEET="${RAW_DIR}/$(basename "~{input_samplesheet}")"
        cp "~{input_samplesheet}" "${LOCAL_SAMPLE_SHEET}"

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

with open(input_csv, "r", newline="") as infile:
    reader = csv.reader(infile)
    rows = list(reader)

if not rows:
    raise RuntimeError("Samplesheet is empty.")

header = rows[0]

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

        echo "Trekker output directory:"
        echo "  ${OUT_DIR}"
        echo ""

        echo "Trekker output files:"

        find \
            "${OUT_DIR}" \
            -type f \
            -print

        echo ""

        OUTPUT_FILE_COUNT=$(find \
            "${OUT_DIR}" \
            -type f \
            | wc -l)

        echo "Trekker output file count: ${OUTPUT_FILE_COUNT}"

        if [[ "${OUTPUT_FILE_COUNT}" -eq 0 ]]; then
            echo ""
            echo "ERROR: Trekker completed but produced no output files."
            echo ""
            echo "Contents of output directory:"
            find \
                "${OUT_DIR}" \
                -print
            exit 1
        fi

        echo ""
        echo "Trekker produced ${OUTPUT_FILE_COUNT} output file(s)."
        echo ""

        DESTINATION="~{output_directory}"
        DESTINATION="${DESTINATION%/}"

        SAMPLE_DESTINATION="${DESTINATION}/${SAMPLE_NAME}"

        echo "============================================================"
        echo "Uploading Trekker outputs"
        echo "============================================================"

        echo "Destination:"
        echo "  ${SAMPLE_DESTINATION}"
        echo ""

        while IFS= read -r FILE; do

            RELATIVE_PATH="${FILE#${OUT_DIR}/}"

            DESTINATION_FILE="${SAMPLE_DESTINATION}/${RELATIVE_PATH}"

            echo "Uploading:"
            echo "  ${FILE}"
            echo "To:"
            echo "  ${DESTINATION_FILE}"
            echo ""

            gcloud storage cp \
                "${FILE}" \
                "${DESTINATION_FILE}"

        done < <(
            find \
                "${OUT_DIR}" \
                -type f \
                -print
        )

        echo ""
        echo "============================================================"
        echo "Trekker output upload complete"
        echo "============================================================"

    >>>

    output {
        Array[File] trekker_outputs = flatten(process_sample.trekker_output_files)
    }

    runtime {
        docker: docker_registry
        cpu: num_cpu
        memory: memory
        bootDiskSizeGB: 25
        disks: "local-disk ${disk_space} HDD"
    }
}