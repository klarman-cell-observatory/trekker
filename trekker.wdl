version 1.0

workflow trekker_sc {
    input {
        File input_samplesheet
        String output_directory
        String memory = "256G"
        Int disk_space = 300
        String docker_registry
    }

    # Read all lines
    Array[String] all_rows = read_lines(input_samplesheet)

    # Remove the header manually by scattering with an index check
    scatter (i in range(length(all_rows))) {
        String line = all_rows[i]

        if (i > 0) {
            call process_sample {
                input:
                    line = line,
                    input_samplesheet = input_samplesheet,
                    output_directory = output_directory,
                    docker_registry = docker_registry,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = disk_space
            }
        }
    }

    output {
        # Array[File?] filtered_outputs = process_sample.result_h5
    }
}

task process_sample {
    input {
        String line
        String output_directory
        String docker_registry
        Int num_cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail
        chmod -R 777 /mnt/disks/cromwell_root/

        FASTQR1_PATH=$(echo "~{line}" | awk -F, '{print $5}')
        FASTQR2_PATH=$(echo "~{line}" | awk -F, '{print $6}')
        BARCODE_PATH=$(echo "~{line}" | awk -F, '{print $4}') #Tile ID
        SAMPLE_NAME=$(echo "~{line}" | awk -F, '{print $1}')
        SC_OUTDIR=$(echo "~{line}" | awk -F, '{print $7}')

        MNT_PATH="/mnt/disks/cromwell_root/"

        echo "Processing sample: $SAMPLE_NAME"
        mkdir -p "$MNT_PATH"/raw
        mkdir -p "$MNT_PATH"/out

        echo "Retrieving file inputs from provided paths in samplesheet"
        gcloud storage cp "~{input_samplesheet}" "$MNT_PATH"/raw/
        gcloud storage cp "$FASTQR1_PATH" "$MNT_PATH"/raw/
        gcloud storage cp "$FASTQR2_PATH" "$MNT_PATH"/raw/
        gcloud storage cp "$BARCODE_PATH" "$MNT_PATH"/raw/
        
        if [["$SC_OUTDIR" == *.h5]]; then
            echo "Cellbender object input detected - pulling and converting Cellbender to Trekker compatible inputs"
            gcloud storage cp "$SC_OUTDIR" "$MNT_PATH"/raw/

            python3 -c '
                
            import pandas as pd
            import scanpy
            import cellbender

            h5_file = "$MNT_PATH"/raw/"$(basename "$SC_OUTDIR")"
            adata = anndata_from_h5(h5_file)
            sc.write_10x_mtx("$MNT_PATH/raw/cellbender_matrix", adata, gex_only=False)

            df = pd.read_csv("$MNT_PATH"/raw/"$(basename "$input_samplesheet")")
            change_list = []
            for i in df['barcode_file','fastq_1', 'fastq_2']:
                i = "$MNT_PATH" + i.split('/')[-1]
                change_list.append(i)

            change_list.append("$MNT_PATH/raw/cellbender_matrix")
            df.loc[df['sample'] == $SAMPLE_NAME, ['barcode_file', 'fastq_1', 'fastq_2', 'sc_outdir']] = change_list
            '

        else
            echo "Running standard analysis"
            gcloud storage cp "$SC_OUTDIR" "$MNT_PATH"/raw/

            python3 -c '
            import pandas as pd

            df = pd.read_csv("$MNT_PATH"/raw/"$(basename "$input_samplesheet")")
            change_list = []
            for i in df['barcode_file','fastq_1', 'fastq_2', 'sc_outdir']:
                i.replace('$SC_OUTDIR', '$MNT_PATH')
                change_list.append(i)

            '
        fi

        bash trekker-v.1.4.11/nuclei_locator_toplevel.sh /mnt/disks/cromwell_root/raw/"$(basename "$input_samplesheet")"
    >>>

    output {
        # File result_h5 = "count_matrix_f100.h5ad"
    }

    runtime {
        docker: docker_registry
        cpu: num_cpu
        memory: memory
        bootDiskSizeGB: 25
        disks: "local-disk ${disk_space} HDD"
        }
    }
