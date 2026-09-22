#!/usr/bin/env bash


#========================== Define Conda Paths =======================
PROFILE_CONDA_PATH=/opt/conda/envs/trekker
source /opt/conda/etc/profile.d/conda.sh

conda activate ${PROFILE_CONDA_PATH}


TREKKEROUT_ROOT=${OUT_DIR}/${ANALYSIS_DATE}_${SAMPLE_ID}/
TREKKEROUT_SAMP=${TREKKEROUT_ROOT}/trekker_${SAMPLE_ID}/
TREKKEROUT_MAIN=${TREKKEROUT_SAMP}/output/
TREKKEROUT_MISC=${TREKKEROUT_SAMP}/misc/
TREKKEROUT_TILE=${TREKKEROUT_MISC}/${TILE_ID}/


mkdir -p ${TREKKEROUT_ROOT}
mkdir -p ${TREKKEROUT_SAMP}
mkdir -p ${TREKKEROUT_MAIN}
mkdir -p ${TREKKEROUT_MISC}
mkdir -p ${TREKKEROUT_TILE}


#==================== Function to Run TrekkerU PIP Converter ====================
run_trekkerU_PIP_converter() {
   fastq_r1="$1"
   sc_outdir="$2"

   fastq_dir=$(dirname "$fastq_r1")
   filename=$(basename "$fastq_r1")
   fastq_prefix="${filename%_R1*}"

   if [[ "$filename" =~ _converted_R1\.fastq\.gz$ ]] && [[ "$sc_outdir" =~ /converted/?$ ]]; then
      echo "fastq_R1 and snRNAseq outputs are both properly reformatted already. Skipping TrekkerU_PIP preprocessing."
   elif [[ "$filename" =~ _converted_R1\.fastq\.gz$ ]] && [[ ! "$sc_outdir" =~ /converted/?$ ]]; then
      echo "Abort the Pipeline. Error: Wrong fastq_1 format. You are using a Trekker R1 reformatted by PIPSeeker. Please modify the 'fastq_1' value in your samplesheet to point to the raw Trekker R1 without reformatting."
      exit 1
   elif [[ ! "$filename" =~ _converted_R1\.fastq\.gz$ ]] && [[ "$sc_outdir" =~ /converted/?$ ]]; then
      echo "Abort the Pipeline. Error: Wrong file format in sc_outdir. You are using a snRNAseq output folder reformatted by trekkerU_PIP_converter. Please modify the 'sc_outdir' value in your samplesheet to point to the raw snRNAseq outputs without reformatting."
      exit 1
   else
      echo "Start TrekkerU_PIP preprocessing..."
      bash "$SCRIPT_DIR/common/TrekkerU_PIP/trekkerU_PIP_converter.sh" \
         "$SCRIPT_DIR" "$fastq_r1" "$fastq_prefix" "$sc_outdir" "V" "$LOG_DIR"
      echo "TrekkerU_PIP preprocessing Done"
      FASTQ_CB="${fastq_dir}/${fastq_prefix}/barcoded_fastqs/${fastq_prefix}_converted_R1.fastq.gz"
      SC_OUTDIR="${sc_outdir}/converted"
   fi
}


if [[ "$SC_PLATFORM" == "TrekkerU_PIP" || "$SC_PLATFORM" == "TrekkerU_PIP_FXT" ]]; then
   run_trekkerU_PIP_converter "$FASTQ_CB" "$SC_OUTDIR"
fi


if [[ "$SC_PLATFORM" == "TrekkerFX_FLEX" ]]; then
   echo "Start check_mtx_status"
   python ${SCRIPT_DIR}/common/check_mtx_status.py \
      ${SC_OUTDIR}/barcodes.tsv.gz >& \
      ${LOG_DIR}/check_mtx_status.log
   echo "check_mtx_status Done"
fi


if [[ "$SC_PLATFORM" == "TrekkerQ_P" || "$SC_PLATFORM" == "TrekkerQ_P_FXT" ]]; then
   if [[ ! -f $SC_OUTDIR/$barcodes || ! -f $SC_OUTDIR/$features || ! -f $SC_OUTDIR/$matrix ]]; then
      if [[ -f "$SC_OUTDIR/count_matrix.mtx" && -f "$SC_OUTDIR/cell_metadata.csv" && -f "$SC_OUTDIR/all_genes.csv" ]] || [[ -f "$SC_OUTDIR/count_matrix.mtx.gz" && -f "$SC_OUTDIR/cell_metadata.csv.gz" && -f "$SC_OUTDIR/all_genes.csv.gz" ]]; then
         echo "Start parse_preprocessing"
         python ${SCRIPT_DIR}/common/TrekkerQ_P/parse_preprocessing.py \
            ${SC_OUTDIR} \
            ${SCRIPT_DIR}/common/TrekkerQ_P >& \
            ${LOG_DIR}/parse_preprocessing.log
         echo "parse_preprocessing Done"
      fi
   fi
fi


echo "Start fastq_parser"
python ${SCRIPT_DIR}/common/fastq_parser.py \
   ${SAMPLE_ID} \
   ${FASTQ_CB} \
   ${FASTQ_TAGS} \
   ${SC_OUTDIR} \
   ${SCRIPT_DIR}/cellbarcode_whitelists/ \
   ${SC_PLATFORM} \
   ${TREKKEROUT_MISC} >& \
   ${LOG_DIR}/fastq_parser.log
echo "fastq_parser Done"


echo "Start splitspatialbarcodes"
python ${SCRIPT_DIR}/common/splitspatialbarcodes.py \
   ${TILE_ID_PATH} \
   ${TREKKEROUT_TILE} >& \
   ${LOG_DIR}/beadbarcode_splitter.log
echo "splitspatialbarcodes Done"


echo "Start bead_matching"
python ${SCRIPT_DIR}/common/bead_matching.py \
   -b ${TREKKEROUT_MISC}/reads_per_SB_${SAMPLE_ID}.txt \
   -i ${TILE_ID} \
   -d ${TREKKEROUT_MISC} \
   -o ${TREKKEROUT_MISC}/matching_result_${SAMPLE_ID}.csv >& \
   ${LOG_DIR}/bead_matching.log
echo "bead_matching Done"


echo "Start spatial"
Rscript --vanilla ${SCRIPT_DIR}/common/spatial.R \
   ${SAMPLE_ID} \
   ${TREKKEROUT_MISC} \
   ${TREKKEROUT_MAIN} \
   ${SC_PLATFORM} \
   ${SUBSAMPLE_UPDATE} \
   ${CORES} >& \
   ${LOG_DIR}/spatial.log
echo "spatial Done"


echo "Start analysis"
Rscript --vanilla ${SCRIPT_DIR}/common/analysis.R \
   ${SAMPLE_ID} \
   ${SC_PLATFORM} \
   ${SCRIPT_DIR} \
   ${SC_OUTDIR} \
   ${TREKKEROUT_MISC} \
   ${TREKKEROUT_MAIN} \
   ${PROFILE_CONDA_PATH} \
   ${SCMULTI_OUTDIR} \
   ${SCMULTI_PREFIX} >& \
   ${LOG_DIR}/analysis.log
echo "analysis Done"


echo "Start mismatch_analysis"
python ${SCRIPT_DIR}/common/extended_qc/mismatch_analysis.py \
   ${SAMPLE_ID} \
   ${TREKKEROUT_MISC} >& \
   ${LOG_DIR}/mismatch_analysis.log
echo "mismatch_analysis Done"


echo "Start genmetrics"
python ${SCRIPT_DIR}/common/genmetrics.py \
   ${SAMPLE_ID} \
   ${TILE_ID} \
   ${SC_PLATFORM} \
   ${LOG_DIR} \
   ${TREKKEROUT_MISC} \
   ${TREKKEROUT_MAIN} >& \
   ${LOG_DIR}/genmetrics.log
echo "genmetrics Done"


echo "Start genreport"
Rscript --vanilla ${SCRIPT_DIR}/common/genreport.R \
   ${SAMPLE_ID} \
   ${TILE_ID} \
   ${SCRIPT_DIR}/common/htmlrender.Rmd \
   ${TREKKEROUT_MAIN} \
   ${SC_PLATFORM} \
   ${TREKKEROUT_MAIN}/${SAMPLE_ID}_summary_metrics.csv \
   ${TREKKEROUT_MAIN}/${SAMPLE_ID}_variable_features_clusters.csv \
   ${TREKKEROUT_MAIN}/${SAMPLE_ID}_variable_features_spatial_moransi.txt \
   ${TREKKEROUT_MAIN}/intermediates/${SAMPLE_ID}_seurat_spatial.rds \
   ${TREKKEROUT_MAIN}/${SAMPLE_ID}_ConfPositioned_seurat_spatial.rds \
   ${TREKKEROUT_MISC}/${SAMPLE_ID}_mismatchfreq.csv >& \
   ${LOG_DIR}/genreport.log
echo "genreport Done"
