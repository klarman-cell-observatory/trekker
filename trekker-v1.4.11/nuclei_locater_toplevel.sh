#!/usr/bin/env bash

set -e

#================== Define Top Level Directory Path For Output ===============
OUT_DIR="/mnt/disks/cromwell_root/out"

#========================== Identify Script Directory ==========================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#======================= Function to check samplesheet Format ================
convert_to_linux_format() {
   local file=$1
   if grep -q $'\r$' "$file"; then
      echo "samplesheet contains carriage return characters. Converting to Linux format..."
      sed -i 's/\r$//' "$file"
      echo "Success in removing carriage return characters. samplesheet is updated to Linux format."
   else
      echo "samplesheet is in Linux format."
   fi
}


#==================== Function to Check snRNAseq Inputs ====================
barcodes="barcodes.tsv.gz"; features="features.tsv.gz"; matrix="matrix.mtx.gz"
check_snRNAseq_inputs() {
   if [[ "$SC_PLATFORM" == "TrekkerQ_P" || "$SC_PLATFORM" == "TrekkerQ_P_FXT" ]]; then
      if [[ -f $SC_OUTDIR/$barcodes && -f $SC_OUTDIR/$features && -f $SC_OUTDIR/$matrix ]]; then
         echo "matrix.tsv.gz, barcodes.tsv.gz, and features.tsv.gz found! Starting Trekker Analysis with the existing files using $PROFILE. If you prefer to NOT use these files, please delete them and retrigger the pipeline."
      elif [[ -f "$SC_OUTDIR/count_matrix.mtx" && -f "$SC_OUTDIR/cell_metadata.csv" && -f "$SC_OUTDIR/all_genes.csv" ]] || [[ -f "$SC_OUTDIR/count_matrix.mtx.gz"  && -f "$SC_OUTDIR/cell_metadata.csv.gz"  && -f "$SC_OUTDIR/all_genes.csv.gz" ]]; then
         echo "snRNAseq outputs located."        
      else
	      echo "Abort the Pipeline. Error: Unable to locate single-nuclei RNAseq pipeline outputs. Please ensure cell_metadata.csv, all_genes.csv, count_matrix.mtx (all compressed or all uncompressed) are in the folder specified by 'sc_outdir' in your samplesheet."
         exit 1
      fi
   elif [[ "$SC_PLATFORM" == "TrekkerU_IL" || "$SC_PLATFORM" == "TrekkerU_IL_FXT" ]]; then
      if ls "$SC_OUTDIR/"*scRNA.filtered.barcodes.tsv.gz "$SC_OUTDIR/"*scRNA.filtered.features.tsv.gz "$SC_OUTDIR/"*scRNA.filtered.matrix.mtx.gz 1> /dev/null 2>&1; then
         echo "filtered snRNAseq outputs located."
         mkdir -p "$SC_OUTDIR/trekkerinterim"
         cp "$SC_OUTDIR/"*scRNA.filtered.barcodes.tsv.gz "$SC_OUTDIR/trekkerinterim/barcodes.tsv.gz"
         cp "$SC_OUTDIR/"*scRNA.filtered.features.tsv.gz "$SC_OUTDIR/trekkerinterim/features.tsv.gz"
         cp "$SC_OUTDIR/"*scRNA.filtered.matrix.mtx.gz "$SC_OUTDIR/trekkerinterim/matrix.mtx.gz"
      elif ls "$SC_OUTDIR/"*scRNA.barcodes.tsv.gz "$SC_OUTDIR/"*scRNA.features.tsv.gz "$SC_OUTDIR/"*scRNA.matrix.mtx.gz 1> /dev/null 2>&1; then
         echo "raw snRNAseq outputs located."
         mkdir -p "$SC_OUTDIR/trekkerinterim"
         cp "$SC_OUTDIR/"*scRNA.barcodes.tsv.gz "$SC_OUTDIR/trekkerinterim/barcodes.tsv.gz"
         cp "$SC_OUTDIR/"*scRNA.features.tsv.gz "$SC_OUTDIR/trekkerinterim/features.tsv.gz"
         cp "$SC_OUTDIR/"*scRNA.matrix.mtx.gz "$SC_OUTDIR/trekkerinterim/matrix.mtx.gz"
      else
         echo "Abort the Pipeline. Error: Unable to locate single-nuclei RNAseq pipeline outputs. Please ensure the filtered files: *scRNA.filtered.barcodes.tsv.gz, *scRNA.filtered.features.tsv.gz, and *scRNA.filtered.matrix.mtx.gz or the raw files: scRNA.barcodes.tsv.gz, scRNA.features.tsv.gz and scRNA.matrix.mtx.gz are in the folder specified by 'sc_outdir' in your samplesheet."
         exit 1
      fi
   else
      if [[ -e $SC_OUTDIR/$barcodes && -e $SC_OUTDIR/$features && -e $SC_OUTDIR/$matrix ]]; then
              echo "snRNAseq outputs located."
      else
         echo "Abort the Pipeline. Error: Unable to locate single-nuclei RNAseq pipeline outputs. Please ensure barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz are in the folder specified by 'sc_outdir' in your samplesheet."
         exit 1
      fi
   fi
}


#================ Function to Check nuclei for TrekkerU_RATAC ===============
check_snRNAseq_atac_barcodes() {
   echo "Checking the number of nuclei in snRNAseq and snATACseq count matrices from the single-nuclei pipeline..."
   exp_cells=$(zcat $SC_OUTDIR/$barcodes | wc -l)
   atac_cells=$(zcat ${SCMULTI_OUTDIR}/atac-barcodes.tsv.gz | wc -l)
   if [[ "$exp_cells" -eq "$atac_cells" ]]; then
      echo "Same number of nuclei found."
   else
      echo "Abort the Pipeline. Error: Different number of nuclei found. Please ensure snRNAseq and snATACseq count matrices have the same number of nuclei."
      exit 1
   fi	 
}


#================ Check if samplesheet Argument is Provided =================
if [ "$#" -ne 1 ]; then
   echo "Abort the Pipeline. Error: Please provide samplesheet.csv"
   echo "Usage: bash $0 <samplesheet.csv>"
   exit 1
fi

# Extract samplesheet (first argument)
SAMPLESHEET="$1"


#================ Check if samplesheet Exists, Check Format & Extract Path ================
SAMPLESHEET=$(realpath "$SAMPLESHEET")

echo "Checking samplesheet..."
if [ ! -f "$SAMPLESHEET" ]; then
   echo "Abort the Pipeline. Error: File '$SAMPLESHEET' does not exist. Please provide a valid path to samplesheet".
   exit 1
else
   echo "'$SAMPLESHEET' exists."
   SAMPLESHEET_DIR=$(dirname "$SAMPLESHEET")
   echo "Checking samplesheet format..."
   convert_to_linux_format "$SAMPLESHEET"
fi


#==================== Function to Check if Input File Path Exists ====================
check_path() {
   local path="$1"
   local file_type="$2"
   if [[ ! -f "$path" ]]; then
      if [[ "$file_type" == "vdj_seurat" ]]; then
         echo "Abort the Pipeline. Error: File '$path' does not exist. Please check your samplesheet and provide a valid value to 'sc_sample' (prefix for vdj seurat file) or ensure that file '${SC_SAMPLE}_Seurat.rds' exists in the path specified by 'scmulti_outdir' in your samplesheet."
         exit 1
      elif [[ "$file_type" == "atac-barcodes" || "$file_type" == "atac-features" || "$file_type" == "atac-matrix" ]]; then
         echo "Abort the Pipeline. Error: File '$path' does not exist. Please check your samplesheet and provide a valid path to the 'scmulti_outdir' value or ensure that the '$file_type' can be located in the path specified by 'scmulti_outdir' of your samplesheet."
	 exit 1
      else
         echo "Abort the Pipeline. Error: File '$path' does not exist. Please check your samplesheet and provide a valid path to the '$file_type' value."
	 exit 1
      fi
   else
      echo "'$path' exists."
   fi
}


#==================== Function to Validate FASTQ Names ====================
check_fastq_name() {
   local path="$1"
   local file_type="$2"
   local file_name
   file_name="$(basename "$path")"

   if [[ "$file_type" == "fastq_1" ]]; then
      if [[ "$file_name" != *R1* || "$file_name" != *.fastq.gz ]]; then
         echo "Abort the Pipeline. Error: Incorrect 'fastq_1' value. File_name: '$file_name' must contain 'R1' and end in '.fastq.gz'. Current value is '$path'. Please check your samplesheet and provide a valid 'fastq_1' value."
         exit 1
      fi
   elif [[ "$file_type" == "fastq_2" ]]; then
      if [[ "$file_name" != *R2* || "$file_name" != *.fastq.gz ]]; then
         echo "Abort the Pipeline. Error: Incorrect 'fastq_2' value. File name: '$file_name' must contain 'R2' and end in '.fastq.gz'. Current value is '$path'. Please check your samplesheet and provide a valid 'fastq_2' value."
         exit 1
      fi
   fi
}


#========================= Function to Handle Errors ========================
error_msg() {
   FAILED_COMMAND="$BASH_COMMAND"
   LOG_FILE="${FAILED_COMMAND##*/}"
   STEP_NAME="${LOG_FILE%.log}"
   echo "Abort the Pipeline. Error: ${STEP_NAME} failed. Please refer to '${LOG_DIR}/${LOG_FILE}' for details."
   exit 1
}


#=================== Set the trap to catch any ERR signal ===================
trap 'error_msg' ERR


#======================= Extract samplesheet Columns =======================
SAMPLE_DATA=$(tail -n +2 "$SAMPLESHEET" | head -n 1)
SAMPLE_ID=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $1}')
SC_SAMPLE=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $2}')
ANALYSIS_DATE=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $3}')
TILE_ID_PATH=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $4}')
TILE_ID="${TILE_ID_PATH##*/}" && TILE_ID="${TILE_ID%_*}"
FASTQ_CB=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $5}')
FASTQ_TAGS=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $6}')
SC_OUTDIR=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $7}')
SC_OUTDIR_PATH="$(dirname "$SC_OUTDIR")"
SC_PLATFORM=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $8}')
PROFILE="$(echo "$SAMPLE_DATA" | awk -F ',' '{print $9}' | sed 's/^ *//;s/ *$//')"
PROFILE_UPDATE=$(echo "$PROFILE" | tr '[:upper:]' '[:lower:]')
SUBSAMPLE="$(echo "$SAMPLE_DATA" | awk -F ',' '{print $10}' | sed 's/^ *//;s/ *$//')"
SUBSAMPLE_UPDATE=$(echo "$SUBSAMPLE" | tr '[:upper:]' '[:lower:]')
CORES=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $11}')
LOG_DIR="$(mkdir -p "$SAMPLESHEET_DIR/log/$SAMPLE_ID" && echo "$SAMPLESHEET_DIR/log/$SAMPLE_ID")"

if [[ "$SC_PLATFORM" == "TrekkerU_RVDJ" || "$SC_PLATFORM" == "TrekkerU_RVDJ_FXT" ]]; then
   SCMULTI_OUTDIR=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $12}')
   SCMULTI_PREFIX=${SC_SAMPLE}
elif [[ "$SC_PLATFORM" == "TrekkerU_RATAC" || "$SC_PLATFORM" == "TrekkerU_RATAC_FXT" ]]; then
      SCMULTI_OUTDIR=$(echo "$SAMPLE_DATA" | awk -F ',' '{print $12}')
      SCMULTI_PREFIX="Na"
else
   SCMULTI_OUTDIR="Na"
   SCMULTI_PREFIX="Na"
fi

#=========================== Initial Input Checks ===========================
echo "Checking Trekker FASTQ R1..."
check_fastq_name "$FASTQ_CB" "fastq_1"
check_path "$FASTQ_CB" "fastq_1"
echo "Checking Trekker FASTQ R2..."
check_fastq_name "$FASTQ_TAGS" "fastq_2"
check_path "$FASTQ_TAGS" "fastq_2"
echo "Checking Trekker tile spatial barcode whitelist..."
check_path "$TILE_ID_PATH" "barcode_file"
echo "Checking nuclei subsampling option..."
if [[ "$SUBSAMPLE_UPDATE" == "yes" || "$SUBSAMPLE_UPDATE" == "no" ]]; then
   echo "Nuclei subsampling option is set to: $SUBSAMPLE_UPDATE"
else
   echo "Abort the Pipeline. Error: The 'subsample' value in the samplesheet must be yes or no (case insensitive). Please check your samplesheet and provide a valid value."
   exit 1  
fi
echo "Checking the number of cores requested..."
if [[ "$CORES" =~ ^[0-9]+$ ]]; then
   echo "Number of cores requested is set to: $CORES"
else
   echo "Abort the Pipeline. Error: The 'cores' value in the samplesheet must be an integer. Please check your samplesheet and provide a valid integer value."
   exit 1  
fi
if [[ "$SC_PLATFORM" == "TrekkerU_RVDJ" || "$SC_PLATFORM" == "TrekkerU_RATAC" || "$SC_PLATFORM" == "TrekkerU_RVDJ_FXT" || "$SC_PLATFORM" == "TrekkerU_RATAC_FXT" ]]; then
   if [ -z "$SCMULTI_OUTDIR" ];then
      echo "Checking scmulti_outdir..."
      echo "Abort the Pipeline. Error: The scmulti_outdir value in the samplesheet is empty. Please check your samplesheet and provide a valid value."
      exit 1
   else
      echo "Checking scmulti_outdir..."
      if [ -d "$SCMULTI_OUTDIR" ]; then
         echo "$SCMULTI_OUTDIR exists."
         if [[ "$SC_PLATFORM" == "TrekkerU_RVDJ" || "$SC_PLATFORM" == "TrekkerU_RVDJ_FXT" ]]; then
            echo "Checking scmulti_outdir/vdj_seurat file..."
            check_path "${SCMULTI_OUTDIR}/${SCMULTI_PREFIX}_Seurat.rds" "vdj_seurat"
         else
            echo "Checking scmulti_outdir/atac files..."
            check_path "${SCMULTI_OUTDIR}/atac-barcodes.tsv.gz" "atac-barcodes"
            check_path "${SCMULTI_OUTDIR}/atac-features.tsv.gz" "atac-features"
            check_path "${SCMULTI_OUTDIR}/atac-matrix.mtx.gz" "atac-matrix"
         fi
      else
         echo "Abort the Pipeline. Error: Unable to locate the multiomic output folder '$SCMULTI_OUTDIR' from the single-nuclei pipeline. Please ensure the path specified by 'scmulti_outdir' in your samplesheet exists." 
         exit 1
      fi
   fi
fi


#====================== Run trekker Pipeline ======================
if [[ "$SC_PLATFORM" == "TrekkerC" || "$SC_PLATFORM" == "TrekkerCX" || "$SC_PLATFORM" == "TrekkerU_C" || "$SC_PLATFORM" == "TrekkerU_CX"  || "$SC_PLATFORM" == "TrekkerU_M" || "$SC_PLATFORM" == "TrekkerFX_FLEX" || 
      "$SC_PLATFORM" == "Trekker5C_C" || "$SC_PLATFORM" == "Trekker5C_CX" || "$SC_PLATFORM" == "TrekkerR" || "$SC_PLATFORM" == "TrekkerU_R" || "$SC_PLATFORM" == "TrekkerU_RATAC" || "$SC_PLATFORM" == "TrekkerU_RVDJ" || 
      "$SC_PLATFORM" == "TrekkerQ_P" || "$SC_PLATFORM" == "TrekkerQ_S" || "$SC_PLATFORM" == "TrekkerU_IL" || "$SC_PLATFORM" == "TrekkerU_PIP" || "$SC_PLATFORM" == "TrekkerSHA_WGA" || "$SC_PLATFORM" == "TrekkerC_FXT" || 
      "$SC_PLATFORM" == "TrekkerCX_FXT" || "$SC_PLATFORM" == "TrekkerU_C_FXT" || "$SC_PLATFORM" == "TrekkerU_CX_FXT" || "$SC_PLATFORM" == "TrekkerU_M_FXT" || "$SC_PLATFORM" == "Trekker5C_C_FXT" || "$SC_PLATFORM" == "Trekker5C_CX_FXT" || 
      "$SC_PLATFORM" == "TrekkerR_FXT" || "$SC_PLATFORM" == "TrekkerU_R_FXT" || "$SC_PLATFORM" == "TrekkerU_RATAC_FXT" || "$SC_PLATFORM" == "TrekkerU_RVDJ_FXT" || "$SC_PLATFORM" == "TrekkerQ_P_FXT" || "$SC_PLATFORM" == "TrekkerQ_S_FXT" || 
      "$SC_PLATFORM" == "TrekkerU_IL_FXT" ||  "$SC_PLATFORM" == "TrekkerU_PIP_FXT" || "$SC_PLATFORM" == "TrekkerSHA_WGA_FXT" || "$SC_PLATFORM" == "TrekkerWGA_FFPE" ]]; then
   
   echo "Checking snRNAseq outputs..." 
   if [ -d $SC_OUTDIR ]; then
      check_snRNAseq_inputs
      if [[ "$SC_PLATFORM" == "TrekkerU_RATAC" || "$SC_PLATFORM" == "TrekkerU_RATAC_FXT" ]]; then
         check_snRNAseq_atac_barcodes
      fi
      echo "Starting Trekker Analysis ($SC_PLATFORM) using $PROFILE"
      if [[ "$PROFILE_UPDATE" == "singularity" ]]; then
         source "$SCRIPT_DIR/nuclei_locater_singularity.sh"
      elif [[ "$PROFILE_UPDATE" == "docker" ]]; then
         source "$SCRIPT_DIR/nuclei_locater_docker.sh"
      elif [[ "$PROFILE_UPDATE" == "conda" ]]; then
         source "$SCRIPT_DIR/nuclei_locater_conda.sh"
      else
         echo "Abort the pipeline. Error: The 'profile' value in the samplesheet must be one of the following: singularity, docker, or conda. Please check your samplesheet and provide a valid value."
         exit 1
      fi
   else
      echo "Abort the Pipeline. Error: Unable to locate single-nuclei RNAseq pipeline output folder. Please ensure the path specified by 'sc_outdir' in your samplesheet exists."
      exit 1
   fi
else
   echo "Abort the pipeline. Error: The 'sc_platform' value in the samplesheet must be one of the following: TrekkerC, TrekkerCX, TrekkerU_C, TrekkerU_CX, TrekkerU_M, TrekkerFX_FLEX, Trekker5C_C, Trekker5C_CX, TrekkerR, TrekkerU_R, TrekkerU_RATAC, TrekkerU_RVDJ, TrekkerQ_P, TrekkerQ_S, TrekkerU_IL, TrekkerU_PIP, TrekkerSHA_WGA, TrekkerU_C_FXT, TrekkerC_FXT, TrekkerU_CX_FXT, TrekkerCX_FXT, TrekkerR_FXT, TrekkerU_R_FXT, TrekkerU_RATAC_FXT, TrekkerU_RVDJ_FXT, TrekkerU_M_FXT, TrekkerQ_S_FXT, TrekkerQ_P_FXT, TrekkerU_IL_FXT, TrekkerU_PIP_FXT, TrekkerSHA_WGA_FXT, TrekkerWGA_FFPE, Trekker5C_CX_FXT, or Trekker5C_C_FXT. Please check your samplesheet and provide a valid value."
fi 

