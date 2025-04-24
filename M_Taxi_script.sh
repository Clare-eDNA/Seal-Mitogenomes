#!/bin/bash
#SBATCH -J MTaxi_test
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH -c 8

source ~/miniforge3/etc/profile.d/conda.sh
source ~/tools/MTaxi/venv/bin/activate
export PATH="$HOME/tools/bedtools2/bin:$PATH"
conda activate bcftools-env

# Define paths and reference files
MTAXI_DIR="/home/adacl33p/tools/MTaxi"
REF1="${MTAXI_DIR}/KT693343-NZFS.fasta"
REF2="${MTAXI_DIR}/NC_008418.1-NZSL.fasta"
POLY1="${REF1}.transv_poly.bed"
POLY2="${REF2}.transv_poly.bed"
FILE_DIR="/projects/sciences/maths_stats/wilcox_group/Clare/Seal-mteDNA/together/renamed/trimmed/qualityControl/mapped/Sorted/"

cd "$MTAXI_DIR" || { echo "❌ Cannot access $MTAXI_DIR"; exit 1; }

for i in {1..40}; do
  sample="S${i}"
  echo "🔄 Processing sample: $sample"

  # Write a sample-specific settings.py
  cat <<EOF > settings.py
import pathlib

sample = "${sample}"
FILE_DIR = "${FILE_DIR}"
MTAXI_FOLD = "${MTAXI_DIR}"

bam1 = f"{FILE_DIR}/${sample}.NZFS.sorted.bam"
bam2 = f"{FILE_DIR}/${sample}.NZSL.sorted.bam"
aligned_sample_to_sp1 = bam1
aligned_sample_to_sp2 = bam2

sp1_ref_file = f"${REF1}"
sp2_ref_file = f"${REF2}"
transv_poly_sp1_file = f"${REF1}.transv_poly.bed"
transv_poly_sp2_file = f"${REF2}.transv_poly.bed"

output = f"{sample}_MTaxi_results.tsv"
out_file = f"{sample}_MTaxi_res.tsv"
debug = False
EOF

  echo "✅ Created settings.py for $sample"
    echo "▶️ Running MTaxi for $sample"
  # Run MTaxi for this sample
  python ~/tools/MTaxi/run_MTaxi.py compare \
    --bam1 "${FILE_DIR}/${sample}.FS.sorted.bam" \
    --bam2 "${FILE_DIR}/${sample}.SL.sorted.bam" \
    --ref1 "${MTAXI_DIR}/${REF1}" \
    --ref2 "${MTAXI_DIR}/${REF2}" \
    --poly1 "${MTAXI_DIR}/${REF1}.transv_poly.bed" \
    --poly2 "${MTAXI_DIR}/${REF2}.transv_poly.bed" \
    --sample "${sample}" \
    --output "${sample}_MTaxi_results.tsv" || echo "❌ CLI mode failed for $sample"

  echo "✅ Finished $sample"
  echo "--------------------------------------"
done

for file in *_MTaxi_results.tsv; do
    echo "📄 File: $file"
    cat "$file"
    echo "--------------------------------------"
done
