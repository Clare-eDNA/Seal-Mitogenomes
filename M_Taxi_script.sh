#!/bin/bash
## note that this is a series of scripts, the first two prep the data from the trimmed and quality controlled data, to sorting/mapping the NZ Fur Seal (target species) and NZ Sea Lion (off target pinniped who overlaps in the Moeraki region)
## Step 4.1: Mapping to Arctocephalus forsteri mitogenome

#!/bin/bash
#SBATCH -J mapping
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH -c 8
module load samtools/1.19.2-eymmh4o
module load bwa/0.7.17-tkijtxm

cd /projects/sciences/maths_stats/wilcox_group/Clare/Seal-mteDNA/together/renamed/trimmed/qualityControl/
ref="/projects/sciences/maths_stats/wilcox_group/Clare/Seal-mteDNA/together/renamed/trimmed/qualityControl/RefmtGenom"
outdir="/projects/sciences/maths_stats/wilcox_group/Clare/Seal-mteDNA/together/renamed/trimmed/qualityControl/mapped"


# Only index if necessary
[ -f "${ref}/KT693343-NZFS.fasta.bwt" ] || bwa index "${ref}/KT693343-NZFS.fasta"
[ -f "${ref}/NC_008418.1-NZSL.fasta.bwt" ] || bwa index "${ref}/NC_008418.1-NZSL.fasta"

for sample in S*1P.fq.gz; do
    base="${sample%_1P.fq.gz}"
    echo "Mapping: ${base}"

    bwa mem -t 8 "${ref}/NC_008418.1-NZSL.fasta" "${base}_1P.fq.gz" "${base}_2P.fq.gz" > "${outdir}/${base}.SL.pe.sam"
    bwa mem -t 8 "${ref}/KT693343-NZFS.fasta" "${base}_1P.fq.gz" "${base}_2P.fq.gz" > "${outdir}/${base}.FS.pe.sam"
done


## Step 4.2: Mapping to Arctocephalus forsteri and sea lion (P hookeri) mitogenome

#!/bin/bash
#SBATCH -J SortingMappingNZ&SL
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH -c 8
set -euo pipefail
module load samtools/1.19.2-eymmh4o
module load bwa/0.7.17-tkijtxm

cd /projects/sciences/maths_stats/wilcox_group/Clare/Seal-mteDNA/together/renamed/trimmed/qualityControl/mapped/

log_error() {
    echo "❌ ERROR at step: $1 for sample: $2" >&2
}

# do sea lions
for samfile in *.SL.pe.sam; do
    base="${samfile%.SL.pe.sam}"
    echo "🔄 Processing: $base"
    if [[ ! -s "${base}.SL.pe.sam" ]]; then
    echo "⚠️ ${base}.SL.pe.sam is empty — skipping."
    continue
fi
    samtools view -b -@ 8 "${base}.SL.pe.sam" | samtools sort -@ 8 -o "${base}.SL.sorted.bam"
    samtools index "${base}.SL.sorted.bam"
    samtools collate -o "${base}.SL.nameCol.bam" "${base}.SL.sorted.bam"
    samtools fixmate -m "${base}.SL.nameCol.bam" "${base}.SL.fixMate.bam"
    samtools sort -o "${base}.SL.pos.bam" "${base}.SL.fixMate.bam"
    samtools markdup -r -s "${base}.SL.pos.bam" "${base}.SL.markDup.bam"
    echo "✅ Done with SL $base"
done

# do fur seals
for samfile in *.FS.pe.sam; do
    base="${samfile%.FS.pe.sam}"
    echo "🔄 Processing: $base"
    if [[ ! -s "${base}.FS.pe.sam" ]]; then
    echo "⚠️ ${base}.FS.pe.sam is empty — skipping."
    continue
fi
    samtools view -b -@ 8 "${base}.FS.pe.sam" | samtools sort -@ 8 -o "${base}.FS.sorted.bam"
    samtools index "${base}.FS.sorted.bam"
    samtools collate -o "${base}.FS.nameCol.bam" "${base}.FS.sorted.bam"
    samtools fixmate -m "${base}.FS.nameCol.bam" "${base}.FS.fixMate.bam"
    samtools sort -o "${base}.FS.pos.bam" "${base}.FS.fixMate.bam"
    samtools markdup -r -s "${base}.FS.pos.bam" "${base}.FS.markDup.bam"
    echo "✅ Done with FS $base"
done

# get stats
echo "📊 Generating flagstat summaries..."
for bam in *.markDup.bam; do
    sample="${bam%.markDup.bam}"
    echo "🔍 Stats for $sample"
    samtools flagstat "$bam" > "${sample}_flagstat.txt"
done

########################################################################################################
############################ 4.3 Test out MTaxi to see about competitive mapping##########################
########################################################################################################

#!/bin/bash
#SBATCH -J MTaxi
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
