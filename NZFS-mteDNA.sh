#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J nzfs-mt-eDNA
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --output nueDNA.%j.out # CHANGE map1 part each run
#SBATCH --error nueDNA.%j.err # CHANGE map1 part each run

## mt-eDNA

cd /nesi/nobackup/uoo02328/clare/Captured-eDNA/mtDNA/NZFS/NewTest2 || exit

## rename files

for filename in *.fastq.gz; do
    [ -f "$filename" ] || continue
    cp "$filename" "${filename//6428P2-??-0-1_/}"
done

## trim off adapters

module load cutadapt/2.10-gimkl-2020a-Python-3.8.2
module load pigz/2.4-GCCcore-7.4.0

for sample in S*.fastq.gz
do
base=$(basename ${sample} .fastq.gz)
echo "${base}"
cutadapt -a AGATCGGAAGAGCA -j 8 -o "${base}".trim.fastq.gz "${sample}"
done

## Quality filtering

module load Trimmomatic/0.39-Java-1.8.0_144

for trimmed in S*R1_.trim.fastq.gz
do
base=$(basename ${trimmed} _R1_.trim.fastq.gz)
echo "${base}"
trimmomatic PE -threads 5 -basein $"{trimmed}" -baseout "${base}".fq.gz SLIDINGWINDOW:10:20 MINLEN:30
done

# we use trimmomatic to use 5 threads for paired-end sequencing to trim F and R reads and output Sample10 (Sample10_1P.fq.gz and Sample10_2P.fq.gz)
# and then we say "OK sliding window" trim off anything that is less than a quality score of 20 within a sliding window of 10
# and THEN we say "Minimum length we want is 30 basepairs"

#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J nzfs-mt-eDNA
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --output mteDNA.%j.out # CHANGE map1 part each run
#SBATCH --error mteDNA.%j.err # CHANGE map1 part each run

## mt-eDNA

cd /nesi/nobackup/uoo02328/clare/Captured-eDNA/mtDNA/NZFS/NewTest2 || exit

## map samples to a genome: Arctocephalus forsteri mitogenome
#nzfs-mtDNAgenome-NC23.fasta

# map
module load cutadapt/2.10-gimkl-2020a-Python-3.8.2
module load pigz/2.4-GCCcore-7.4.0
module load BWA/0.7.17-gimkl-2017a
module load SAMtools/1.10-GCC-9.2.0

bwa index nzfs-mtDNAgenome-NC23.fasta

for sample in S*1P.fq.gz
do
base=$(basename ${sample} _1P.fq.gz)
echo "${base}"
bwa mem -t 4 nzfs-mtDNAgenome-NC23.fasta "${base}"_1P.fq.gz "${base}"_2P.fq.gz > "${base}".pe.sam
samtools view -b -@ 4 "${base}".pe.sam | samtools sort -@ 4 -o "${base}".sorted.bam
samtools index "${base}".sorted.bam
samtools collate -o "${base}".nameCol.bam "${base}".sorted.bam
samtools fixmate -m "${base}".nameCol.bam "${base}".fixMate.bam
samtools sort -o "${base}".pos.bam "${base}".fixMate.bam
samtools markdup -r -s "${base}".pos.bam "${base}".markDup.bam
done

## alternatively  for removing duplicates


#picard MarkDuplicates -I S10.bait.sorted.bam -O marked_dups.bam -M marked_dups_metrics.txt -REMOVE_SEQUENCING_DUPLICATES true


#!/bin/bash -e

#SBATCH -A uoo02328
#SBATCH --job-name      BLAST
#SBATCH --time          24:00:00  # ~10 CPU minutes / MB blastn query vs nt
#SBATCH --mem           48G
#SBATCH -c 8
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --output blast.%j.out # CHANGE map1 part each run
#SBATCH --error blast.%j.err # CHANGE map1 part each run

module load BLAST/2.10.0-GCC-9.2.0
module load BLASTDB/2021-01
module load SAMtools/1.10-GCC-9.2.0

#for file in *markDup.bam
#do
# base=$(basename ${file} .markDup.bam)
# samtools fasta "${file}" > "${base}".fasta
#done

# bam to fasta
samtools fasta S10.bait.markDup.bam > S10.bait.mapped.fasta


# This script takes one argument, the FASTA file of query sequences.
FORMAT="6 qseqid qstart qend qseq sseqid sgi sacc sstart send staxids sscinames stitle length evalue bitscore"
BLASTOPTS="-evalue 0.05 -max_target_seqs 10"
BLASTAPP=blastn
DB=nt
#BLASTAPP=blastx
#DB=nr

for file in *fasta
do
base=$(basename ${file} .fasta)
$BLASTAPP $BLASTOPTS -db $DB -query "${file}" -outfmt "$FORMAT" -out "${base}".$DB.$BLASTAPP -num_threads $SLURM_CPUS_PER_TASK
done
