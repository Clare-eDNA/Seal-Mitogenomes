#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J nzfs-mt-eDNA
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --output mt-eDNA.%j.out # CHANGE map1 part each run
#SBATCH --error mt-eDNA.%j.err # CHANGE map1 part each run

## mt-eDNA

cd /nesi/nobackup/uoo02328/clare/Captured-eDNA/mtDNA/NZFS/NewTest3/onlyMapped || exit

## map samples to a genome: Arctocephalus forsteri mitogenome
#KT693367.fasta

# get some stats and prepare for coverage plot

module load BWA/0.7.17-gimkl-2017a
module load SAMtools/1.10-GCC-9.2.0


for sample in S*onlyMapped.bam
do
base=$(basename ${sample} .onlyMapped.bam)
echo "${base}"
echo "sample" "${base}" >> flagstat.txt
samtools flagstat "${sample}" >> flagstat.txt
samtools depth "${sample}" > "${base}".coverage
done

CF *
MOE 1
OHP *
SFB 3



S16
S1
S21
S22
S23
S24
S25
S26
S27
S28
S2
S34
S38
S3
S4
S5
S6
S7
S8 

