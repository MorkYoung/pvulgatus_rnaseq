#!/bin/bash
#SBATCH --job-name=rna_seq
#SBATCH --account=fc_wolflab
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --array=1-22

# SILVA 138 + Rfam v14.1 clustered rRNA refs (SortMeRNA v4.3.4 default db)
SORTMERNA_REF="../refs/sortmerna/smr_v4.3_default_db.fasta"

SAMPLE_FILE="../rna_seq/samples.csv"
# file organization
#treatment,rna_seq_path,genome_prefix,strain
SAMPLE_ROW=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_FILE}")
SAMPLE_ID=$(echo "${SAMPLE_ROW}" | cut -d ',' -f 1)
SAMPLE_RNA_SEQ_PATH=$(echo "${SAMPLE_ROW}" | cut -d ',' -f 2)
SAMPLE_RNA_SEQ_PATH='../rna_seq/'$SAMPLE_RNA_SEQ_PATH
SAMPLE_GENOME_PREFIX=$(echo "${SAMPLE_ROW}" | cut -d ',' -f 3)
SAMPLE_STRAIN=$(echo "${SAMPLE_ROW}" | cut -d ',' -f 4)


GENOME_FASTA="../ref_genomes/${SAMPLE_GENOME_PREFIX}_genomic.fna"
OUTDIR="../rna_seq/vulgatus_processed/${SAMPLE_ID}"
mkdir -p "${OUTDIR}"

source /global/scratch/projects/fc_wolflab/software/miniforge3/etc/profile.d/conda.sh
conda activate /global/scratch/projects/fc_wolflab/software/miniforge3/envs/mark_rnaseq

### TrimGalore
trim_galore --cores 4 --quality 20 --length 20 \
  "${SAMPLE_RNA_SEQ_PATH}" -o "${OUTDIR}"
TRIMMED=$(ls "${OUTDIR}"/*trimmed*.fq* | head -1)

### SortMeRNA
sortmerna --ref "${SORTMERNA_REF}" \
  --reads "${TRIMMED}" \
  --other "${OUTDIR}/${SAMPLE_ID}.norRNA.fastq.gz" \
  --aligned "${OUTDIR}/${SAMPLE_ID}.rrna.fastq.gz" \
  --fastx --log

### bowtie2
bowtie2 --threads 8 --very-sensitive-local \
  -x "${GENOME_FASTA%.fna}" \
  -U "${OUTDIR}/${SAMPLE_ID}.norRNA.fastq.gz" \
  2> "${OUTDIR}/${SAMPLE_ID}.bowtie2.log" \
  | samtools view -b -F 4 -F 256 -q 10 \
  | samtools sort -@ 4 -o "${OUTDIR}/${SAMPLE_ID}.mapped.bam"
samtools index "${OUTDIR}/${SAMPLE_ID}.mapped.bam"

### clean up 

gzip "${OUTDIR}/${SAMPLE_ID}.norRNA.fastq.gz
rm "${OUTDIR}/${SAMPLE_ID}.rrna.fastq.gz"
rm $TRIMMED