#!/bin/bash

# Root directories
ROOT_DIR="$(pwd)"
INPUT_DIR="${ROOT_DIR}/input"
OUTPUT_DIR="${ROOT_DIR}/output"
GENOME_DIR="${ROOT_DIR}/genomes"
GENOME_INDEX="${GENOME_DIR}/INDEX/GRCh38"
GTF_FILE="${GENOME_DIR}/GTF/Homo_sapiens.GRCh38.109.gtf"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}/fastqc_results/raw"
mkdir -p "${OUTPUT_DIR}/fastqc_results/trimmed"
mkdir -p "${OUTPUT_DIR}/trimmomatic"
mkdir -p "${OUTPUT_DIR}/hisat2"
mkdir -p "${OUTPUT_DIR}/htseq"

# Path to Trimmomatic
TRIMMOMATIC_JAR="/home/satim/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTER_FILE="/home/satim/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"


# Loop through all FASTQ files
for R1 in "${INPUT_DIR}"/*_1.fastq.gz; do
    SAMPLE_NAME=$(basename "$R1" _1.fastq.gz)
    R2="${INPUT_DIR}/${SAMPLE_NAME}_2.fastq.gz"

    if [ ! -f "$R2" ]; then
        echo "Paired file for $R1 not found. Skipping..."
        continue
    fi

    echo "Processing sample: $SAMPLE_NAME"

    # Run FastQC on raw data
    echo "Running FastQC on raw FASTQ files..."
    fastqc "$R1" "$R2" -o "${OUTPUT_DIR}/fastqc_results/raw"

    # File paths for trimmed files
    TRIM_R1="${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_trimmed_1.fastq"
    TRIM_R2="${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_trimmed_2.fastq"
    TRIM_U1="${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_unpaired_1.fastq"
    TRIM_U2="${OUTPUT_DIR}/trimmomatic/${SAMPLE_NAME}_unpaired_2.fastq"

    # Trimming with Trimmomatic
    echo "Trimming reads with Trimmomatic..."
    java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
        "$R1" "$R2" \
        "$TRIM_R1" "$TRIM_U1" "$TRIM_R2" "$TRIM_U2" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Run FastQC on trimmed files
    echo "Running FastQC on trimmed FASTQ files..."
    fastqc "$TRIM_R1" "$TRIM_R2" -o "${OUTPUT_DIR}/fastqc_results/trimmed"

    # (Optional) Compress trimmed files
    gzip -f "$TRIM_R1"
    gzip -f "$TRIM_R2"
    TRIM_R1="${TRIM_R1}.gz"
    TRIM_R2="${TRIM_R2}.gz"

    # Align with HISAT2
    echo "Aligning reads with HISAT2..."
    SAM_FILE="${OUTPUT_DIR}/hisat2/${SAMPLE_NAME}.sam"
    hisat2 -x "$GENOME_INDEX" -1 "$TRIM_R1" -2 "$TRIM_R2" -S "$SAM_FILE"

    # Convert and sort BAM files
    echo "Converting and sorting BAM files..."
    BAM_FILE="${OUTPUT_DIR}/hisat2/${SAMPLE_NAME}.bam"
    SORTED_BAM="${OUTPUT_DIR}/hisat2/${SAMPLE_NAME}_sorted.bam"
    samtools view -b "$SAM_FILE" > "$BAM_FILE"
    samtools sort "$BAM_FILE" -o "$SORTED_BAM"
    #samtools index "$SORTED_BAM"
    rm "$SAM_FILE"  # Clean up

    # Count reads with HTSeq
    echo "Counting reads with HTSeq..."
    COUNTS_FILE="${OUTPUT_DIR}/htseq/${SAMPLE_NAME}_counts.txt"
    htseq-count -f bam -r pos -s no -i gene_id "$SORTED_BAM" "$GTF_FILE" > "$COUNTS_FILE"

    echo "Sample $SAMPLE_NAME processed successfully."
done

echo "All samples processed."

