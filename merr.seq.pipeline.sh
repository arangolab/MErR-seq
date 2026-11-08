#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --array=0-7%8
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --job-name=merrseq
#SBATCH --output=merrseq_%A_%a.out
#SBATCH --error=merrseq_%A_%a.err

set -euo pipefail

echo "========== PIPELINE START =========="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
date

# -----------------------------
# Input / output directories
# -----------------------------
# ***For new project: change BASE path and number of arrays ***

BASE="/projects/b1042/Arangolab/Daniel/MErRseq"

INPUT_DIR="$BASE/fastq"
OUTPUT_CUTADAPT="$BASE/cutadapt"
OUTPUT_UMI="$BASE/umi_tools/umi_extract"
OUTPUT_47S="$BASE/Bowtie/47S"
OUTPUT_BIGWIG="$BASE/Bigwig"

mkdir -p "$OUTPUT_CUTADAPT" "$OUTPUT_UMI" "$OUTPUT_47S" "$OUTPUT_BIGWIG"

# -----------------------------
# Sample Selection
# -----------------------------
FILE_EXT="_1.fq.gz"

mapfile -d '' files < <( 
    find -L "$INPUT_DIR" -maxdepth 1 -type f -name "*${FILE_EXT}" -print0 | sort -z
)

echo "Found ${#files[@]} R1 files"

if [[ "$SLURM_ARRAY_TASK_ID" -ge "${#files[@]}" ]]; then
    echo "ERROR: Task ID $SLURM_ARRAY_TASK_ID exceeds number of files (${#files[@]})"
    exit 1
fi

R1="${files[$SLURM_ARRAY_TASK_ID]}"
R2="${R1%_1.fq.gz}_2.fq.gz"   


if [[ ! -f "$R2" ]]; then
    echo "ERROR: Matching R2 not found"
    echo "R1: $R1"
    echo "Expected R2: $R2"
    exit 1
fi

BASENAME=$(basename "$R1")
PREFIX="${BASENAME%_1.fq.gz}"     # A1_CKDL250030627-1A_2333KWLT3_L3 (log only)
SAMPLE="${BASENAME%%_*}"          # A1 (used for ALL output names)

echo "Processing sample: $SAMPLE ($PREFIX)"
echo "R1 = $R1"
echo "R2 = $R2"

# -----------------------------
# Cutadapt
# -----------------------------
CUTADAPT_R1="$OUTPUT_CUTADAPT/${SAMPLE}_cutadapt_1.fq.gz"
CUTADAPT_R2="$OUTPUT_CUTADAPT/${SAMPLE}_cutadapt_2.fq.gz"
CUTADAPT_METRICS="$OUTPUT_CUTADAPT/${SAMPLE}_cutadapt.metrics"

echo "---------- CUTADAPT ----------"

if [[ -s "$CUTADAPT_R1" && -s "$CUTADAPT_R2" ]]; then
    echo "Cutadapt outputs already exist for $SAMPLE"
    echo "Skipping cutadapt"
else
    module purge
    module load cutadapt/4.2

    if cutadapt \
        -j 8 \
        --match-read-wildcards \
        --times 3 \
        -e 0.1 \
        -O 3 \
        --nextseq-trim=10 \
        -m 50 \
        -a AGATCGGAAGAGCACA \
        -a ACACGAGAAGGCTAGA \
        -a TGTGCTCTTCCGATCT \
        -a TCTAGCCTTCTCGTGT \
        -a ACGACGCTCTTCCGATCT \
        -a TCTAGCCTTCTCGCAGCA \
        -a AGATCGGAAGAGCGTCGT \
        -a TGCTGCGAGAAGGCTAGA \
        -a G{10} \
        -a A{10} \
        -A AGATCGGAAGAGCACA \
        -A ACACGAGAAGGCTAGA \
        -A TGTGCTCTTCCGATCT \
        -A TCTAGCCTTCTCGTGT \
        -A ACGACGCTCTTCCGATCT \
        -A TCTAGCCTTCTCGCAGCA \
        -A AGATCGGAAGAGCGTCGT \
        -A TGCTGCGAGAAGGCTAGA \
        -A G{10} \
        -A A{10} \
        -o "$CUTADAPT_R1" \
        -p "$CUTADAPT_R2" \
        "$R1" "$R2" \
        > "$CUTADAPT_METRICS"
    then
        echo "Cutadapt succeeded"
    else
        echo "ERROR: Cutadapt failed"
        exit 1
    fi

    if [[ ! -s "$CUTADAPT_R1" || ! -s "$CUTADAPT_R2" ]]; then
        echo "ERROR: Cutadapt output files missing"
        exit 1
    fi
fi

# -----------------------------
# UMI-tools extract
# -----------------------------
echo "---------- UMI EXTRACT ----------"

R1_UMI="$OUTPUT_UMI/${SAMPLE}_umi_1.fq.gz"
R2_UMI="$OUTPUT_UMI/${SAMPLE}_umi_2.fq.gz"

if [[ -s "$R1_UMI" && -s "$R2_UMI" ]]; then
    echo "UMI outputs already exist for $SAMPLE"
    echo "Skipping umi_tools extract"
else
    set +u
    eval "$(conda shell.bash hook)"
    conda activate /projects/p31544/envs/umi-env
    set -u

    if umi_tools extract \
        --extract-method=string \
        --bc-pattern=NNNNNNNNNN \
        --stdin="$CUTADAPT_R1" \
        --read2-in="$CUTADAPT_R2" \
        --stdout="$R1_UMI" \
        --read2-out="$R2_UMI"
    then
        echo "UMI extraction succeeded"
    else
        echo "ERROR: UMI extraction failed"
        exit 1
    fi

    if [[ ! -s "$R1_UMI" || ! -s "$R2_UMI" ]]; then
        echo "ERROR: UMI output files missing"
        exit 1
    fi

    set +u
    conda deactivate
    set -u
fi

# ============================================================
# Bowtie2 rRNA ALIGNMENT - 47S
# ============================================================
echo "---------- BOWTIE2 rRNA ALIGNMENT (47S) ----------"

INDEX_47S="/projects/p31544/genome_indexes/Bowtie/rRNA_processing/47S/47S"

SAM_47S="$OUTPUT_47S/${SAMPLE}_47S.sam"
BAM_47S="$OUTPUT_47S/${SAMPLE}_47S_sorted.bam"
METRICS_47S="$OUTPUT_47S/${SAMPLE}_47S_metrics.txt"
LOG_47S="$OUTPUT_47S/${SAMPLE}_47S_bowtie2.log"

if [[ -s "$BAM_47S" && -s "${BAM_47S}.bai" ]]; then
    echo "47S alignment already complete for $SAMPLE -> skipping"
else
    module purge
    module load bowtie2/2.5.5

    if bowtie2 \
        -p 8 -q --local \
        --no-unal \
        -x "$INDEX_47S" \
        -1 "$R1_UMI" -2 "$R2_UMI" \
        -S "$SAM_47S" \
        --met-file "$METRICS_47S" \
        2> "$LOG_47S"
    then
        echo "Bowtie2 47S succeeded"
        cat "$LOG_47S"
    else
        echo "ERROR: Bowtie2 47S failed"
        cat "$LOG_47S"
        exit 1
    fi

    module purge
    module load samtools/1.6

    if samtools sort -@ 8 -o "$BAM_47S" "$SAM_47S" && samtools index "$BAM_47S"; then
        echo "47S sort + index succeeded"
    else
        echo "ERROR: 47S sort/index failed"
        exit 1
    fi

    rm -f "$SAM_47S"
fi

# ============================================================
# UMI DEDUP with UMICollapse (replaces umi_tools dedup)
# ============================================================
# Strategy (memory-safe for rRNA, where all reads pile on one contig):
#   1. Take R1 of properly paired, primary alignments
#   2. Deduplicate R1 in single-end mode (R1 position + UMI)
#      -> no --paired, so UMICollapse never buffers mates in memory
#   3. Keep BOTH mates of every surviving R1 (filter by read name)
# Output stays coordinate-sorted, so it only needs indexing.


# UMICollapse settings

# Conda env created ONCE beforehand (not inside the array job):
#   module load python-miniconda3/4.12.0
#   conda create -p /projects/p31544/envs/umicollapse -c conda-forge -c bioconda \
#       --override-channels --strict-channel-priority umicollapse -y
# samtools comes from the cluster module (samtools/1.6), not from this env.
UMICOLLAPSE_ENV="/projects/p31544/envs/umicollapse"

# Java heap for UMICollapse
UMICOLLAPSE_HEAP="$(( ${SLURM_MEM_PER_NODE:-131072} * 85 / 100 ))m"

# Edit distance for merging UMIs:
#   1 = directional clustering, tolerates 1 sequencing error in the 10-nt UMI (recommended)
#   0 = exact match only (equivalent to your previous umi_tools --method=unique)
UMI_EDIT_DIST=0

# Usage: run_umicollapse <input_sorted.bam> <output_dedup.bam> <label>
# Writes a coordinate-sorted, indexed dedup BAM (both mates) plus a log and a
# before/after fragment-count summary (replaces --output-stats).

run_umicollapse () {
    local IN_BAM="$1"
    local OUT_BAM="$2"
    local LABEL="$3"
    local STEM="${OUT_BAM%.bam}"
    local R1_BAM="${STEM}.R1.tmp.bam"
    local R1_DEDUP="${STEM}.R1dedup.tmp.bam"
    local KEEP="${STEM}.keep_names.tmp.txt"
    local LOG="${STEM}.umicollapse.log"
    local SUMMARY="${STEM}.summary.txt"

    echo "---------- UMI DEDUP ($LABEL) ----------"

    if [[ -s "$OUT_BAM" && -s "${OUT_BAM}.bai" ]]; then
        echo "$LABEL dedup BAM already exists, skipping"
        return 0
    fi

    module purge
    module load python-miniconda3/4.12.0
    module load samtools/1.6
    set +u
    eval "$(conda shell.bash hook)"
    conda activate "$UMICOLLAPSE_ENV"
    set -u

    # --- 1. R1 of properly paired primary alignments ---
    #   -f 67    = paired (1) + proper pair (2) + first in pair (64)
    #   -F 0x90C = drop unmapped, mate unmapped, secondary, supplementary
    #   (replaces --remove-unpaired / --remove-chimeric)
    if samtools view -@ 8 -b -f 67 -F 0x90C -o "$R1_BAM" "$IN_BAM" \
        && samtools index "$R1_BAM"
    then
        echo "$LABEL R1 extracted"
    else
        echo "ERROR: $LABEL R1 extraction failed"
        exit 1
    fi

    # --- 2. Single-end UMI dedup on R1 ---
    # JAVA_TOOL_OPTIONS sets the JVM heap regardless of the wrapper script
    if JAVA_TOOL_OPTIONS="-Xmx${UMICOLLAPSE_HEAP}" umicollapse bam \
        -i "$R1_BAM" \
        -o "$R1_DEDUP" \
        --umi-sep _ \
        --algo dir \
        -k "$UMI_EDIT_DIST" \
        --two-pass \
        > "$LOG" 2>&1
    then
        echo "UMICollapse $LABEL succeeded"
    else
        echo "ERROR: UMICollapse $LABEL failed (see $LOG)"
        cat "$LOG"
        exit 1
    fi

    # --- 3. Names of surviving fragments ---
    samtools view "$R1_DEDUP" | cut -f1 > "$KEEP"

    if [[ ! -s "$KEEP" ]]; then
        echo "ERROR: no reads survived $LABEL dedup"
        exit 1
    fi

    # --- 4. Keep both mates of surviving fragments (order preserved = still sorted) ---
    if samtools view -@ 8 -h -f 3 -F 0x90C "$IN_BAM" \
        | awk 'BEGIN{FS="\t"} NR==FNR {keep[$1]; next} /^@/ || ($1 in keep)' "$KEEP" - \
        | samtools view -@ 8 -b -o "$OUT_BAM" - \
        && samtools index "$OUT_BAM"
    then
        echo "$LABEL dedup BAM (both mates) written + indexed"
    else
        echo "ERROR: $LABEL mate filtering failed"
        exit 1
    fi

    # --- Lightweight stats: fragments (R1) before vs after ---
    local BEFORE AFTER
    BEFORE=$(samtools view -c "$R1_BAM")
    AFTER=$(wc -l < "$KEEP")
    {
        echo -e "sample\tset\tfragments_before\tfragments_after\tfraction_kept"
        awk -v s="$SAMPLE" -v l="$LABEL" -v b="$BEFORE" -v a="$AFTER" \
            'BEGIN{printf "%s\t%s\t%d\t%d\t%.4f\n", s, l, b, a, (b>0 ? a/b : 0)}'
    } > "$SUMMARY"
    cat "$SUMMARY"

    rm -f "$R1_BAM" "${R1_BAM}.bai" "$R1_DEDUP" "$KEEP"

    set +u
    conda deactivate
    set -u
}

DEDUP_BAM_47S="$OUTPUT_47S/${SAMPLE}_47S_dedup.bam"

run_umicollapse "$BAM_47S"   "$DEDUP_BAM_47S"   "47S"

# ============================================================
# R1 / R2 SPLIT + BIGWIG - 47S
# ============================================================
echo "---------- R1/R2 + BIGWIG (47S) ----------"

for READ in R1 R2; do

    if [[ "$READ" == "R1" ]]; then FLAG=64; else FLAG=128; fi

    READ_BAM_47S="$OUTPUT_47S/${SAMPLE}_47S_dedup_${READ}.bam"
    BIGWIG_47S="$OUTPUT_BIGWIG/${SAMPLE}_47S_${READ}.bw"

    # --- extract R1 or R2 from dedup BAM ---
    if [[ -s "$READ_BAM_47S" && -s "${READ_BAM_47S}.bai" ]]; then
        echo "47S $READ BAM already exists, skipping"
    else
        module purge
        module load samtools/1.6

        if samtools view -@ 8 -b -f "$FLAG" -o "$READ_BAM_47S" "$DEDUP_BAM_47S" \
            && samtools index "$READ_BAM_47S"
        then
            echo "47S $READ extracted + indexed"
        else
            echo "ERROR: 47S $READ extraction failed"
            exit 1
        fi
    fi

    # --- bigwig ---
    if [[ -s "$BIGWIG_47S" ]]; then
        echo "47S $READ bigwig already exists, skipping"
    else
        module purge
        module load deeptools/3.5.6

        if bamCoverage \
            -b "$READ_BAM_47S" \
            -o "$BIGWIG_47S" \
            --binSize 1 \
            -p 8 \
            --normalizeUsing RPKM

        then
            echo "47S $READ bigwig created"
        else
            echo "ERROR: 47S $READ bigwig failed"
            exit 1
        fi
    fi

done

echo "========== PIPELINE END =========="
date