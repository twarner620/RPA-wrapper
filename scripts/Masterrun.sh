#!/usr/bin/env bash

###################MODIFY AS PER THE BATCH RUNS ###################################
ref='human-genome/non-altchromosomes.fa'  # Fasta file for the reference genome
range='human-genome/nucl.bed'             # BED file of ranges in the nucleus
bed='ribodataset/bedfolder/'              # Folder of bed files with extension .bed
order='config/order'                      # Order file for sample plotting sequence
outdir='run1'                             # Output directory (change per run to keep results separate)

###################DO NOT MODIFY BELOW THIS LINE ###################################
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
scripts="$SCRIPT_DIR/../RPA"

source "$SCRIPT_DIR/Heatmapwrapper.sh"

# Force Python to flush stdout after every write so log files stay current
export PYTHONUNBUFFERED=1

# Resolve input paths to absolute so they remain valid after cd
[[ "$ref"   != /* ]] && ref="$PWD/$ref"
[[ "$range" != /* ]] && range="$PWD/$range"
[[ "$bed"   != /* ]] && bed="$PWD/$bed"
[[ "$order" != /* ]] && order="$PWD/$order"

# bg_freq output is identical for the same ref+range — compute once and reuse across runs.
# All three strand variants write to distinct filenames so they can be computed in parallel.
BGFREQ_DIR="$PWD/bg_freq"
RANGE_BASE="$(basename "$range" .bed)"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

# Create all output and log directories up front so bg_freq subshells can log immediately
mkdir -p "$outdir/log" "$outdir/both/log" "$outdir/same/log" "$outdir/opp/log"

(
    exec > "$outdir/log/bg_freq_$TIMESTAMP.log" 2>&1
    if [[ ! -f "$BGFREQ_DIR/$RANGE_BASE.di.freq" ]]; then
        echo "[$(date)] Computing bg_freq"
        bg_freq "$scripts" "$ref" "$range"
    else
        echo "[$(date)] bg_freq already computed, skipping"
    fi
) &
(
    exec > "$outdir/log/bg_freq_ss_$TIMESTAMP.log" 2>&1
    if [[ ! -f "$BGFREQ_DIR/${RANGE_BASE}_same.di.freq" ]]; then
        echo "[$(date)] Computing bg_freq_ss"
        bg_freq_ss "$scripts" "$ref" "$range"
    else
        echo "[$(date)] bg_freq_ss already computed, skipping"
    fi
) &
(
    exec > "$outdir/log/bg_freq_os_$TIMESTAMP.log" 2>&1
    if [[ ! -f "$BGFREQ_DIR/${RANGE_BASE}_opp.di.freq" ]]; then
        echo "[$(date)] Computing bg_freq_os"
        bg_freq_os "$scripts" "$ref" "$range"
    else
        echo "[$(date)] bg_freq_os already computed, skipping"
    fi
) &

wait  # all bg_freq variants must finish before downstream steps

# Run all three strand workflows in parallel.
# Each subdirectory gets an absolute symlink to the shared bg_freq/ cache so
# the Heatmapwrapper functions find their background files at ./bg_freq/ as expected.
(
    exec > "$outdir/both/log/$TIMESTAMP.log" 2>&1
    cd "$outdir/both"
    ln -sfn "$BGFREQ_DIR" bg_freq
    echo "[$(date)] [both] Starting sample_freq"
    sample_freq "$scripts" "$ref" "$range" "$bed"
    echo "[$(date)] [both] Starting norm_freq"
    norm_freq "$scripts" "$ref" "$range" "$bed"
    echo "[$(date)] [both] Starting resort_plot"
    resort_plot "$scripts" "$ref" "$range" "$bed" "$order"
    echo "[$(date)] [both] Done"
) &

(
    exec > "$outdir/same/log/$TIMESTAMP.log" 2>&1
    cd "$outdir/same"
    ln -sfn "$BGFREQ_DIR" bg_freq
    echo "[$(date)] [same] Starting sample_freq_ss"
    sample_freq_ss "$scripts" "$ref" "$range" "$bed"
    echo "[$(date)] [same] Starting norm_freq_ss"
    norm_freq_ss "$scripts" "$ref" "$range" "$bed"
    echo "[$(date)] [same] Starting resort_plot_ss"
    resort_plot_ss "$scripts" "$ref" "$range" "$bed" "$order"
    echo "[$(date)] [same] Done"
) &

(
    exec > "$outdir/opp/log/$TIMESTAMP.log" 2>&1
    cd "$outdir/opp"
    ln -sfn "$BGFREQ_DIR" bg_freq
    echo "[$(date)] [opp] Starting sample_freq_os"
    sample_freq_os "$scripts" "$ref" "$range" "$bed"
    echo "[$(date)] [opp] Starting norm_freq_os"
    norm_freq_os "$scripts" "$ref" "$range" "$bed"
    echo "[$(date)] [opp] Starting resort_plot_os"
    resort_plot_os "$scripts" "$ref" "$range" "$bed" "$order"
    echo "[$(date)] [opp] Done"
) &

wait
echo "[$(date)] All workflows complete"
