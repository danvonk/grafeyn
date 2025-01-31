#!/bin/bash

EXAMPLE_CIRCUITS_DIR="../ExampleCircuits"
OUTPUT_DIR="$EXAMPLE_CIRCUITS_DIR/Outputs"
TMP_OUTPUT_FILE="$OUTPUT_DIR/tmp_output.txt"
SCHEDULERS=("naive" "gfq" "gnb" "nbsv")
INFO_FILE="$OUTPUT_DIR/specifications.txt"

# Save system and feyman specifications
CPU_MODEL=$(grep -m 1 "model name" /proc/cpuinfo | cut -d: -f2 | sed 's/^ *//')
CPU_CORES=$(nproc)
TOTAL_RAM=$(free -h | awk '/^Mem:/ {print $2}')
USED_RAM=$(free -h | awk '/^Mem:/ {print $3}')
GPU_INFO=$(lspci | grep -i --color 'vga\|3d\|2d')
OS_INFO=$(lsb_release -d | cut -f2-)
{
  echo "System Specifications"
  echo "====================="
  echo "Operating System: $OS_INFO"
  echo "CPU Model: $CPU_MODEL"
  echo "CPU Cores: $CPU_CORES"
  echo "Total RAM: $TOTAL_RAM"
  echo "Used RAM: $USED_RAM"
  echo "Graphics Card: $GPU_INFO"
  echo ""
  echo "Feyman Specifications"
  echo "====================="
  echo "Default Values:"
  echo "  Default Dense-Threshold: 0.25"
  echo "  Default Pull-Threshold: 0.8"
  echo "  Default Disable-Gate-Fusion: false"
  echo "  Default Simulator: parallel"
  echo "  Default Threads to use: 1"
  echo "  Default Block size for parallelism: 10000"
  echo "Changed Values:"
  echo "  Scheduler: all different types"
  for scheduler in "${SCHEDULERS[@]}"; do
    echo "    $scheduler"
  done
} > "$INFO_FILE"
echo "System and feyman specifications saved to $INFO_FILE"
echo ""

# Run with different Schedulers
mkdir -p "$OUTPUT_DIR"

echo "Starting Scheduler Comparison."

for file in "$EXAMPLE_CIRCUITS_DIR"/*; do
  if [ -f "$file" ]; then
    filename=$(basename "$file")

    echo "Processing: $file"

    for scheduler in "${SCHEDULERS[@]}"; do
      echo "Scheduler used: $scheduler"
      output_file="$OUTPUT_DIR/output_${filename}_$scheduler.txt"

      cargo run --package feynsum-rust --bin feynsum-rust -- --input "$file" --output "$TMP_OUTPUT_FILE" --scheduler "$scheduler" > "$output_file" 2>&1

      echo "Done Processing. Output in: $output_file"
    done

    echo ""
  fi
done

echo "Scheduler Comparison ended."
