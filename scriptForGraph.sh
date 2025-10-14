#!/bin/bash

SRC="main.cpp"
EXE="main"
GRID_FILE="./gridGenration/grids.txt"
CSV_FILE="./outputs/results.csv"

mkdir -p "$(dirname "$CSV_FILE")"

echo "Compiling $SRC..."
g++ -std=c++17 "$SRC" -o "$EXE"
if [ $? -ne 0 ]; then
    echo "❌ Compilation failed!"
    exit 1
fi

# Write CSV header
echo "index,n,p,solutionCost,timeTaken(s)" > "$CSV_FILE"

echo "Running $EXE for indexes 0 to 94..."
for i in $(seq 0 94); do
    START_TIME=$(date +%s.%N)

    OUTPUT=$(./"$EXE" "$GRID_FILE" "$i")

    END_TIME=$(date +%s.%N)
    TIME_TAKEN=$(echo "$END_TIME - $START_TIME" | bc)

    # Parse n, p, and SolutionCost using grep and regex
    N_VAL=$(echo "$OUTPUT" | grep -oP 'n=\K[0-9]+')
    P_VAL=$(echo "$OUTPUT" | grep -oP 'p=\K[0-9.]+')
    COST_VAL=$(echo "$OUTPUT" | grep -oP 'SolutionCost\s*:\s*\K[0-9.]+')

    # Append to CSV
    echo "$i,$N_VAL,$P_VAL,$COST_VAL,$TIME_TAKEN" >> "$CSV_FILE"

    echo "Run #$i → n=$N_VAL p=$P_VAL cost=$COST_VAL time=${TIME_TAKEN}s"
done

echo "✅ All runs complete. CSV saved to: $CSV_FILE"
