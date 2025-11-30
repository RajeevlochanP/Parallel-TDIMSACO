SEQ_COMPILER="g++"
MPI_COMPILER="mpic++"
FLAGS="-std=c++17"

GRID_FILE="./gridGenration/grids.txt"
RAW_DATA_FILE="raw_data.tmp"
OUTPUT_CSV="final_results.csv"

GRID_START=0
GRID_END=99

CONFIGS=(
    "1 2"
    "2 4"
    "5 10"
    "10 20"
)
# ===============================================

# 1. Compile Codes
echo "--- Compiling ---"
$SEQ_COMPILER main.cpp -o main $FLAGS
if [ $? -ne 0 ]; then echo "Error compiling main.cpp"; exit 1; fi

$MPI_COMPILER mainMPI.cpp -o mainMPI $FLAGS # MPI usually doesn't need -std=c++17, but add if needed
if [ $? -ne 0 ]; then echo "Error compiling mainMPI.cpp"; exit 1; fi

echo "Compilation Successful."

# 2. Initialize Raw Data File (Headers)
# We store raw runs here first: N, Ants, Iters, TimeSeq, TimeMPI
echo "N,Ants,Iters,Time_Seq,Time_MPI" > $RAW_DATA_FILE

echo "--- Running Experiments (Grids $GRID_START to $GRID_END) ---"

# 3. Main Loops
for config in "${CONFIGS[@]}"; do
    # Split the config string "1 2" into variables
    set -- $config
    ants=$1
    iters=$2

    echo "Running Config: Ants=$ants, Iters=$iters..."

    for (( grid_id=$GRID_START; grid_id<=$GRID_END; grid_id++ )); do
        
        # --- Run Sequential ---
        # Capture output
        out_seq=$(./main "$GRID_FILE" "$grid_id" "$ants" "$iters")
        
        # Extract N (Grid Size) using grep and cut
        # Looks for "n=100" -> cuts delimiter "=" -> takes 2nd field "100"
        n_val=$(echo "$out_seq" | grep -o "n=[0-9]*" | cut -d= -f2)
        
        # Extract Time
        t_seq=$(echo "$out_seq" | grep -o "Time taken: [0-9.]*" | cut -d: -f2 | xargs)

        # --- Run MPI ---
        # Note: -np 4 is hardcoded here, change if needed
        out_mpi=$(mpirun -np 4 ./mainMPI "$GRID_FILE" "$grid_id" "$ants" "$iters")
        
        # Extract Time MPI
        t_mpi=$(echo "$out_mpi" | grep -o "Time taken: [0-9.]*" | cut -d: -f2 | xargs)

        # Safety Check: If empty (crashed), set to 0
        if [ -z "$t_seq" ]; then t_seq=0; fi
        if [ -z "$t_mpi" ]; then t_mpi=0; fi
        if [ -z "$n_val" ]; then n_val=-1; fi

        # Append to raw file
        echo "$n_val,$ants,$iters,$t_seq,$t_mpi" >> $RAW_DATA_FILE
    done
done

echo "--- Aggregating & Averaging Results ---"

# 4. Aggregate with AWK
# This script reads the raw file, groups by (N, Ants, Iters), sums times, and counts occurrences.
# Finally, it prints the Average.

awk -F, 'NR>1 {
    # Create a unique key for grouping
    key = $1 "," $2 "," $3
    
    # Accumulate sums and counts
    sum_seq[key] += $4
    sum_mpi[key] += $5
    count[key] += 1
}
END {
    print "N_Value,Ants,Iters,Avg_Time_Seq,Avg_Time_MPI"
    for (k in count) {
        if (count[k] > 0) {
            printf "%s,%.6f,%.6f\n", k, sum_seq[k]/count[k], sum_mpi[k]/count[k]
        }
    }
}' $RAW_DATA_FILE | sort -t, -k1n > $OUTPUT_CSV 
# 'sort -t, -k1n' sorts the CSV numerically by the first column (N_Value)

# Cleanup
rm $RAW_DATA_FILE

echo "Done! Results saved to $OUTPUT_CSV"