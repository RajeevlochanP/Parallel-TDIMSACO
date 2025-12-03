#!/bin/bash

# Test Indices
indices=(500)

# 1. Compilation
echo "--- Compiling ---"
g++ -std=c++17 main.cpp -o main
g++ -std=c++17 -fopenmp mainOpenMP.cpp -o mainOpenMP
g++ -std=c++17 -fopenmp mainOpenMPTest.cpp -o mainOpenMPTest
mpic++ mainMPI.cpp -o mainMPI
mpic++ -fopenmp MPI_OpenMP.cpp -o MPI_OpenMP

# List of executables to run
arr=(./main ./mainOpenMP ./mainOpenMPTest ./mainMPI ./MPI_OpenMP)

gridFile="./gridGenration/grids.txt"
np=10

echo ""
echo "--- Starting Benchmark on ${#indices[@]} grids ---"

# Loop through each executable in the array
for exe in "${arr[@]}"; do
    
    # Create empty arrays for this specific executable
    times=()
    costs=()

    # Loop through each grid index
    for idx in "${indices[@]}"; do
        
        # Logic: If the executable name has "MPI", use mpirun. Otherwise, run directly.
        if [[ "$exe" == *"MPI"* ]]; then
            cmd="mpirun -np $np $exe $gridFile $idx"
        else
            cmd="$exe $gridFile $idx"
        fi

        # Run the command
        # echo "Running: $cmd" # Uncomment for debugging
        output=$($cmd)

        # Parse Output
        exec_time=$(echo "$output" | grep "Execution time:" | awk -F': ' '{print $2}')
        sol_cost=$(echo "$output" | grep "SolutionCost :" | awk -F': ' '{print $2}')

        # Store results (Default to 0 if parsing fails)
        if [[ -n "$exec_time" ]]; then times+=($exec_time); else times+=(0); fi
        if [[ -n "$sol_cost" ]]; then costs+=($sol_cost); else costs+=(0); fi
    done

    # Print results formatted like Python lists for easy copy-pasting
    # We clean the exe name for the variable (remove ./ and .out)
    clean_name=$(basename "$exe" .out)
    
    echo "------------------------------------------------"
    echo "Results for $clean_name"
    echo "${clean_name}_times = [$(IFS=,; echo "${times[*]}")]"
    echo "${clean_name}_costs = [$(IFS=,; echo "${costs[*]}")]"

done