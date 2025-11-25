#!/bin/bash

# indices=(11,30,71,91,118,141,157,187)
indices=(15 16 22 29 30);  

exe="./mainMPI.out"
gridFile="./gridGenration/grids.txt"

np=4

times=()
costs=()

echo "Starting Benchmark on ${#indices[@]} grids using $np processes..."

for idx in "${indices[@]}"; do
    output=$(mpirun -np $np $exe $gridFile $idx)

    exec_time=$(echo "$output" | grep "Execution time:" | awk -F': ' '{print $2}')
    sol_cost=$(echo "$output" | grep "SolutionCost :" | awk -F': ' '{print $2}')

    if [[ -n "$exec_time" && -n "$sol_cost" ]]; then 
        times+=($exec_time)
        costs+=($sol_cost)
    else
        echo "Index $idx       : Error occurred (Data not found)"
    fi
done
echo "execution_times = [$(IFS=,; echo "${times[*]}")]"
echo "solution_costs  = [$(IFS=,; echo "${costs[*]}")]"