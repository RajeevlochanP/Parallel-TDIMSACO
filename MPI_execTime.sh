
indices=(500)


echo "Compiling"
g++ -std=c++17 main.cpp -o main
g++ -std=c++17 -fopenmp mainOpenMP.cpp -o mainOpenMP
g++ -std=c++17 -fopenmp mainOpenMPTest.cpp -o mainOpenMPTest
mpic++ mainMPI.cpp -o mainMPI
mpic++ -fopenmp MPI_OpenMP.cpp -o MPI_OpenMP


arr=(./main ./mainOpenMP ./mainOpenMPTest ./mainMPI ./MPI_OpenMP)

gridFile="./gridGenration/grids.txt"
np=10

echo ""
echo "Starting Benchmark on ${#indices[@]} grids"


for exe in "${arr[@]}"; do
    

    times=()
    costs=()


    for idx in "${indices[@]}"; do
        
        if [[ "$exe" == *"MPI"* ]]; then
            cmd="mpirun -np $np $exe $gridFile $idx"
        else
            cmd="$exe $gridFile $idx"
        fi

        output=$($cmd)

        exec_time=$(echo "$output" | grep "Execution time " | awk -F': ' '{print $2}')
        sol_cost=$(echo "$output" | grep "SolutionCost " | awk -F': ' '{print $2}')

  
        if [[ -n "$exec_time" ]]; then times+=($exec_time); else times+=(0); fi
        if [[ -n "$sol_cost" ]]; then costs+=($sol_cost); else costs+=(0); fi
    done

    clean_name=$(basename "$exe" .out)
    
    echo "Results for $clean_name"
    echo "${clean_name}_times = [$(IFS=,; echo "${times[*]}")]"
    echo "${clean_name}_costs = [$(IFS=,; echo "${costs[*]}")]"

done