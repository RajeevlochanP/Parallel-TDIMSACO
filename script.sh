g++ -std=c++17 main.cpp -o main
g++ -std=c++17 -fopenmp mainOpenMP.cpp -o mainOpenMP
g++ -std=c++17 -fopenmp mainOpenMPTest.cpp -o mainOpenMPTest
mpic++ mainMPI.cpp -o mainMPI
mpic++ -fopenmp MPI_OpenMP.cpp -o MPI_OpenMP
echo "Running Tests:"
./main ./gridGenration/grids.txt 500 10 20 | grep "Time taken"
./mainOpenMP ./gridGenration/grids.txt 500 10 20 | grep "Time taken"
./mainOpenMPTest ./gridGenration/grids.txt 500 10 20 | grep "Time taken"
mpirun -np 10 ./mainMPI ./gridGenration/grids.txt 500 10 20 | grep "Time taken"
mpirun -np 10 ./MPI_OpenMP ./gridGenration/grids.txt 500 10 20 | grep "Time taken"
