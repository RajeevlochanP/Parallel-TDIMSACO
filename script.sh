g++ -std=c++17 -O2 -fopenmp mainOpenMP.cpp -o mainOpenMP
g++ -std=c++17 -O2 -fopenmp main.cpp -o main
./main
./mainOpenMP

python3 gridGeneration.py input.csv grids.txt --mode random_walk