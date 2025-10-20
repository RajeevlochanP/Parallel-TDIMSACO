# g++ -std=c++17 main.cpp -o main
# g++ -std=c++17 -fopenmp mainOpenMP.cpp -o mainOpenMP
g++ -std=c++17 -fopenmp mainOpenMPTest.cpp -o mainOpenMPTest
# ./main ./gridGenration/grids.txt 91
# ./mainOpenMP ./gridGenration/grids.txt 91
./mainOpenMPTest ./gridGenration/grids.txt 91

# python3 gridGeneration.py input.csv grids.txt --mode random_walk