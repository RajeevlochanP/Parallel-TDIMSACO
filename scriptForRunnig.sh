#!/bin/bash

SRC="main.cpp"
EXE="main"
GRID_FILE="./gridGenration/grids.txt"
NUM_RUNS=3

echo "Compiling $SRC with -pg for gprof..."
g++ -pg -std=c++17 "$SRC" -o "$EXE"
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

for i in $(seq 1 $NUM_RUNS); do
    INDEX=$i
    echo "Run #$i → index = $INDEX"

    ./"$EXE" "$GRID_FILE" "$INDEX" > "./outputs/output_$i.txt"

    gprof "$EXE" gmon.out > "./profilings/file${i}.txt"
    echo "Profiling result saved to file${i}.txt"

    rm gmon.out
done

echo "✅ All runs complete. Profiling results: file1.txt, file2.txt, file3.txt"

Script