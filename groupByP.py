import subprocess
import re
import csv
import sys
import math
from collections import defaultdict

SOURCE_FILES = [
    {
        'compiler': 'g++',
        'src': 'main.cpp', 
        'out': './main', 
        'flags': ['-std=c++17'],
        'run_cmd': [] 
    },
    {
        'compiler': 'mpic++',
        'src': 'mainMPI.cpp', 
        'out': './mainMPI', 
        'flags': [],
        'run_cmd': ['mpirun', '-np', '4'] 
    },
]

TEST_CONFIGS = [
    (10, 20)
]

NUM_P_GROUPS = 5 

GRID_FILE_PATH = "./gridGenration/grids.txt"
GRID_RANGE = range(0, 99) 
OUTPUT_CSV = "groupByPMPI.csv"

def compile_code():
    for file_info in SOURCE_FILES:
        cmd = [file_info['compiler']] + file_info['flags'] + [file_info['src'], "-o", file_info['out']]
        print(f"Compiling {file_info['src']}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error compiling {file_info['src']}:\n{result.stderr}")
            sys.exit(1)
    print("Compilation successful\n")

def parse_output(output_str):
    data = {}
    p_match = re.search(r"p=([\d.]+)", output_str)
    if p_match:
        data['p'] = float(p_match.group(1))
    
    time_match = re.search(r"(?:Execution time|Time taken):\s*([\d.]+)", output_str)
    if time_match:
        data['time'] = float(time_match.group(1))
        
    return data

def get_bin_range(p_value, num_groups):
    if num_groups <= 0: return (0.0, 1.0)
    step = 0.5 / num_groups 
    idx = int(p_value / step)
    if idx >= num_groups:
        idx = num_groups - 1
    lower = idx * step
    upper = (idx + 1) * step
    return (round(lower, 3), round(upper, 3))

def main():
    compile_code()
    
    aggregated_data = defaultdict(lambda: defaultdict(list))
    
    step_size = 1.0 / NUM_P_GROUPS
    print(f"Running for grids {GRID_RANGE.start} to {GRID_RANGE.stop - 1}")
    
    for grid_id in GRID_RANGE:
        
        for (ants, iters) in TEST_CONFIGS:
            
            for file_info in SOURCE_FILES:
                exe = file_info['out']
                filename = file_info['src']
                run_prefix = file_info['run_cmd']
                
                try:
                    args = [exe, GRID_FILE_PATH, str(grid_id), str(ants), str(iters)]

                    full_cmd = run_prefix + args
                    result = subprocess.run(full_cmd, capture_output=True, text=True)
                    parsed = parse_output(result.stdout)
                    
                    if 'p' not in parsed or 'time' not in parsed:
                        print(f"Warning: Could not parse output for {filename} on Grid {grid_id}")
                        continue
                    
                    bin_range = get_bin_range(parsed['p'], NUM_P_GROUPS)
                    config_key = (filename, ants, iters)
                    aggregated_data[bin_range][config_key].append(parsed['time'])
                    
                except Exception as e:
                    print(f"Error running {filename} on Grid {grid_id}: {e}")

        if grid_id % 10 == 0:
            print(f"Processed up to grid {grid_id}")

    print("\nProcessing data & writing CSV")
    sorted_bins = sorted(aggregated_data.keys(), key=lambda x: x[0])
    headers = ["P_Range_Start", "P_Range_End", "Range_Label"]
    column_keys = []
    
    for (ants, iters) in TEST_CONFIGS:
        for f in SOURCE_FILES:
            fname = f['src']
            headers.append(f"{fname} (A={ants}, I={iters})")
            column_keys.append((fname, ants, iters))
    
    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        
        for bin_tuple in sorted_bins:
            start, end = bin_tuple
            label = f"{start:.2f}-{end:.2f}"
            
            row = [start, end, label]
            
            for key in column_keys:
                times = aggregated_data[bin_tuple].get(key, [])
                if times:
                    total_time = sum(times)
                    row.append(f"{total_time:.6f}")
                else:
                    row.append("0.000000")
            
            writer.writerow(row)
            
    print(f"Successfully saved results to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()