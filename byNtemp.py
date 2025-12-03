import subprocess
import re
import csv
import sys
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
    }
]

TEST_CONFIGS = [
    (1, 2),
    (2, 4),
    (5, 10),
    (10, 20)
]

GRID_FILE_PATH = "./gridGenration/grids.txt"
GRID_RANGE = range(0, 71)
OUTPUT_CSV = "groupByN_MPI.csv"

def compile_code():
    print("Compiling Sources ")
    for file_info in SOURCE_FILES:
        cmd = [file_info['compiler']] + file_info['flags'] + [file_info['src'], "-o", file_info['out']]
        print(f"Compiling {file_info['src']}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error compiling {file_info['src']}:\n{result.stderr}")
            sys.exit(1)
    print("Compilation done\n")

def parse_output(output_str):
    data = {}
    
    n_match = re.search(r"n=(\d+)", output_str)
    if n_match:
        data['n'] = int(n_match.group(1))
    
    time_match = re.search(r"(?:Execution time|Time taken):\s*([\d.]+)", output_str)
    if time_match:
        data['time'] = float(time_match.group(1))
        
    sol_match = re.search(r"Solution\s*:\s*(.*)", output_str)
    if sol_match:
        data['path_str'] = sol_match.group(1).strip()
        
    return data


def main():
    compile_code()
    
    aggregated_data = defaultdict(lambda: defaultdict(list))
    
    print(f"Running experiments for grids {GRID_RANGE.start} to {GRID_RANGE.stop - 1} ")

    for grid_id in GRID_RANGE:
        current_n = None
        for (ants, iters) in TEST_CONFIGS:
            
            current_paths = {} 
            
            for file_info in SOURCE_FILES:
                exe = file_info['out']
                filename = file_info['src']
                run_prefix = file_info['run_cmd'] 
                
                try:
                    args = [exe, GRID_FILE_PATH, str(grid_id), str(ants), str(iters)]
                    
                    full_cmd = run_prefix + args
                    
                    result = subprocess.run(full_cmd, capture_output=True, text=True)
                    
                    output = result.stdout
                    parsed = parse_output(output)
                    
                    if 'n' not in parsed or 'time' not in parsed:
                        print(f"Could not parse output for {filename} on Grid {grid_id}")
                        continue
                    if current_n is None:
                        current_n = parsed['n']
                    elif current_n != parsed['n']:
                        print(f"Mismatch in grid size for grid {grid_id} (Expected {current_n}, got {parsed['n']})")
                    
                    config_key = (filename, ants, iters)
                    aggregated_data[current_n][config_key].append(parsed['time'])
                    current_paths[filename] = parsed.get('path_str', "")
                    
                except Exception as e:
                    print(f"Error running {filename} on Grid {grid_id} with Params ({ants},{iters}): {e}")

            if len(SOURCE_FILES) > 1:
                reference_file = SOURCE_FILES[0]['src']
                ref_path = current_paths.get(reference_file)
                
                for f_info in SOURCE_FILES:
                    fname = f_info['src']
                    other_path = current_paths.get(fname)
                    
                    if ref_path and other_path and ref_path != other_path:
                        print(f"Path mismatch on grid {grid_id} [Ants={ants}, Iter={iters}] between {reference_file} and {fname}")

        if grid_id % 10 == 0:
            print(f"Processed up to grid {grid_id}")

    print("\n Processing Data & Writing CSV")
    
    n_values = sorted(aggregated_data.keys())
    
    headers = ["N_Value"]
    column_keys = []
    
    for (ants, iters) in TEST_CONFIGS:
        for f in SOURCE_FILES:
            fname = f['src']
            headers.append(f"{fname} (A={ants}, I={iters})")
            column_keys.append((fname, ants, iters))
    
    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        
        for n in n_values:
            row = [n]
            for key in column_keys:
                times = aggregated_data[n].get(key, [])
                
                if times:
                    avg_time = sum(times) / len(times)
                    row.append(f"{avg_time:.6f}")
                else:
                    row.append("0.000000")
            writer.writerow(row)
            
    print(f"Successfully saved results to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()