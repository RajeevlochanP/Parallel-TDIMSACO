import subprocess
import re
import csv
import sys
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================

SOURCE_FILES = [
    {
        'compiler': 'g++',
        'src': 'main.cpp', 
        'out': './main', 
        'flags': ['-std=c++17']
    },
    {
        'compiler': 'g++',
        'src': 'mainOpenMPTest.cpp', 
        'out': './mainOpenMPTest', 
        'flags': ['-std=c++17', '-fopenmp']
    },
    {
        'compiler': 'g++',
        'src': 'mainOpenMP.cpp', 
        'out': './mainOpenMP', 
        'flags': ['-std=c++17', '-fopenmp']
    }
]

# Define the pairs of (Number of Ants, Number of Iterations) you want to test
# Format: (ants, iterations)
TEST_CONFIGS = [
    (1, 2),
    (2, 4),
    (5, 10),
    (10, 20)
]

GRID_FILE_PATH = "./gridGenration/grids.txt"
GRID_RANGE = range(0, 100)  # Example: 0 to 499
OUTPUT_CSV = "groupByN.csv"

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def compile_code():
    print("--- Compiling Sources ---")
    for file_info in SOURCE_FILES:
        cmd = [file_info['compiler']] + file_info['flags'] + [file_info['src'], "-o", file_info['out']]
        print(f"Compiling {file_info['src']}...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error compiling {file_info['src']}:\n{result.stderr}")
            sys.exit(1)
    print("Compilation successful.\n")

def parse_output(output_str):
    """
    Extracts n, time, and solution path from the C++ output.
    """
    data = {}
    
    # Regex to find N (Size)
    n_match = re.search(r"n=(\d+)", output_str)
    if n_match:
        data['n'] = int(n_match.group(1))
    
    # Regex to find Time
    time_match = re.search(r"Time taken:\s*([\d.]+)", output_str)
    if time_match:
        data['time'] = float(time_match.group(1))
        
    # Regex to find Solution String: (0,0),(3,3)...
    sol_match = re.search(r"Solution\s*:\s*(.*)", output_str)
    if sol_match:
        data['path_str'] = sol_match.group(1).strip()
        
    return data

# ==========================================
# MAIN EXECUTION
# ==========================================

def main():
    compile_code()
    
    # Structure: results[n][(filename, ants, iterations)] = [time1, time2, ...]
    # We use a tuple key to distinguish results for specific parameters
    aggregated_data = defaultdict(lambda: defaultdict(list))
    
    print(f"--- Running Experiments for Grids {GRID_RANGE.start} to {GRID_RANGE.stop - 1} ---")

    for grid_id in GRID_RANGE:
        current_n = None
        
        # Iterate through every parameter configuration (Ants, Iterations)
        for (ants, iters) in TEST_CONFIGS:
            
            current_paths = {} # To compare paths ONLY within this specific config
            
            for file_info in SOURCE_FILES:
                exe = file_info['out']
                filename = file_info['src']
                
                try:
                    # Execute C++ program: ./exe grid_file grid_id ants iterations
                    # This maps to argv[0], argv[1], argv[2], argv[3], argv[4]
                    cmd = [exe, GRID_FILE_PATH, str(grid_id), str(ants), str(iters)]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    output = result.stdout
                    parsed = parse_output(output)
                    
                    if 'n' not in parsed or 'time' not in parsed:
                        print(f"Warning: Could not parse output for {filename} on Grid {grid_id}")
                        continue
                    
                    # Consistency Check: Ensure 'n' is consistent across files/params
                    if current_n is None:
                        current_n = parsed['n']
                    elif current_n != parsed['n']:
                        print(f"CRITICAL ERROR: Mismatch in grid size for Grid {grid_id} (Expected {current_n}, got {parsed['n']})")
                    
                    # Store Time with unique key combining File + Params
                    config_key = (filename, ants, iters)
                    aggregated_data[current_n][config_key].append(parsed['time'])
                    
                    # Store Path for comparison
                    current_paths[filename] = parsed.get('path_str', "")
                    
                except Exception as e:
                    print(f"Error running {filename} on Grid {grid_id} with Params ({ants},{iters}): {e}")

            # --- Path Consistency Check (Per Parameter Set) ---
            # We only compare if we have more than 1 file running
            if len(SOURCE_FILES) > 1:
                reference_file = SOURCE_FILES[0]['src']
                ref_path = current_paths.get(reference_file)
                
                for f_info in SOURCE_FILES:
                    fname = f_info['src']
                    other_path = current_paths.get(fname)
                    
                    # Only compare if both ran successfully
                    if ref_path and other_path and ref_path != other_path:
                        print(f"Warning: Path mismatch on Grid {grid_id} [Ants={ants}, Iter={iters}] between {reference_file} and {fname}")

        if grid_id % 10 == 0:
            print(f"Processed up to Grid {grid_id}...")

    # ==========================================
    # CALCULATE AVERAGES AND WRITE CSV
    # ==========================================
    print("\n--- Processing Data & Writing CSV ---")
    
    # Get all unique N values sorted
    n_values = sorted(aggregated_data.keys())
    
    # Generate Headers Dynamically
    # Column format: "Filename (A=X, I=Y)"
    headers = ["N_Value"]
    column_keys = []
    
    # We iterate configs then files to group columns logically
    for (ants, iters) in TEST_CONFIGS:
        for f in SOURCE_FILES:
            fname = f['src']
            header_str = f"{fname} (A={ants}, I={iters})"
            headers.append(header_str)
            column_keys.append((fname, ants, iters)) # Store key for data retrieval
    
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
                    row.append("N/A")
            writer.writerow(row)
            
    print(f"Successfully saved results to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()