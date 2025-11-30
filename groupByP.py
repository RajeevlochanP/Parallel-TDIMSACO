import subprocess
import re
import csv
import sys
import math
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

# Pairs of (Number of Ants, Number of Iterations)
TEST_CONFIGS = [
    (10, 20)
]

# --- NEW SETTING ---
# How many groups do you want between 0.0 and 1.0?
# Example: 5 groups -> ranges of 0.2 (0-0.2, 0.2-0.4...)
# Example: 10 groups -> ranges of 0.1 (0-0.1, 0.1-0.2...)
NUM_P_GROUPS = 5 

GRID_FILE_PATH = "./gridGenration/grids.txt"
GRID_RANGE = range(0, 50) 
OUTPUT_CSV = "groupByP.csv"

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
    data = {}
    
    # Regex to find P (Probability)
    p_match = re.search(r"p=([\d.]+)", output_str)
    if p_match:
        data['p'] = float(p_match.group(1))
    
    # Regex to find Time
    time_match = re.search(r"Time taken:\s*([\d.]+)", output_str)
    if time_match:
        data['time'] = float(time_match.group(1))
        
    return data

def get_bin_range(p_value, num_groups):
    """
    Calculates the start and end of the group this P value belongs to.
    """
    if num_groups <= 0: return (0.0, 1.0)
    
    step = 0.5 / num_groups
    
    # Calculate which index (0, 1, 2...) this P belongs to
    # e.g., p=0.3, step=0.2 -> index 1
    idx = int(p_value / step)
    
    # Edge case: if p=1.0, it might calculate index=num_groups, pull it back to last bin
    if idx >= num_groups:
        idx = num_groups - 1
        
    lower = idx * step
    upper = (idx + 1) * step
    
    # Return formatted tuple for sorting/display
    return (round(lower, 3), round(upper, 3))

# ==========================================
# MAIN EXECUTION
# ==========================================

def main():
    compile_code()
    
    # Structure: results[range_tuple][(filename, ants, iterations)] = [time1, time2, ...]
    aggregated_data = defaultdict(lambda: defaultdict(list))
    
    step_size = 1.0 / NUM_P_GROUPS
    print(f"--- Running Experiments for Grids {GRID_RANGE.start} to {GRID_RANGE.stop - 1} ---")
    print(f"--- Grouping into {NUM_P_GROUPS} sets (Step size: {step_size:.2f}) ---")

    for grid_id in GRID_RANGE:
        
        for (ants, iters) in TEST_CONFIGS:
            
            for file_info in SOURCE_FILES:
                exe = file_info['out']
                filename = file_info['src']
                
                try:
                    cmd = [exe, GRID_FILE_PATH, str(grid_id), str(ants), str(iters)]
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    parsed = parse_output(result.stdout)
                    
                    if 'p' not in parsed or 'time' not in parsed:
                        print(f"Warning: Could not parse output for {filename} on Grid {grid_id}")
                        continue
                    
                    # Calculate Bin
                    bin_range = get_bin_range(parsed['p'], NUM_P_GROUPS)
                    
                    # Store Time
                    config_key = (filename, ants, iters)
                    aggregated_data[bin_range][config_key].append(parsed['time'])
                    
                except Exception as e:
                    print(f"Error running {filename} on Grid {grid_id}: {e}")

        if grid_id % 10 == 0:
            print(f"Processed up to Grid {grid_id}...")

    # ==========================================
    # SUMMARIZE AND WRITE CSV
    # ==========================================
    print("\n--- Processing Data & Writing CSV ---")
    
    # Sort bins by lower bound
    sorted_bins = sorted(aggregated_data.keys(), key=lambda x: x[0])
    
    # Headers
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
                    # SUM execution times
                    total_time = sum(times)
                    row.append(f"{total_time:.6f}")
                else:
                    row.append("0.000000")
            
            writer.writerow(row)
            
    print(f"Successfully saved results to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()