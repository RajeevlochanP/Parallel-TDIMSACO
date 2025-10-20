#!/usr/bin/env python3
"""
generate_grids_from_csv.py

Usage:
    python generate_grids_from_csv.py input.csv grids.txt [--mode monotonic|random_walk] [--seed BASE_SEED]

CSV format: no header, each row: n,p  (example: 10,0.35)
Outputs grids.txt containing multiple grids in parseable format.

Grid file format (per-grid):
GRID <index> n=<n> p=<p> seed=<seed>
<row0: space-separated 0/1>
...
<rowN-1>

Example:
GRID 0 n=5 p=0.35 seed=123456789
0 1 0 0 1
0 0 0 1 0
...
"""
import sys
import csv
import random
import argparse
import secrets
import hashlib
import os

def generate_grid(n, p, seed=None, mode='random_walk', max_walk_attempts=1000):
    assert n >= 1
    assert 0.0 <= p <= 1.0
    if seed is not None:
        random.seed(seed)

    grid = [[1 if random.random() < p else 0 for _ in range(n)] for _ in range(n)]

    def mark_path(path):
        for (r, c) in path:
            grid[r][c] = 0

    if mode == 'monotonic':
        moves = ['R'] * (n - 1) + ['D'] * (n - 1)
        random.shuffle(moves)
        r = c = 0
        path = [(r, c)]
        for m in moves:
            if m == 'R':
                c += 1
            else:
                r += 1
            path.append((r, c))
        mark_path(path)

    elif mode == 'random_walk':
        targets = (n - 1, n - 1)
        directions = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        attempt = 0
        success = False
        while attempt < max_walk_attempts and not success:
            attempt += 1
            r = c = 0
            visited = {(r, c)}
            path = [(r, c)]
            while (r, c) != targets:
                nbrs = []
                for dr, dc in directions:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < n and 0 <= nc < n and (nr, nc) not in visited:
                        nbrs.append((nr, nc))
                if not nbrs:
                    break
                nbrs.sort(key=lambda x: abs(x[0]-targets[0]) + abs(x[1]-targets[1]))
                k = min(3, len(nbrs))
                choice = random.choice(nbrs[:k])
                r, c = choice
                visited.add((r, c))
                path.append((r, c))
            if (r, c) == targets:
                mark_path(path)
                success = True
        if not success:
            return generate_grid(n, p, seed=seed, mode='monotonic')
    else:
        raise ValueError("mode must be 'monotonic' or 'random_walk'")

    grid[0][0] = 0
    grid[n-1][n-1] = 0
    return grid

def main():
    parser = argparse.ArgumentParser(description="Generate grids from CSV.")
    parser.add_argument("input_csv")
    parser.add_argument("output_grids", help="output grids text file (e.g. grids.txt)")
    parser.add_argument("--mode", choices=["monotonic", "random_walk"], default="monotonic")
    parser.add_argument("--seed", type=int, default=None, help="optional base seed (64-bit int). If omitted a high-quality seed is chosen.")
    args = parser.parse_args()

    if args.seed is None:
        seed_base = int.from_bytes(hashlib.sha256(os.urandom(64)).digest()[:8], "big")
    else:
        seed_base = args.seed

    with open(args.input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rows = [row for row in reader if row and len(row) >= 2]

    if not rows:
        print("No rows found in input CSV.")
        return

    with open(args.output_grids, "w") as outf:
        for idx, row in enumerate(rows):
            try:
                n = int(row[0])
                p = float(row[1])
            except Exception as e:
                print(f"Skipping row {idx} ({row}): parse error {e}")
                continue
            per_seed = (seed_base + idx) & ((1<<63)-1)
            grid = generate_grid(n, p, seed=per_seed, mode=args.mode)
            outf.write(f"GRID {idx} n={n} p={p} seed={per_seed}\n")
            for r in range(n):
                outf.write(" ".join(str(x) for x in grid[r]) + "\n")
            outf.write("\n")
    print(f"Wrote {len(rows)} grids to {args.output_grids} with base seed {seed_base}")

if __name__ == "__main__":
    main()
