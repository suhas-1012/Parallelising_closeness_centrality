#!/bin/bash
set -e

GRAPH_FILE=${1:-biconnected_graph.csv}
THREADS=${2:-16}

# Use associative array or just an ordered array Since bash 4+ we have arrays
FILES=(
    "01_seq_naive_bfs.cpp"
    "02_seq_multisource_bfs.cpp"
    "03_seq_bcc_reduced.cpp"
    "04_par_simple_mimd.cpp"
    "05_par_mimd_msbfs.cpp"
    "06_par_level_sync.cpp"
    "07_dynamic_shukla.cpp"
    "08_novel_bcc_spmm.cpp"
    "09_novel_vdbcc.cpp"
)

echo "Compiling..."
for f in "${FILES[@]}"; do
    bin_name="${f%.cpp}.bin"
    echo "  -> compiling $f to $bin_name"
    # we use -fopenmp for all, as 04, 05, 06, 08 need it
    g++ -std=c++17 "$f" -o "$bin_name" -fopenmp
done
echo "Compilation complete."
echo "=================================================="
echo "Comparing Performance on Graph: $GRAPH_FILE"
echo "Threads (where applicable): $THREADS"
echo "=================================================="

for f in "${FILES[@]}"; do
    bin_name="${f%.cpp}.bin"
    echo "--- Running $f ---"
    
    # 04, 05, 08, 09 might use threads. Others might ignore but we pass it anyway
    ./"$bin_name" "$GRAPH_FILE" "$THREADS"
    echo ""
done

echo "=================================================="
echo "All runs complete."
