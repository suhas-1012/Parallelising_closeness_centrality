// sequential_cc.cpp — Phase 1: Sequential Closeness Centrality (Baseline)
// ======================================================================
// This is the SEQUENTIAL baseline. No parallelism at all.
// It implements:
//   1. Naive BFS (one BFS per source)
//   2. Multi-Source BFS using 64-bit bitmasks (bit-parallelism, NOT SIMD)
//   3. BCC decomposition + graph reduction (structural optimization)
//
// Compile: g++ -O2 -std=c++17 -o seq_cc sequential_cc.cpp
// Usage:   ./seq_cc [graph_file]
//          ./seq_cc              (generates random test graph)
//
// Graph file format: first line "n m", then m lines "u v" (0-indexed, undirected)
// ======================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <iomanip>
#include <stack>
#include <set>
#include <unordered_set>
#include <cstdint>

using namespace std;

// ============================================================
// Timer utility
// ============================================================
struct Timer {
    chrono::high_resolution_clock::time_point start;
    Timer() : start(chrono::high_resolution_clock::now()) {}
    double elapsed_ms() const {
        auto now = chrono::high_resolution_clock::now();
        return chrono::duration<double, milli>(now - start).count();
    }
};

// ============================================================
// CSR (Compressed Sparse Row) Graph — cache-friendly layout
// ============================================================
struct CSRGraph {
    int n;                       // number of nodes
    long long m;                 // number of edges (undirected, counted once)
    vector<long long> rowPtr;    // rowPtr[v]..rowPtr[v+1]-1 = neighbors of v
    vector<int> colIdx;          // colIdx[i] = neighbor node ID
    
    CSRGraph() : n(0), m(0) {}
    
    // Build from adjacency list
    CSRGraph(int n, const vector<vector<int>>& adj) : n(n) {
        rowPtr.resize(n + 1, 0);
        // Count degrees
        for (int u = 0; u < n; u++) {
            rowPtr[u + 1] = adj[u].size();
        }
        // Prefix sum
        for (int u = 0; u < n; u++) {
            rowPtr[u + 1] += rowPtr[u];
        }
        m = rowPtr[n] / 2; // each edge counted twice in undirected
        
        // Fill column indices
        colIdx.resize(rowPtr[n]);
        vector<long long> offset(n, 0);
        for (int u = 0; u < n; u++) {
            for (int v : adj[u]) {
                colIdx[rowPtr[u] + offset[u]++] = v;
            }
        }
    }
    
    // Iterate over neighbors of u
    int degree(int u) const { return (int)(rowPtr[u + 1] - rowPtr[u]); }
    const int* neighbors_begin(int u) const { return &colIdx[rowPtr[u]]; }
    const int* neighbors_end(int u) const { return &colIdx[rowPtr[u + 1]]; }
    
    // Read from file: first line "n m", then m lines "u v"
    static CSRGraph readFromFile(const string& filename) {
        ifstream fin(filename);
        if (!fin) {
            cerr << "Error: Cannot open " << filename << endl;
            exit(1);
        }
        int n, m;
        fin >> n >> m;
        vector<vector<int>> adj(n);
        for (int i = 0; i < m; i++) {
            int u, v;
            fin >> u >> v;
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
        return CSRGraph(n, adj);
    }
    
    // Generate a random connected graph for testing
    static CSRGraph generateRandom(int n, int m) {
        vector<vector<int>> adj(n);
        srand(42);
        // Spanning tree first (ensures connectivity)
        for (int i = 1; i < n; i++) {
            int parent = rand() % i;
            adj[i].push_back(parent);
            adj[parent].push_back(i);
        }
        // Add remaining random edges
        set<pair<int,int>> existing;
        for (int u = 0; u < n; u++)
            for (int v : adj[u])
                existing.insert({min(u,v), max(u,v)});
        
        int added = n - 1;
        while (added < m) {
            int u = rand() % n;
            int v = rand() % n;
            if (u == v) continue;
            auto e = make_pair(min(u,v), max(u,v));
            if (existing.count(e)) continue;
            existing.insert(e);
            adj[u].push_back(v);
            adj[v].push_back(u);
            added++;
        }
        return CSRGraph(n, adj);
    }
};

// ============================================================
// Method 1: Naive Sequential BFS — O(n * (n + m))
//   The simplest possible baseline. One BFS per source.
// ============================================================
vector<double> naiveSequentialCC(const CSRGraph& g) {
    int n = g.n;
    vector<double> cc(n, 0.0);
    vector<int> dist(n);
    
    for (int s = 0; s < n; s++) {
        // BFS from source s
        fill(dist.begin(), dist.end(), -1);
        queue<int> q;
        dist[s] = 0;
        q.push(s);
        
        long long totalDist = 0;
        int reachable = 0;
        
        while (!q.empty()) {
            int u = q.front(); q.pop();
            const int* begin = g.neighbors_begin(u);
            const int* end   = g.neighbors_end(u);
            for (const int* it = begin; it != end; ++it) {
                int v = *it;
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    totalDist += dist[v];
                    reachable++;
                    q.push(v);
                }
            }
        }
        
        if (totalDist > 0) {
            cc[s] = (double)reachable / (double)totalDist;
        }
    }
    
    return cc;
}

// ============================================================
// Method 2: Multi-Source BFS using 64-bit bitmasks
//
// NOTE ON TERMINOLOGY:
//   This is BIT-PARALLELISM (a data structure trick), NOT SIMD.
//   - We pack 64 BFS traversals into a single uint64_t word
//   - Bitwise OR propagates frontiers for all 64 sources at once
//   - popcount counts how many sources discovered a given node
//   - This uses STANDARD scalar CPU instructions (OR, AND, POPCNT)
//   - There are NO SIMD intrinsics (no SSE, AVX, etc.)
//
// Why it's faster:
//   - Reading adjacency data from memory is the bottleneck
//   - Naive BFS reads the graph n times (once per source)
//   - MS-BFS reads the graph ⌈n/64⌉ times (64 sources per pass)
//   - This is a ~64× reduction in memory traffic
//   - Each "step" does 64 BFS frontier propagations with one OR
//
// The symmetry trick for standard closeness centrality:
//   For undirected graphs, d(u,v) = d(v,u).
//   When source 's' discovers destination 'v' at level 'l':
//     sumDist[s] += l         (row-wise: accumulate for source)
//     sumDist[v] += l         (column-wise: accumulate for destination)
//   Both are equivalent! We use column-wise because:
//     sumDist[v] += popcount(frontier_at_v) * level
//   This is a single operation per node per level, regardless of
//   how many sources are active.
// ============================================================
vector<double> multiSourceBFS_CC(const CSRGraph& g) {
    int n = g.n;
    vector<double> sumDist(n, 0.0);
    
    const int BATCH = 64; // bits in uint64_t
    
    // Temporary arrays — allocated once, reused per batch
    vector<uint64_t> frontier(n);
    vector<uint64_t> visited(n);
    vector<uint64_t> nextFrontier(n);
    
    for (int batchStart = 0; batchStart < n; batchStart += BATCH) {
        int batchEnd = min(batchStart + BATCH, n);
        int batchSize = batchEnd - batchStart;
        
        // Reset
        fill(frontier.begin(), frontier.end(), 0ULL);
        fill(visited.begin(), visited.end(), 0ULL);
        
        // Initialize: source i (within batch) has bit i set
        for (int i = 0; i < batchSize; i++) {
            int src = batchStart + i;
            frontier[src] = (1ULL << i);
            visited[src]  = (1ULL << i);
        }
        
        int level = 0;
        bool anyActive = true;
        
        while (anyActive) {
            level++;
            fill(nextFrontier.begin(), nextFrontier.end(), 0ULL);
            anyActive = false;
            
            // Frontier expansion: for each node u in the frontier,
            // propagate its bits to all neighbors
            for (int u = 0; u < n; u++) {
                uint64_t f = frontier[u];
                if (f == 0ULL) continue; // skip inactive nodes
                
                const int* begin = g.neighbors_begin(u);
                const int* end   = g.neighbors_end(u);
                for (const int* it = begin; it != end; ++it) {
                    nextFrontier[*it] |= f; // bitwise OR
                }
            }
            
            // Filter out visited nodes and accumulate distances
            for (int v = 0; v < n; v++) {
                nextFrontier[v] &= ~visited[v]; // remove already-visited
                if (nextFrontier[v] != 0ULL) {
                    anyActive = true;
                    visited[v] |= nextFrontier[v];
                    
                    // Symmetry trick:
                    // popcount = number of sources that just discovered v
                    // Each of those sources is at distance 'level' from v
                    // So sumDist[v] += popcount * level
                    int cnt = __builtin_popcountll(nextFrontier[v]);
                    sumDist[v] += (double)cnt * level;
                }
            }
            
            swap(frontier, nextFrontier);
        }
    }
    
    // Convert sum-of-distances to closeness centrality
    // CC(v) = (n - 1) / sumDist(v)
    vector<double> cc(n, 0.0);
    for (int v = 0; v < n; v++) {
        if (sumDist[v] > 0.0) {
            cc[v] = (double)(n - 1) / sumDist[v];
        }
    }
    
    return cc;
}

// ============================================================
// BCC Decomposition (Tarjan's Algorithm)
// ============================================================
struct BCCDecomposition {
    int n;
    vector<vector<int>> bccs;       // bccs[i] = list of nodes in BCC i
    vector<bool> isArticulation;     // isArticulation[v]
    vector<vector<int>> nodeBCCs;    // nodeBCCs[v] = BCCs containing v
    
    // Internal
    vector<int> disc, low, par;
    int timer_;
    stack<pair<int,int>> edgeStack;
    
    void extractBCC(int u, int v) {
        vector<int> bcc;
        while (!edgeStack.empty()) {
            auto [eu, ev] = edgeStack.top();
            edgeStack.pop();
            bcc.push_back(eu);
            bcc.push_back(ev);
            if ((eu == u && ev == v) || (eu == v && ev == u)) break;
        }
        sort(bcc.begin(), bcc.end());
        bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
        if (!bcc.empty()) {
            bccs.push_back(bcc);
        }
    }
    
    void dfs(const CSRGraph& g, int u) {
        disc[u] = low[u] = timer_++;
        int children = 0;
        
        const int* begin = g.neighbors_begin(u);
        const int* end   = g.neighbors_end(u);
        for (const int* it = begin; it != end; ++it) {
            int v = *it;
            if (disc[v] == -1) {
                children++;
                par[v] = u;
                edgeStack.push({u, v});
                dfs(g, v);
                low[u] = min(low[u], low[v]);
                
                // Articulation point check
                if ((par[u] == -1 && children > 1) ||
                    (par[u] != -1 && low[v] >= disc[u])) {
                    isArticulation[u] = true;
                    extractBCC(u, v);
                }
            } else if (v != par[u] && disc[v] < disc[u]) {
                edgeStack.push({u, v});
                low[u] = min(low[u], disc[v]);
            }
        }
    }
    
    BCCDecomposition(const CSRGraph& g) : n(g.n), isArticulation(n, false),
                                           nodeBCCs(n), disc(n, -1),
                                           low(n, 0), par(n, -1), timer_(0) {
        for (int i = 0; i < n; i++) {
            if (disc[i] == -1) {
                dfs(g, i);
                // Remaining edges form last BCC of this component
                if (!edgeStack.empty()) {
                    vector<int> bcc;
                    while (!edgeStack.empty()) {
                        auto [u, v] = edgeStack.top();
                        edgeStack.pop();
                        bcc.push_back(u);
                        bcc.push_back(v);
                    }
                    sort(bcc.begin(), bcc.end());
                    bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
                    bccs.push_back(bcc);
                }
            }
        }
        
        // Build node → BCC mapping
        for (int i = 0; i < (int)bccs.size(); i++) {
            for (int v : bccs[i]) {
                nodeBCCs[v].push_back(i);
            }
        }
    }
    
    void printStats() const {
        int numArt = count(isArticulation.begin(), isArticulation.end(), true);
        int maxBCC = 0, minBCC = n;
        long long totalBCCnodes = 0;
        for (const auto& bcc : bccs) {
            maxBCC = max(maxBCC, (int)bcc.size());
            minBCC = min(minBCC, (int)bcc.size());
            totalBCCnodes += bcc.size();
        }
        cout << "  BCCs: " << bccs.size() 
             << " | Articulation points: " << numArt
             << " | Largest BCC: " << maxBCC 
             << " | Smallest BCC: " << minBCC << endl;
    }
};

// ============================================================
// Chain Detection (degree-2 non-articulation nodes)
// ============================================================
int countChainNodes(const CSRGraph& g, const BCCDecomposition& bcc) {
    int count = 0;
    for (int v = 0; v < g.n; v++) {
        if (g.degree(v) == 2 && !bcc.isArticulation[v]) {
            count++;
        }
    }
    return count;
}

// ============================================================
// R3 Redundant Node Detection (degree-3 in a K4)
// ============================================================
int countR3Nodes(const CSRGraph& g) {
    int count = 0;
    for (int v = 0; v < g.n; v++) {
        if (g.degree(v) != 3) continue;
        
        const int* nb = g.neighbors_begin(v);
        int a = nb[0], b = nb[1], c = nb[2];
        
        // Check if a-b, b-c, a-c edges exist
        auto hasEdge = [&](int u, int w) -> bool {
            const int* begin = g.neighbors_begin(u);
            const int* end   = g.neighbors_end(u);
            for (const int* it = begin; it != end; ++it)
                if (*it == w) return true;
            return false;
        };
        
        if (hasEdge(a, b) && hasEdge(b, c) && hasEdge(a, c)) {
            count++;
        }
    }
    return count;
}

// ============================================================
// Validation: compare two CC vectors
// ============================================================
bool validateResults(const vector<double>& cc1, const vector<double>& cc2,
                     const string& name1, const string& name2) {
    double maxErr = 0.0;
    int worstNode = -1;
    for (int i = 0; i < (int)cc1.size(); i++) {
        double err = abs(cc1[i] - cc2[i]);
        if (err > maxErr) {
            maxErr = err;
            worstNode = i;
        }
    }
    cout << "  Max error (" << name1 << " vs " << name2 << "): " 
         << scientific << setprecision(2) << maxErr;
    if (maxErr < 1e-9) {
        cout << "  ✓ MATCH" << endl;
        return true;
    } else {
        cout << "  ✗ MISMATCH at node " << worstNode << endl;
        return false;
    }
}

// ============================================================
// main
// ============================================================
int main(int argc, char* argv[]) {
    CSRGraph g;
    
    if (argc > 1) {
        cout << "Reading graph from: " << argv[1] << endl;
        g = CSRGraph::readFromFile(argv[1]);
    } else {
        int n = 2000, m = 8000;
        cout << "No input file. Generating random graph: n=" << n << ", m=" << m << endl;
        g = CSRGraph::generateRandom(n, m);
    }
    
    cout << "Graph: " << g.n << " nodes, " << g.m << " edges\n" << endl;
    
    // --- Structural Analysis ---
    cout << "=== Structural Analysis ===" << endl;
    {
        Timer t;
        BCCDecomposition bcc(g);
        double ms = t.elapsed_ms();
        bcc.printStats();
        int chainNodes = countChainNodes(g, bcc);
        int r3Nodes = countR3Nodes(g);
        int removable = chainNodes + r3Nodes;
        cout << "  Chain nodes: " << chainNodes 
             << " | R3 nodes: " << r3Nodes
             << " | Total removable: " << removable 
             << " (" << fixed << setprecision(1) 
             << 100.0 * removable / g.n << "%)" << endl;
        cout << "  BCC decomposition time: " << fixed << setprecision(1) << ms << " ms\n" << endl;
    }
    
    // --- Method 1: Naive Sequential BFS ---
    cout << "=== Method 1: Naive Sequential BFS ===" << endl;
    vector<double> cc_naive;
    double naive_ms;
    {
        Timer t;
        cc_naive = naiveSequentialCC(g);
        naive_ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << naive_ms << " ms" << endl;
    }
    
    // --- Method 2: Multi-Source BFS (bit-parallelism, batch=64) ---
    cout << "\n=== Method 2: Multi-Source BFS (bit-parallelism, batch=64) ===" << endl;
    vector<double> cc_msbfs;
    double msbfs_ms;
    {
        Timer t;
        cc_msbfs = multiSourceBFS_CC(g);
        msbfs_ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << msbfs_ms << " ms" << endl;
    }
    
    // --- Validation ---
    cout << "\n=== Validation ===" << endl;
    validateResults(cc_naive, cc_msbfs, "Naive", "MS-BFS");
    
    // --- Speedup ---
    cout << "\n=== Performance Summary ===" << endl;
    cout << "  Naive BFS:       " << fixed << setprecision(1) << naive_ms << " ms" << endl;
    cout << "  MS-BFS (b=64):   " << fixed << setprecision(1) << msbfs_ms << " ms";
    if (naive_ms > 0 && msbfs_ms > 0)
        cout << "  (speedup: " << setprecision(2) << naive_ms / msbfs_ms << "x)";
    cout << endl;
    
    // --- Top-10 nodes ---
    cout << "\n=== Top 10 Nodes by Closeness Centrality ===" << endl;
    vector<int> idx(g.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc_naive[a] > cc_naive[b]; });
    for (int i = 0; i < min(10, g.n); i++) {
        cout << "  Node " << setw(5) << idx[i] 
             << ": CC = " << fixed << setprecision(6) << cc_naive[idx[i]] << endl;
    }
    
    return 0;
}
