// parallel_cc.cpp — Phase 2: MIMD Parallel Closeness Centrality (OpenMP)
// ======================================================================
// MIMD = Multiple Instruction, Multiple Data
// Each thread independently runs its own BFS traversals on different
// source nodes. Threads share the graph (read-only) but write to
// separate output buffers → no synchronization overhead.
//
// Parallelism model:
//   ┌───────────────────────────────────────────────────────┐
//   │  MIMD (OpenMP threads)                                │
//   │  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐    │
//   │  │ Thread 0 │ │ Thread 1 │ │ Thread 2 │ │ Thread 3│   │
//   │  │ Sources  │ │ Sources  │ │ Sources  │ │ Sources │   │
//   │  │  0..k    │ │  k..2k   │ │ 2k..3k   │ │ 3k..n  │   │
//   │  │          │ │          │ │          │ │         │   │
//   │  │ MS-BFS   │ │ MS-BFS   │ │ MS-BFS   │ │ MS-BFS  │   │
//   │  │ batch=64 │ │ batch=64 │ │ batch=64 │ │ batch=64│   │
//   │  └─────────┘ └─────────┘ └─────────┘ └─────────┘    │
//   │  Each thread uses bit-parallelism (64 BFS per word)   │
//   │  to reduce memory traffic within its assigned range.  │
//   └───────────────────────────────────────────────────────┘
//
// No SIMD intrinsics are used. All operations are standard scalar
// CPU instructions (OR, AND, NOT, POPCNT).
//
// Compile: g++ -O2 -std=c++17 -fopenmp -o par_cc parallel_cc.cpp
// Usage:   ./par_cc [graph_file] [num_threads]
//          ./par_cc                     (random graph, auto threads)
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
#include <set>
#include <cstdint>
#include <omp.h>

using namespace std;

// ============================================================
// Timer
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
// CSR Graph
// ============================================================
struct CSRGraph {
    int n;
    long long m;
    vector<long long> rowPtr;
    vector<int> colIdx;
    
    CSRGraph() : n(0), m(0) {}
    
    CSRGraph(int n, const vector<vector<int>>& adj) : n(n) {
        rowPtr.resize(n + 1, 0);
        for (int u = 0; u < n; u++)
            rowPtr[u + 1] = adj[u].size();
        for (int u = 0; u < n; u++)
            rowPtr[u + 1] += rowPtr[u];
        m = rowPtr[n] / 2;
        colIdx.resize(rowPtr[n]);
        vector<long long> off(n, 0);
        for (int u = 0; u < n; u++)
            for (int v : adj[u])
                colIdx[rowPtr[u] + off[u]++] = v;
    }
    
    int degree(int u) const { return (int)(rowPtr[u + 1] - rowPtr[u]); }
    const int* neighbors_begin(int u) const { return &colIdx[rowPtr[u]]; }
    const int* neighbors_end(int u) const { return &colIdx[rowPtr[u + 1]]; }
    
    static CSRGraph readFromFile(const string& filename) {
        ifstream fin(filename);
        if (!fin) { cerr << "Error: Cannot open " << filename << endl; exit(1); }
        int n, m;
        fin >> n >> m;
        vector<vector<int>> adj(n);
        for (int i = 0; i < m; i++) {
            int u, v; fin >> u >> v;
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
        return CSRGraph(n, adj);
    }
    
    static CSRGraph generateRandom(int n, int m) {
        vector<vector<int>> adj(n);
        srand(42);
        for (int i = 1; i < n; i++) {
            int p = rand() % i;
            adj[i].push_back(p);
            adj[p].push_back(i);
        }
        set<pair<int,int>> existing;
        for (int u = 0; u < n; u++)
            for (int v : adj[u])
                existing.insert({min(u,v), max(u,v)});
        int added = n - 1;
        while (added < m) {
            int u = rand() % n, v = rand() % n;
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
// Sequential Naive BFS — baseline for correctness validation
// ============================================================
vector<double> naiveSequentialCC(const CSRGraph& g) {
    int n = g.n;
    vector<double> cc(n, 0.0);
    
    for (int s = 0; s < n; s++) {
        vector<int> dist(n, -1);
        queue<int> q;
        dist[s] = 0;
        q.push(s);
        long long totalDist = 0;
        int reachable = 0;
        
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (const int* it = g.neighbors_begin(u); it != g.neighbors_end(u); ++it) {
                if (dist[*it] == -1) {
                    dist[*it] = dist[u] + 1;
                    totalDist += dist[*it];
                    reachable++;
                    q.push(*it);
                }
            }
        }
        if (totalDist > 0) cc[s] = (double)reachable / totalDist;
    }
    return cc;
}

// ============================================================
// PARALLEL METHOD 1: Simple MIMD — one thread per BFS source
//   Each thread runs independent BFS traversals.
//   No shared writes, no synchronization, no atomics.
//   This is the simplest MIMD parallelization.
//
//   Parallelism: O(n) independent tasks → near-perfect scaling
// ============================================================
vector<double> parallelNaiveCC(const CSRGraph& g, int nThreads) {
    int n = g.n;
    vector<double> cc(n, 0.0);
    
    #pragma omp parallel num_threads(nThreads)
    {
        // Each thread has its own dist array — no sharing needed
        vector<int> dist(n);
        queue<int> q;
        
        #pragma omp for schedule(dynamic, 1)
        for (int s = 0; s < n; s++) {
            fill(dist.begin(), dist.end(), -1);
            dist[s] = 0;
            q.push(s);
            long long totalDist = 0;
            int reachable = 0;
            
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (const int* it = g.neighbors_begin(u); it != g.neighbors_end(u); ++it) {
                    if (dist[*it] == -1) {
                        dist[*it] = dist[u] + 1;
                        totalDist += dist[*it];
                        reachable++;
                        q.push(*it);
                    }
                }
            }
            if (totalDist > 0) cc[s] = (double)reachable / totalDist;
        }
    }
    return cc;
}

// ============================================================
// PARALLEL METHOD 2: MIMD + Bit-Parallelism
//   Combines two levels of parallelism:
//     Level 1 (MIMD):  OpenMP threads process different source batches
//     Level 2 (Bit):   Each thread does 64 BFS traversals per batch
//
//   Total parallelism = num_threads × 64
//
//   Each thread:
//     - Gets a range of source nodes
//     - Processes them in batches of 64 using bitmask MS-BFS
//     - Writes to thread-local sumDist array
//     - After all batches: reduce thread-local arrays into global
//
//   NO SIMD. NO GPU. Pure multi-core CPU parallelism.
// ============================================================
vector<double> parallelMSBFS_CC(const CSRGraph& g, int nThreads) {
    int n = g.n;
    vector<double> globalSumDist(n, 0.0);
    
    const int BATCH = 64;
    int numBatches = (n + BATCH - 1) / BATCH;
    
    #pragma omp parallel num_threads(nThreads)
    {
        // Thread-local storage — NO sharing between threads
        vector<double> localSumDist(n, 0.0);
        vector<uint64_t> frontier(n);
        vector<uint64_t> visited(n);
        vector<uint64_t> nextFrontier(n);
        
        // Each thread processes a subset of batches
        #pragma omp for schedule(dynamic, 1)
        for (int batch = 0; batch < numBatches; batch++) {
            int batchStart = batch * BATCH;
            int batchEnd = min(batchStart + BATCH, n);
            int batchSize = batchEnd - batchStart;
            
            // Reset arrays for this batch
            fill(frontier.begin(), frontier.end(), 0ULL);
            fill(visited.begin(), visited.end(), 0ULL);
            
            // Initialize sources
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
                
                // Frontier expansion (sequential within thread)
                for (int u = 0; u < n; u++) {
                    uint64_t f = frontier[u];
                    if (f == 0ULL) continue;
                    for (const int* it = g.neighbors_begin(u); it != g.neighbors_end(u); ++it) {
                        nextFrontier[*it] |= f;
                    }
                }
                
                // Update visited and accumulate distances
                for (int v = 0; v < n; v++) {
                    nextFrontier[v] &= ~visited[v];
                    if (nextFrontier[v] != 0ULL) {
                        anyActive = true;
                        visited[v] |= nextFrontier[v];
                        int cnt = __builtin_popcountll(nextFrontier[v]);
                        localSumDist[v] += (double)cnt * level;
                    }
                }
                
                swap(frontier, nextFrontier);
            }
        }
        
        // Reduce: merge thread-local results into global
        // This is a CRITICAL SECTION but happens only once at the end
        #pragma omp critical
        {
            for (int v = 0; v < n; v++) {
                globalSumDist[v] += localSumDist[v];
            }
        }
    }
    
    // Convert to CC
    vector<double> cc(n, 0.0);
    for (int v = 0; v < n; v++) {
        if (globalSumDist[v] > 0.0)
            cc[v] = (double)(n - 1) / globalSumDist[v];
    }
    return cc;
}

// ============================================================
// PARALLEL METHOD 3: MIMD + Bit-Parallelism + Better Reduction
//   Same as Method 2, but uses omp reduction-style merging
//   instead of a critical section. This avoids the serial
//   bottleneck at the end.
// ============================================================
vector<double> parallelMSBFS_CC_v2(const CSRGraph& g, int nThreads) {
    int n = g.n;
    
    const int BATCH = 64;
    int numBatches = (n + BATCH - 1) / BATCH;
    
    // Allocate per-thread sumDist arrays
    vector<vector<double>> threadSumDist(nThreads, vector<double>(n, 0.0));
    
    #pragma omp parallel num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
        vector<uint64_t> frontier(n);
        vector<uint64_t> visited(n);
        vector<uint64_t> nextFrontier(n);
        
        #pragma omp for schedule(dynamic, 1)
        for (int batch = 0; batch < numBatches; batch++) {
            int batchStart = batch * BATCH;
            int batchEnd = min(batchStart + BATCH, n);
            int batchSize = batchEnd - batchStart;
            
            fill(frontier.begin(), frontier.end(), 0ULL);
            fill(visited.begin(), visited.end(), 0ULL);
            
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
                
                for (int u = 0; u < n; u++) {
                    uint64_t f = frontier[u];
                    if (f == 0ULL) continue;
                    for (const int* it = g.neighbors_begin(u); it != g.neighbors_end(u); ++it) {
                        nextFrontier[*it] |= f;
                    }
                }
                
                for (int v = 0; v < n; v++) {
                    nextFrontier[v] &= ~visited[v];
                    if (nextFrontier[v] != 0ULL) {
                        anyActive = true;
                        visited[v] |= nextFrontier[v];
                        int cnt = __builtin_popcountll(nextFrontier[v]);
                        threadSumDist[tid][v] += (double)cnt * level;
                    }
                }
                
                swap(frontier, nextFrontier);
            }
        }
    }
    
    // Parallel reduction of per-thread arrays
    vector<double> cc(n, 0.0);
    #pragma omp parallel for num_threads(nThreads)
    for (int v = 0; v < n; v++) {
        double totalSumDist = 0.0;
        for (int t = 0; t < nThreads; t++) {
            totalSumDist += threadSumDist[t][v];
        }
        if (totalSumDist > 0.0)
            cc[v] = (double)(n - 1) / totalSumDist;
    }
    return cc;
}

// ============================================================
// PARALLEL METHOD 4: MIMD + Level-Synchronous Parallel BFS
//   This parallelizes WITHIN each BFS level (intra-BFS parallelism).
//   Multiple threads cooperate on expanding a single BFS frontier.
//   
//   Why: For large-diameter graphs, the MS-BFS approach has many
//   levels, each of which is sequential. This method makes each
//   level parallel across threads.
//
//   Uses: Two-phase approach
//     Phase A: Parallel frontier expansion (threads write to nextFrontier)
//     Phase B: Parallel visited update + distance accumulation
//   
//   Synchronization: Implicit barrier at end of each #pragma omp for
// ============================================================
vector<double> parallelLevelSyncCC(const CSRGraph& g, int nThreads) {
    int n = g.n;
    vector<double> sumDist(n, 0.0);
    
    const int BATCH = 64;
    
    // For this method, we process batches sequentially but
    // parallelize each BFS level within the batch
    vector<uint64_t> frontier(n);
    vector<uint64_t> visited(n);
    vector<uint64_t> nextFrontier(n);
    
    for (int batchStart = 0; batchStart < n; batchStart += BATCH) {
        int batchEnd = min(batchStart + BATCH, n);
        int batchSize = batchEnd - batchStart;
        
        fill(frontier.begin(), frontier.end(), 0ULL);
        fill(visited.begin(), visited.end(), 0ULL);
        
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
            
            // Phase A: Parallel frontier expansion
            // Each thread handles a chunk of nodes
            // WARNING: nextFrontier[v] may be written by multiple threads
            //   for different u's. We use atomic OR to handle this.
            #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64)
            for (int u = 0; u < n; u++) {
                uint64_t f = frontier[u];
                if (f == 0ULL) continue;
                for (const int* it = g.neighbors_begin(u); it != g.neighbors_end(u); ++it) {
                    int v = *it;
                    // Atomic OR — necessary because multiple threads
                    // may write to the same nextFrontier[v]
                    #pragma omp atomic
                    nextFrontier[v] |= f;
                }
            }
            
            // Phase B: Parallel update (no conflicts — each v is independent)
            #pragma omp parallel for num_threads(nThreads) schedule(static) reduction(||:anyActive)
            for (int v = 0; v < n; v++) {
                nextFrontier[v] &= ~visited[v];
                if (nextFrontier[v] != 0ULL) {
                    anyActive = true;
                    visited[v] |= nextFrontier[v];
                    int cnt = __builtin_popcountll(nextFrontier[v]);
                    sumDist[v] += (double)cnt * level;
                }
            }
            
            swap(frontier, nextFrontier);
        }
    }
    
    vector<double> cc(n, 0.0);
    for (int v = 0; v < n; v++) {
        if (sumDist[v] > 0.0)
            cc[v] = (double)(n - 1) / sumDist[v];
    }
    return cc;
}

// ============================================================
// Validation
// ============================================================
bool validateResults(const vector<double>& cc1, const vector<double>& cc2,
                     const string& name1, const string& name2) {
    double maxErr = 0.0;
    for (int i = 0; i < (int)cc1.size(); i++)
        maxErr = max(maxErr, abs(cc1[i] - cc2[i]));
    cout << "  Max error (" << name1 << " vs " << name2 << "): "
         << scientific << setprecision(2) << maxErr;
    if (maxErr < 1e-9) { cout << "  ✓ MATCH" << endl; return true; }
    else { cout << "  ✗ MISMATCH" << endl; return false; }
}

// ============================================================
// main
// ============================================================
int main(int argc, char* argv[]) {
    CSRGraph g;
    int nThreads = omp_get_max_threads();
    
    if (argc > 1) {
        cout << "Reading graph from: " << argv[1] << endl;
        g = CSRGraph::readFromFile(argv[1]);
    } else {
        int n = 2000, m = 8000;
        cout << "No input file. Generating random graph: n=" << n << ", m=" << m << endl;
        g = CSRGraph::generateRandom(n, m);
    }
    
    if (argc > 2) {
        nThreads = atoi(argv[2]);
    }
    
    cout << "Graph: " << g.n << " nodes, " << g.m << " edges" << endl;
    cout << "Threads: " << nThreads << "\n" << endl;
    
    // --- Sequential baseline ---
    cout << "=== Sequential Naive BFS (baseline) ===" << endl;
    vector<double> cc_seq;
    double seq_ms;
    {
        Timer t;
        cc_seq = naiveSequentialCC(g);
        seq_ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << seq_ms << " ms" << endl;
    }
    
    // --- Parallel Method 1: Simple MIMD ---
    cout << "\n=== Parallel Method 1: Simple MIMD (thread per source) ===" << endl;
    {
        Timer t;
        auto cc = parallelNaiveCC(g, nThreads);
        double ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << ms << " ms";
        cout << "  (speedup: " << setprecision(2) << seq_ms / ms << "x)" << endl;
        validateResults(cc_seq, cc, "Sequential", "Par-Naive");
    }
    
    // --- Parallel Method 2: MIMD + Bit-Parallelism (critical section) ---
    cout << "\n=== Parallel Method 2: MIMD + MS-BFS (critical section merge) ===" << endl;
    {
        Timer t;
        auto cc = parallelMSBFS_CC(g, nThreads);
        double ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << ms << " ms";
        cout << "  (speedup: " << setprecision(2) << seq_ms / ms << "x)" << endl;
        validateResults(cc_seq, cc, "Sequential", "Par-MSBFS-v1");
    }
    
    // --- Parallel Method 3: MIMD + Bit-Parallelism (parallel reduction) ---
    cout << "\n=== Parallel Method 3: MIMD + MS-BFS (parallel reduction) ===" << endl;
    {
        Timer t;
        auto cc = parallelMSBFS_CC_v2(g, nThreads);
        double ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << ms << " ms";
        cout << "  (speedup: " << setprecision(2) << seq_ms / ms << "x)" << endl;
        validateResults(cc_seq, cc, "Sequential", "Par-MSBFS-v2");
    }
    
    // --- Parallel Method 4: Level-Synchronous ---
    cout << "\n=== Parallel Method 4: Level-Synchronous Parallel BFS ===" << endl;
    {
        Timer t;
        auto cc = parallelLevelSyncCC(g, nThreads);
        double ms = t.elapsed_ms();
        cout << "  Time: " << fixed << setprecision(1) << ms << " ms";
        cout << "  (speedup: " << setprecision(2) << seq_ms / ms << "x)" << endl;
        validateResults(cc_seq, cc, "Sequential", "Par-LevelSync");
    }
    
    // --- Scaling analysis ---
    cout << "\n=== Thread Scaling Analysis (Method 3: MIMD + MS-BFS) ===" << endl;
    cout << "  " << setw(8) << "Threads" << setw(12) << "Time(ms)" 
         << setw(10) << "Speedup" << setw(12) << "Efficiency" << endl;
    cout << "  " << string(42, '-') << endl;
    
    int maxT = min(nThreads, 16);
    for (int t = 1; t <= maxT; t *= 2) {
        Timer timer;
        auto cc = parallelMSBFS_CC_v2(g, t);
        double ms = timer.elapsed_ms();
        double speedup = seq_ms / ms;
        double efficiency = speedup / t * 100.0;
        cout << "  " << setw(8) << t 
             << setw(12) << fixed << setprecision(1) << ms
             << setw(10) << setprecision(2) << speedup << "x"
             << setw(11) << setprecision(1) << efficiency << "%" << endl;
    }
    if (maxT != nThreads) {
        Timer timer;
        auto cc = parallelMSBFS_CC_v2(g, nThreads);
        double ms = timer.elapsed_ms();
        double speedup = seq_ms / ms;
        double efficiency = speedup / nThreads * 100.0;
        cout << "  " << setw(8) << nThreads
             << setw(12) << fixed << setprecision(1) << ms
             << setw(10) << setprecision(2) << speedup << "x"
             << setw(11) << setprecision(1) << efficiency << "%" << endl;
    }
    
    return 0;
}
