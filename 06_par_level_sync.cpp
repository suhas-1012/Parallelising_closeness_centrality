/*
 * 06_par_level_sync.cpp
 * ---------------------
 * Parallel Level-Synchronous Pull-Based Multi-Source BFS
 *
 * Algorithm: Multi-source BFS with 64-bit packing, but parallelized
 *            at each BFS level using a PULL (bottom-up) approach.
 *            Instead of frontier nodes pushing bits to their neighbors
 *            (which causes write conflicts → needs atomics),
 *            each unvisited target node PULLS bits from its neighbors.
 *            Since each thread writes to its own target node v, there are
 *            no write conflicts and no atomics/locks needed.
 *
 *            This is fine-grained level-synchronous MIMD: all threads
 *            cooperate on the SAME BFS level, with a barrier between levels.
 *            CC(v) = (n-1) / Σ d(v,u)
 *
 * Parallelism: MIMD — level-synchronous with pull direction (lock-free)
 * Complexity:  Time O(diameter × E / T),  Space O(V)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <set>
#include <cstdint>
#include <omp.h>

using namespace std;

struct Graph {
    int n;
    vector<vector<int>> adj;
    Graph() : n(0) {}
    Graph(int n) : n(n), adj(n) {}
    void addEdge(int u, int v) { adj[u].push_back(v); adj[v].push_back(u); }

    static Graph readFromFile(const string& f) {
        ifstream fin(f);
        int n, m; fin >> n >> m;
        Graph g(n);
        for (int i = 0; i < m; i++) { int u, v; fin >> u >> v; g.addEdge(u, v); }
        return g;
    }

    static Graph generateRandom(int n, int m) {
        Graph g(n);
        srand(42);
        for (int i = 1; i < n; i++) { int p = rand() % i; g.addEdge(i, p); }
        set<pair<int,int>> ex;
        for (int u = 0; u < n; u++)
            for (int v : g.adj[u]) ex.insert({min(u,v), max(u,v)});
        int added = n - 1;
        while (added < m) {
            int u = rand() % n, v = rand() % n;
            if (u == v) continue;
            auto e = make_pair(min(u,v), max(u,v));
            if (ex.count(e)) continue;
            ex.insert(e); g.addEdge(u, v); added++;
        }
        return g;
    }
};

vector<double> naive_cc(const Graph& g) {
    int n = g.n;
    vector<double> cc(n, 0.0);
    for (int s = 0; s < n; s++) {
        vector<int> dist(n, -1);
        queue<int> q;
        dist[s] = 0; q.push(s);
        long long td = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : g.adj[u])
                if (dist[v] == -1) { dist[v] = dist[u]+1; td += dist[v]; q.push(v); }
        }
        if (td > 0) cc[s] = (double)(n - 1) / td;
    }
    return cc;
}

/*
 * parallel_levelsync_cc:
 *   Multi-source BFS with 64-bit packing.
 *   Each BFS level is parallelized across threads using PULL direction:
 *     - Each thread processes a set of target nodes v
 *     - For each v, it reads (pulls) frontier bits from all neighbors
 *     - Each thread writes ONLY to its own v → no write conflicts, no atomics
 *     - Barrier between levels ensures correctness
 */
vector<double> parallel_levelsync_cc(const Graph& g, int nThreads) {
    int n = g.n;
    vector<double> sumDist(n, 0.0);
    const int BATCH = 64;
    vector<uint64_t> frontier(n), visited(n), nextF(n);

    for (int bs = 0; bs < n; bs += BATCH) {
        int be = min(bs + BATCH, n);
        int sz = be - bs;

        fill(frontier.begin(), frontier.end(), 0ULL);
        fill(visited.begin(), visited.end(), 0ULL);

        for (int i = 0; i < sz; i++) {
            frontier[bs + i] = (1ULL << i);
            visited[bs + i]  = (1ULL << i);
        }

        int level = 0;
        bool active = true;

        while (active) {
            level++;
            active = false;

            /*
             * PULL direction: each thread processes target nodes v.
             * For each v, pull frontier bits from neighbors (reads only).
             * Write to nextF[v] — each v is processed by exactly one thread.
             * No write conflicts → no atomics needed.
             */
            #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64) reduction(||:active)
            for (int v = 0; v < n; v++) {
                uint64_t bits = 0ULL;
                for (int u : g.adj[v]) {
                    bits |= frontier[u];   // PULL: read from neighbor's frontier
                }
                bits &= ~visited[v];       // mask already-visited sources
                nextF[v] = bits;

                if (bits != 0ULL) {
                    active = true;
                    visited[v] |= bits;
                    sumDist[v] += (double)__builtin_popcountll(bits) * level;
                }
            }

            swap(frontier, nextF);
        }
    }

    vector<double> cc(n, 0.0);
    for (int v = 0; v < n; v++)
        if (sumDist[v] > 0.0) cc[v] = (double)(n-1) / sumDist[v];
    return cc;
}

int main(int argc, char* argv[]) {
    Graph g;
    int nThreads = omp_get_max_threads();
    if (argc > 1) g = Graph::readFromFile(argv[1]);
    else g = Graph::generateRandom(2000, 8000);
    if (argc > 2) nThreads = atoi(argv[2]);

    cout << "Method: Parallel Level-Synchronous Pull-Based BFS" << endl;
    cout << "Threads: " << nThreads << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_seq = naive_cc(g);
    auto t1 = chrono::high_resolution_clock::now();
    double seq_ms = chrono::duration<double, milli>(t1 - t0).count();

    t0 = chrono::high_resolution_clock::now();
    auto cc_par = parallel_levelsync_cc(g, nThreads);
    t1 = chrono::high_resolution_clock::now();
    double par_ms = chrono::duration<double, milli>(t1 - t0).count();

    double maxErr = 0;
    for (int i = 0; i < g.n; i++) maxErr = max(maxErr, abs(cc_seq[i] - cc_par[i]));

    cout << "Sequential time: " << fixed << setprecision(1) << seq_ms << " ms" << endl;
    cout << "LevelSync time:  " << fixed << setprecision(1) << par_ms << " ms" << endl;
    cout << "Speedup: " << setprecision(2) << seq_ms / par_ms << "x" << endl;
    cout << "Max error: " << scientific << setprecision(2) << maxErr;
    cout << (maxErr < 1e-9 ? "  MATCH" : "  MISMATCH") << endl;

    return 0;
}
