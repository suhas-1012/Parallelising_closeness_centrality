/*
 * 05_par_mimd_msbfs.cpp
 * ---------------------
 * Parallel MIMD + Multi-Source BFS (bit-parallelism)
 *
 * Algorithm: Combines MIMD thread-level parallelism with 64-bit
 *            bit-packing.  Each thread independently processes a
 *            batch of 64 sources using bitwise OR-semiring SpMM.
 *            Threads do NOT share frontier state — pure MIMD.
 *            CC(v) = (n-1) / Σ d(v,u)
 *
 * Parallelism: MIMD — each thread owns a batch of 64 BFS sources
 * Complexity:  Time O(V × E / (64 × T)),  Space O(V × T)
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
        long long td = 0; int reach = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : g.adj[u])
                if (dist[v] == -1) { dist[v] = dist[u]+1; td += dist[v]; reach++; q.push(v); }
        }
        if (td > 0) cc[s] = (double)(n - 1) / td;
    }
    return cc;
}

vector<double> parallel_msbfs_cc(const Graph& g, int nThreads) {
    int n = g.n;
    const int BATCH = 64;
    int numBatches = (n + BATCH - 1) / BATCH;

    vector<vector<double>> threadSumDist(nThreads, vector<double>(n, 0.0));

    #pragma omp parallel num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
        vector<uint64_t> frontier(n), visited(n), nextF(n);

        #pragma omp for schedule(dynamic, 1)
        for (int batch = 0; batch < numBatches; batch++) {
            int bs = batch * BATCH;
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
                fill(nextF.begin(), nextF.end(), 0ULL);
                active = false;

                for (int u = 0; u < n; u++) {
                    if (frontier[u] == 0ULL) continue;
                    for (int v : g.adj[u])
                        nextF[v] |= frontier[u];
                }

                for (int v = 0; v < n; v++) {
                    nextF[v] &= ~visited[v];
                    if (nextF[v] != 0ULL) {
                        active = true;
                        visited[v] |= nextF[v];
                        threadSumDist[tid][v] += (double)__builtin_popcountll(nextF[v]) * level;
                    }
                }

                swap(frontier, nextF);
            }
        }
    }

    vector<double> cc(n, 0.0);
    #pragma omp parallel for num_threads(nThreads)
    for (int v = 0; v < n; v++) {
        double total = 0.0;
        for (int t = 0; t < nThreads; t++)
            total += threadSumDist[t][v];
        if (total > 0.0)
            cc[v] = (double)(n - 1) / total;
    }
    return cc;
}

int main(int argc, char* argv[]) {
    Graph g;
    int nThreads = omp_get_max_threads();
    if (argc > 1) g = Graph::readFromFile(argv[1]);
    else g = Graph::generateRandom(2000, 8000);
    if (argc > 2) nThreads = atoi(argv[2]);

    cout << "Method: Parallel MIMD + Multi-Source BFS (bit-parallelism)" << endl;
    cout << "Threads: " << nThreads << " | Batch: 64 BFS per word" << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_seq = naive_cc(g);
    auto t1 = chrono::high_resolution_clock::now();
    double seq_ms = chrono::duration<double, milli>(t1 - t0).count();

    t0 = chrono::high_resolution_clock::now();
    auto cc_par = parallel_msbfs_cc(g, nThreads);
    t1 = chrono::high_resolution_clock::now();
    double par_ms = chrono::duration<double, milli>(t1 - t0).count();

    double maxErr = 0;
    for (int i = 0; i < g.n; i++) maxErr = max(maxErr, abs(cc_seq[i] - cc_par[i]));

    cout << "Sequential time: " << fixed << setprecision(1) << seq_ms << " ms" << endl;
    cout << "Parallel time:   " << fixed << setprecision(1) << par_ms << " ms" << endl;
    cout << "Speedup: " << setprecision(2) << seq_ms / par_ms << "x" << endl;
    cout << "Max error: " << scientific << setprecision(2) << maxErr;
    cout << (maxErr < 1e-9 ? "  MATCH" : "  MISMATCH") << endl;

    cout << "\nScaling:" << endl;
    for (int t = 1; t <= nThreads; t *= 2) {
        t0 = chrono::high_resolution_clock::now();
        parallel_msbfs_cc(g, t);
        t1 = chrono::high_resolution_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();
        cout << "  " << t << " threads: " << fixed << setprecision(1) << ms
             << " ms  speedup=" << setprecision(2) << seq_ms/ms << "x" << endl;
    }

    return 0;
}
