/*
 * 04_par_simple_mimd.cpp
 * ----------------------
 * Parallel MIMD Naive BFS — one thread per source node
 *
 * Algorithm: Each OpenMP thread independently picks a source node and runs
 *            a full BFS from that source.  No shared state between threads
 *            (each thread has its own dist[] and queue).  Pure MIMD.
 *            CC(v) = (n-1) / Σ d(v,u)
 *
 * Parallelism: MIMD — OpenMP parallel for with dynamic scheduling
 * Complexity:  Time O(V × (V + E) / T),  Space O(V × T)
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

vector<double> parallel_naive_cc(const Graph& g, int nThreads) {
    int n = g.n;
    vector<double> cc(n, 0.0);

    #pragma omp parallel num_threads(nThreads)
    {
        vector<int> dist(n);
        queue<int> q;

        #pragma omp for schedule(dynamic, 1)
        for (int s = 0; s < n; s++) {
            fill(dist.begin(), dist.end(), -1);
            dist[s] = 0; q.push(s);
            long long td = 0; int reach = 0;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int v : g.adj[u])
                    if (dist[v] == -1) { dist[v] = dist[u]+1; td += dist[v]; reach++; q.push(v); }
            }
            if (td > 0) cc[s] = (double)(n - 1) / td;
        }
    }
    return cc;
}

int main(int argc, char* argv[]) {
    Graph g;
    int nThreads = omp_get_max_threads();
    if (argc > 1) g = Graph::readFromFile(argv[1]);
    else g = Graph::generateRandom(2000, 8000);
    if (argc > 2) nThreads = atoi(argv[2]);

    cout << "Method: Parallel Simple MIMD (thread per source)" << endl;
    cout << "Threads: " << nThreads << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_seq = naive_cc(g);
    auto t1 = chrono::high_resolution_clock::now();
    double seq_ms = chrono::duration<double, milli>(t1 - t0).count();

    t0 = chrono::high_resolution_clock::now();
    auto cc_par = parallel_naive_cc(g, nThreads);
    t1 = chrono::high_resolution_clock::now();
    double par_ms = chrono::duration<double, milli>(t1 - t0).count();

    double maxErr = 0;
    for (int i = 0; i < g.n; i++) maxErr = max(maxErr, abs(cc_seq[i] - cc_par[i]));

    cout << "Sequential time: " << fixed << setprecision(1) << seq_ms << " ms" << endl;
    cout << "Parallel time:   " << fixed << setprecision(1) << par_ms << " ms" << endl;
    cout << "Speedup: " << setprecision(2) << seq_ms / par_ms << "x" << endl;
    cout << "Max error: " << scientific << setprecision(2) << maxErr;
    cout << (maxErr < 1e-9 ? "  MATCH" : "  MISMATCH") << endl;

    return 0;
}
