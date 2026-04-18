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
};

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

            #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64) reduction(||:active)
            for (int v = 0; v < n; v++) {
                uint64_t bits = 0ULL;
                for (int u : g.adj[v]) {
                    bits |= frontier[u];
                }
                bits &= ~visited[v];
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
    if (argc > 1) {
        g = Graph::readFromFile(argv[1]);
    }
    else {
        cout << "Usage: " << argv[0] << " <graph_file>" << endl;
        return 1;
    }

    if (argc > 2) nThreads = atoi(argv[2]);
    omp_set_num_threads(nThreads);

    cout << "Method: Parallel Level-Synchronous Pull-Based BFS" << endl;
    cout << "Threads: " << nThreads << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_par = parallel_levelsync_cc(g, nThreads);
    auto t1 = chrono::high_resolution_clock::now();
    double par_ms = chrono::duration<double, milli>(t1 - t0).count();

    cout << "LevelSync time:  " << fixed << setprecision(1) << par_ms << " ms" << endl;
    vector<int> idx(g.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc_par[a] > cc_par[b]; });
    cout << "Top 10:" << endl;
    for (int i = 0; i < min(10, g.n); i++)
        cout << "  Node " << idx[i] << ": " << fixed << setprecision(64) << cc_par[idx[i]] << endl;
    ofstream values("a/6.csv");
    for (int i = 0; i < g.n; i++)
        values << cc_par[i] << "\n";
    ofstream csv("experiment.csv", ios::app);
    csv<<"6"<<","<<par_ms<<"\n";

    return 0;
}
