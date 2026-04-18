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

vector<double> multisource_bfs_cc(const Graph& g) {
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
                    sumDist[v] += (double)__builtin_popcountll(nextF[v]) * level;
                }
            }
            swap(frontier, nextF);
        }
    }

    vector<double> cc(n, 0.0);
    for (int v = 0; v < n; v++)
        if (sumDist[v] > 0.0) cc[v] = (double)(n - 1) / sumDist[v];
    return cc;
}

int main(int argc, char* argv[]) {
    Graph g;
    if (argc > 1) {
        g = Graph::readFromFile(argv[1]);
    }
    else {
        std::cout << "Usage: " << argv[0] << " <graph_file>" << std::endl;
        return 1;
    }

    cout << "Method: Sequential Multi-Source BFS (bit-parallelism, batch=64)" << endl;
    cout << "Formula: CC(v) = (n-1) / sumDist[v], where sumDist[v] += popcount(frontier[v]) * level" << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_ms = multisource_bfs_cc(g);
    auto t1 = chrono::high_resolution_clock::now();
    double ms_ms = chrono::duration<double, milli>(t1 - t0).count();

    cout << "Time: " << fixed << setprecision(1) << ms_ms << " ms" << endl;
    vector<int> idx(g.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc_ms[a] > cc_ms[b]; });
    cout << "Top 10:" << endl;
    for (int i = 0; i < min(10, g.n); i++)
        cout << "  Node " << idx[i] << ": " << fixed << setprecision(64) << cc_ms[idx[i]] << endl;
    ofstream values("a/2.csv");
    for (int i = 0; i < g.n; i++)
        values << cc_ms[i] << endl;
    ofstream csv("experiment.csv", ios::app);
    csv<<"2"<<","<<ms_ms<<"\n";

    return 0;
}
