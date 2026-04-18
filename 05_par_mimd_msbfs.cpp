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
    if (argc > 2) {
        g = Graph::readFromFile(argv[1]);
        nThreads = atoi(argv[2]);
    }
    else {
        cout << "Usage: " << argv[0] << " <graph_file> <num_threads>" << endl;
        return 1;
    }
    omp_set_num_threads(nThreads);
    int n = g.n;

    cout << "Method: Parallel MIMD + Multi-Source BFS (bit-parallelism)" << endl;
    cout << "Threads: " << nThreads << " | Batch: 64 BFS per word" << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_par = parallel_msbfs_cc(g, nThreads);
    auto t1 = chrono::high_resolution_clock::now();
    double par_ms = chrono::duration<double, milli>(t1 - t0).count();


    cout << "Parallel time:   " << fixed << setprecision(1) << par_ms << " ms" << endl;
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc_par[a] > cc_par[b]; });
    cout << "Top 10:" << endl;
    for (int i = 0; i < min(10, n); i++)
        cout << "  Node " << idx[i] << ": " << fixed << setprecision(64) << cc_par[idx[i]] << endl;
    ofstream values("a/5.csv");
    for (int i = 0; i < n; i++)
        values << cc_par[i] << endl;
    ofstream csv("experiment.csv", ios::app);
    csv<<"5"<<","<<par_ms<<"\n";
    

    return 0;
}
