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
};

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
    // if (argc > 1) g = Graph::readFromFile(argv[1]);
    // else g = Graph::generateRandom(2000, 8000);
    // if (argc > 2) nThreads = atoi(argv[2]);
    if (argc > 2) {
        g = Graph::readFromFile(argv[1]);
        nThreads = atoi(argv[2]);
    }
    else {
        std::cout << "Usage: " << argv[0] << " <graph_file> <num_threads>" << std::endl;
        return 1;
    }
    cout << "Method: Parallel Simple MIMD (thread per source)" << endl;
    cout << "Threads: " << nThreads << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_par = parallel_naive_cc(g, nThreads);
    auto t1 = chrono::high_resolution_clock::now();
    double par_ms = chrono::duration<double, milli>(t1 - t0).count();

    cout << "Parallel time:   " << fixed << setprecision(1) << par_ms << " ms" << endl;
    vector<int> idx(g.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc_par[a] > cc_par[b]; });
    cout << "Top 10:" << endl;
    for (int i = 0; i < min(10, g.n); i++)
        cout << "  Node " << idx[i] << ": " << fixed << setprecision(64) << cc_par[idx[i]] << endl;
    ofstream values("a/4.csv");
    for (int i = 0; i < g.n; i++)
        values << cc_par[i] << "\n";
    ofstream csv("experiment.csv", ios::app);
    csv << "4," << par_ms << "\n";
    csv.close();
    return 0;
}
