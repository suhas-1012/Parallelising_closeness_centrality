#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <set>

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

/*
 * closeness_centrality:
 *   For each source s = 0..n-1:
 *     BFS from s → compute dist[v] for all v
 *     totalDist = Σ dist[v]
 *     CC(s) = (n-1) / totalDist
 */
vector<double> closeness_centrality(const Graph& g) {
    int n = g.n;
    vector<double> cc(n, 0.0);

    for (int s = 0; s < n; s++) {
        vector<int> dist(n, -1);
        queue<int> q;
        dist[s] = 0;
        q.push(s);
        long long totalDist = 0;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : g.adj[u]) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    totalDist += dist[v];
                    q.push(v);
                }
            }
        }
        if (totalDist > 0) cc[s] = (double)(n - 1) / (double)totalDist;
    }
    return cc;
}

int main(int argc, char* argv[]) {
    Graph g;
    if (argc > 1) g = Graph::readFromFile(argv[1]);
    else g = Graph::generateRandom(2000, 8000);

    cout << "Method: Sequential Naive BFS" << endl;
    cout << "Formula: CC(v) = (n-1) / sum_u d(v,u)" << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc = closeness_centrality(g);
    auto t1 = chrono::high_resolution_clock::now();
    double ms = chrono::duration<double, milli>(t1 - t0).count();

    cout << "Time: " << fixed << setprecision(1) << ms << " ms" << endl;

    vector<int> idx(g.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc[a] > cc[b]; });
    cout << "Top 10:" << endl;
    for (int i = 0; i < min(10, g.n); i++)
        cout << "  Node " << idx[i] << ": " << fixed << setprecision(6) << cc[idx[i]] << endl;

    return 0;
}
