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
#include <cmath>
#include <omp.h>

using namespace std;

struct Graph {
    int n;
    vector<vector<int>> adj;
    Graph() : n(0) {}
    Graph(int n) : n(n), adj(n) {}

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void removeEdge(int u, int v) {
        adj[u].erase(find(adj[u].begin(), adj[u].end(), v));
        adj[v].erase(find(adj[v].begin(), adj[v].end(), u));
    }

    bool hasEdge(int u, int v) const {
        for (int w : adj[u]) if (w == v) return true;
        return false;
    }

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

struct DistanceOracle {
    int n;
    vector<vector<int>> dist;
    vector<double> sumDist;
    vector<double> cc;

    DistanceOracle(const Graph& g) : n(g.n), dist(n, vector<int>(n, -1)),
                                      sumDist(n, 0.0), cc(n, 0.0) {
        #pragma omp parallel
        {
            queue<int> q;  //thread-local queue

            #pragma omp for schedule(dynamic, 1)
            for (int s = 0; s < n; s++) {
                dist[s][s] = 0;
                q.push(s);
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    for (int v : g.adj[u]) {
                        if (dist[s][v] == -1) {
                            dist[s][v] = dist[s][u] + 1;
                            q.push(v);
                        }
                    }
                }
                double sd = 0;
                for (int v = 0; v < n; v++)
                    if (dist[s][v] > 0) sd += dist[s][v];
                sumDist[s] = sd;
                if (sd > 0) cc[s] = (double)(n-1) / sd;
            }
        }
    }
    void bfsFrom(const Graph& g, int s) {
        fill(dist[s].begin(), dist[s].end(), -1);
        dist[s][s] = 0;
        queue<int> q;
        q.push(s);
        sumDist[s] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : g.adj[u]) {
                if (dist[s][v] == -1) {
                    dist[s][v] = dist[s][u] + 1;
                    q.push(v);
                }
            }
        }
        for (int v = 0; v < n; v++)
            if (dist[s][v] > 0) sumDist[s] += dist[s][v];
        if (sumDist[s] > 0) cc[s] = (double)(n-1) / sumDist[s];
        else cc[s] = 0;
    }
};

struct EdgeUpdate {
    int u, v;
    bool isInsertion;
};

int dynamicUpdate(Graph& g, DistanceOracle& oracle, const vector<EdgeUpdate>& batch) {
    int n = g.n;
    int totalAffected = 0;

    for (auto& e : batch) {
        if (!e.isInsertion) continue;

        // affected sources using current (consistent) oracle
        vector<int> affectedList;
        for (int s = 0; s < n; s++) {
            if (oracle.dist[s][e.u] >= 0 && oracle.dist[s][e.v] >= 0 &&
                abs(oracle.dist[s][e.u] - oracle.dist[s][e.v]) > 1) {
                affectedList.push_back(s);
            }
        }

        g.addEdge(e.u, e.v);

        //e-BFS from affected sources (MIMD parallel)
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < (int)affectedList.size(); i++) {
            oracle.bfsFrom(g, affectedList[i]);
        }
        totalAffected += (int)affectedList.size();
    }


    return totalAffected;
}

int main(int argc, char* argv[]) {
    Graph g;
    if (argc > 1) g = Graph::readFromFile(argv[1]);
    else g = Graph::generateRandom(1000, 4000);

    int nThreads = omp_get_max_threads();
    cout << "Method: Dynamic CC (Shukla-style affected source filtering, MIMD parallel)" << endl;
    cout << "Threads: " << nThreads << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    DistanceOracle oracle(g);
    auto t1 = chrono::high_resolution_clock::now();
    double init_ms = chrono::duration<double, milli>(t1 - t0).count();
    cout << "Initial CC time: " << fixed << setprecision(1) << init_ms << " ms" << endl;

    srand(123);
    int numBatches = 5;
    int batchSize = 10;

    for (int b = 0; b < numBatches; b++) {
        vector<EdgeUpdate> batch;
        for (int i = 0; i < batchSize; i++) {
            int u = rand() % g.n;
            int v = rand() % g.n;
            while (v == u || g.hasEdge(u, v)) { v = rand() % g.n; }
            batch.push_back({u, v, true});
        }

        t0 = chrono::high_resolution_clock::now();
        int affected = dynamicUpdate(g, oracle, batch);
        t1 = chrono::high_resolution_clock::now();
        double update_ms = chrono::duration<double, milli>(t1 - t0).count();

        DistanceOracle fullRecomp(g);
        double maxErr = 0;
        for (int i = 0; i < g.n; i++)
            maxErr = max(maxErr, abs(oracle.cc[i] - fullRecomp.cc[i]));

        cout << "Batch " << b+1 << ": " << batchSize << " insertions, "
             << affected << "/" << g.n << " sources re-BFS'd, "
             << fixed << setprecision(1) << update_ms << " ms, "
             << "full-recomp=" << init_ms << " ms, "
             << "error=" << scientific << setprecision(2) << maxErr
             << (maxErr < 1e-9 ? " MATCH" : " MISMATCH") << endl;

        oracle = fullRecomp;
    }

    return 0;
}
