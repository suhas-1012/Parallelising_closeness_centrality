#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <set>
#include <stack>
#include <climits>

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

    bool hasEdge(int u, int v) const {
        for (int w : adj[u]) if (w == v) return true;
        return false;
    }
};

struct BCCDecomposition {
    int n, timer_;
    vector<int> disc, low, par;
    vector<bool> isArt;
    vector<vector<int>> bccs;
    stack<pair<int,int>> st;

    BCCDecomposition(const Graph& g) : n(g.n), timer_(0), disc(n,-1),
        low(n,0), par(n,-1), isArt(n, false) {
        for (int i = 0; i < n; i++)
            if (disc[i] == -1) { dfs(g, i); flushStack(); }
    }

    void flushStack() {
        if (st.empty()) return;
        vector<int> b;
        while (!st.empty()) {
            auto [u,v] = st.top(); st.pop();
            b.push_back(u); b.push_back(v);
        }
        sort(b.begin(), b.end());
        b.erase(unique(b.begin(), b.end()), b.end());
        bccs.push_back(b);
    }

    void extractBCC(int u, int v) {
        vector<int> b;
        while (!st.empty()) {
            auto [eu,ev] = st.top(); st.pop();
            b.push_back(eu); b.push_back(ev);
            if ((eu == u && ev == v) || (eu == v && ev == u)) break;
        }
        sort(b.begin(), b.end());
        b.erase(unique(b.begin(), b.end()), b.end());
        if (!b.empty()) bccs.push_back(b);
    }

    void dfs(const Graph& g, int u) {
        disc[u] = low[u] = timer_++;
        int children = 0;
        for (int v : g.adj[u]) {
            if (disc[v] == -1) {
                children++;
                par[v] = u;
                st.push({u, v});
                dfs(g, v);
                low[u] = min(low[u], low[v]);
                if ((par[u] == -1 && children > 1) || (par[u] != -1 && low[v] >= disc[u])) {
                    isArt[u] = true;
                    extractBCC(u, v);
                }
            } else if (v != par[u] && disc[v] < disc[u]) {
                st.push({u, v});
                low[u] = min(low[u], disc[v]);
            }
        }
    }
};

vector<int> findR3(const Graph& g) {
    vector<int> result;
    for (int v = 0; v < g.n; v++) {
        if ((int)g.adj[v].size() != 3) continue;
        int a = g.adj[v][0], b = g.adj[v][1], c = g.adj[v][2];
        if (g.hasEdge(a,b) && g.hasEdge(b,c) && g.hasEdge(a,c))
            result.push_back(v);
    }
    return result;
}

vector<int> findR4(const Graph& g) {
    vector<int> result;
    for (int v = 0; v < g.n; v++) {
        if ((int)g.adj[v].size() != 4) continue;
        int a = g.adj[v][0], b = g.adj[v][1], c = g.adj[v][2], d = g.adj[v][3];
        int perm[3][4] = {{a,b,c,d},{a,b,d,c},{a,c,b,d}};
        for (auto& p : perm) {
            if (g.hasEdge(p[0],p[1]) && g.hasEdge(p[1],p[2]) &&
                g.hasEdge(p[2],p[3]) && g.hasEdge(p[3],p[0]) &&
                !g.hasEdge(p[0],p[2]) && !g.hasEdge(p[1],p[3])) {
                result.push_back(v);
                break;
            }
        }
    }
    return result;
}

vector<double> bcc_reduced_cc(const Graph& g) {
    int n = g.n;
    BCCDecomposition bcc(g);
    vector<int> r3 = findR3(g);
    vector<int> r4 = findR4(g);

    vector<bool> redundant(n, false);
    for (int v : r3) redundant[v] = true;
    for (int v : r4) redundant[v] = true;

    int numRedundant = (int)(r3.size() + r4.size());
    int numArt = count(bcc.isArt.begin(), bcc.isArt.end(), true);

    cout << "BCCs: " << bcc.bccs.size() << endl;
    cout << "Articulation points: " << numArt << endl;
    cout << "R3 nodes: " << r3.size() << endl;
    cout << "R4 nodes: " << r4.size() << endl;
    cout << "Total redundant: " << numRedundant << " (" 
         << fixed << setprecision(1) << 100.0*numRedundant/n << "%)" << endl;

    vector<double> sumDist(n, 0.0);

    for (int s = 0; s < n; s++) {
        if (redundant[s]) continue;
        vector<int> dist(n, -1);
        queue<int> q;
        dist[s] = 0; q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : g.adj[u])
                if (dist[v] == -1) { dist[v] = dist[u]+1; q.push(v); }
        }
        for (int v = 0; v < n; v++) {
            if (v != s && dist[v] > 0) {
                sumDist[s] += dist[v];
                if (redundant[v]) sumDist[v] += dist[v];
            }
        }
    }

    for (int i = 0; i < (int)r3.size(); i++) {
        int ri = r3[i];
        for (int j = i+1; j < (int)r3.size(); j++) {
            int rj = r3[j];
            bool areNeighbors = g.hasEdge(ri, rj);
            if (areNeighbors) {
                sumDist[ri] += 1; sumDist[rj] += 1;
            } else {
                int minD = INT_MAX;
                for (int ni : g.adj[ri]) {
                    if (redundant[ni]) continue;
                    for (int nj : g.adj[rj]) {
                        if (redundant[nj]) continue;
                        if (ni == nj) { minD = min(minD, 2); continue; }
                        vector<int> dist(n, -1);
                        queue<int> bq;
                        dist[ni] = 0; bq.push(ni);
                        while (!bq.empty()) {
                            int u = bq.front(); bq.pop();
                            if (u == nj) break;
                            for (int w : g.adj[u])
                                if (dist[w] == -1) { dist[w] = dist[u]+1; bq.push(w); }
                        }
                        if (dist[nj] >= 0) minD = min(minD, dist[nj] + 2);
                    }
                }
                if (minD < INT_MAX) { sumDist[ri] += minD; sumDist[rj] += minD; }
            }
        }
        for (int j = 0; j < (int)r4.size(); j++) {
            int rj = r4[j];
            bool areNeighbors = g.hasEdge(ri, rj);
            if (areNeighbors) { sumDist[ri] += 1; sumDist[rj] += 1; }
            else {
                int minD = 2;
                for (int ni : g.adj[ri])
                    for (int nj : g.adj[rj])
                        if (!redundant[ni] && !redundant[nj] && ni == nj) minD = min(minD, 2);
                sumDist[ri] += minD; sumDist[rj] += minD;
            }
        }
    }

    for (int i = 0; i < (int)r4.size(); i++) {
        for (int j = i+1; j < (int)r4.size(); j++) {
            int ri = r4[i], rj = r4[j];
            if (g.hasEdge(ri, rj)) { sumDist[ri] += 1; sumDist[rj] += 1; }
            else {
                int minD = 2;
                for (int ni : g.adj[ri])
                    for (int nj : g.adj[rj])
                        if (!redundant[ni] && !redundant[nj] && ni == nj) minD = min(minD, 2);
                sumDist[ri] += minD; sumDist[rj] += minD;
            }
        }
    }

    vector<double> cc(n, 0.0);
    for (int v = 0; v < n; v++)
        if (sumDist[v] > 0) cc[v] = (double)(n-1) / sumDist[v];
    return cc;
}


int main(int argc, char* argv[]) {
    Graph g;
    if (argc > 1) {
        g = Graph::readFromFile(argv[1]);
    }
    else {
        cout << "Usage: " << argv[0] << " <graph_file>" << endl;
        return 1;
    }

    cout << "Method: Sequential BCC + R3/R4 Graph Reduction" << endl;
    cout << "Formula: CC(v) = (n-1) / sum_u d(v,u)" << endl;
    cout << "Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_bcc = bcc_reduced_cc(g);
    auto t1 = chrono::high_resolution_clock::now();
    double bcc_ms = chrono::duration<double, milli>(t1 - t0).count();

    cout << "BCC-reduced time: " << fixed << setprecision(1) << bcc_ms << " ms" << endl;
    cout << "BFS calls saved: " << (int)(findR3(g).size() + findR4(g).size()) << " / " << g.n << endl;

    vector<int> idx(g.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) { return cc_bcc[a] > cc_bcc[b]; });
    cout << "Top 10:" << endl;
    for (int i = 0; i < min(10, g.n); i++)
        cout << "  Node " << idx[i] << ": " << fixed << setprecision(64) << cc_bcc[idx[i]] << endl;
    ofstream values("a/3.csv");
    for (int i = 0; i < g.n; i++)
        values << cc_bcc[i] << endl;
    ofstream csv("experiment.csv", ios::app);
    csv<<"3"<<","<<bcc_ms<<"\n";
    csv.close();

    return 0;
}
