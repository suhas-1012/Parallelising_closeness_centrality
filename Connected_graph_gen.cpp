#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cstdlib>
using namespace std;

class Graph {
public:
    int n;
    vector<vector<int>> adj;

    Graph(int n) {
        this->n = n;
        adj.resize(n);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

Graph generateRandom(int n, int m) {
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

void saveCSV(Graph &g, int m) {
    ofstream file("connected_graph.csv");

    file << g.n << " " << m << "\n";

    set<pair<int,int>> printed;

    for (int u = 0; u < g.n; u++) {
        for (int v : g.adj[u]) {
            auto e = make_pair(min(u,v), max(u,v));
            if (!printed.count(e)) {
                file << e.first << " " << e.second << "\n";
                printed.insert(e);
            }
        }
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <n> <m>\n";
        return 1;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    if (m < n - 1) {
        cout << "Error: For connected graph, m must be >= n-1\n";
        return 1;
    }

    Graph g = generateRandom(n, m);
    saveCSV(g, m);

    return 0;
}