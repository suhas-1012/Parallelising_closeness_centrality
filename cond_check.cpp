#include <bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: ./check graph.csv\n";
        return 0;
    }

    ifstream fin(argv[1]);
    if (!fin) {
        cout << "Error opening file\n";
        return 0;
    }

    int n, m;
    fin >> n >> m;

    vector<vector<int>> adj(n);
    set<pair<int,int>> edges;

    bool hasDuplicate = false;
    bool hasSelfLoop = false;

    // Read edges
    for (int i = 0; i < m; i++) {
        int u, v;
        fin >> u >> v;

        if (u == v) {
            hasSelfLoop = true;
            continue;
        }

        // normalize (undirected)
        if (u > v) swap(u, v);

        if (edges.count({u, v})) {
            hasDuplicate = true;
        } else {
            edges.insert({u, v});
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }

    fin.close();

    // BFS for connectivity
    vector<bool> visited(n, false);
    queue<int> q;

    int start = 0;
    q.push(start);
    visited[start] = true;

    int visitedCount = 1;

    while (!q.empty()) {
        int u = q.front(); q.pop();

        for (int v : adj[u]) {
            if (!visited[v]) {
                visited[v] = true;
                q.push(v);
                visitedCount++;
            }
        }
    }

    bool isConnected = (visitedCount == n);

    // Output results
    cout << "===== Graph Check Report =====\n";

    cout << "Nodes: " << n << "\n";
    cout << "Edges (input): " << m << "\n";
    cout << "Unique edges: " << edges.size() << "\n\n";

    if (hasSelfLoop){
        cout << "Self-loops detected\n";
    } else {
        cout << "No self-loops\n";
    }
    if (hasDuplicate){
        cout << "Duplicate edges detected\n";
    } else {
        cout << "No duplicate edges\n";
    }
    if (isConnected){
        cout << "Graph is CONNECTED\n";
    } else {
        cout << "Graph is NOT connected\n";
    }

    return 0;
}