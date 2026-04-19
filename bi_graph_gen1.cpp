#include <bits/stdc++.h>
using namespace std;

set<pair<int,int>> edges;

// Add undirected edge (normalized)
void add_edge(int u, int v) {
    if (u == v) return;
    if (u > v) swap(u, v);
    edges.insert({u, v});
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: ./bi n m\n";
        return 0;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    if (n < 6) {
        cout << "Need n >= 6\n";
        return 0;
    }

    // Max possible edges in undirected graph
    long long max_edges = 1LL * n * (n - 1) / 2;
    if (m > max_edges) {
        cout << "Too many edges requested. Max = " << max_edges << "\n";
        return 0;
    }

    edges.clear();

    int k = max(2, n / 2000); // number of BCC blocks
    int curr = 0;

    vector<vector<int>> blocks;

    // Step 1: split nodes into blocks
    for (int i = 0; i < k; i++) {
        int remaining = n - curr;
        int sz = (i == k - 1) ? remaining : max(3, remaining / (k - i));

        vector<int> block;
        for (int j = 0; j < sz; j++)
            block.push_back(curr++);

        blocks.push_back(block);
    }

    // Step 2: make each block a cycle (=> BCC)
    for (auto &block : blocks) {
        int sz = block.size();
        for (int i = 0; i < sz; i++)
            add_edge(block[i], block[(i + 1) % sz]);
    }

    // Step 3: connect blocks (articulation points)
    for (int i = 0; i < (int)blocks.size() - 1; i++) {
        int u = blocks[i].back();
        int v = blocks[i + 1].front();
        add_edge(u, v);
    }

    // Step 4: add extra edges ONLY inside blocks
    srand(42);

    while ((int)edges.size() < m) {
        int b = rand() % blocks.size();
        auto &block = blocks[b];

        int u = block[rand() % block.size()];
        int v = block[rand() % block.size()];

        add_edge(u, v);
    }

    // Output
    ofstream fout("biconnected_graph.csv");
    fout << n << " " << edges.size() << "\n";

    for (auto &e : edges)
        fout << e.first << " " << e.second << "\n";

    fout.close();

    cout << "Graph generated: biconnected_graph.csv\n";
    cout << "Nodes: " << n << ", Edges: " << edges.size() << "\n";

    return 0;
}