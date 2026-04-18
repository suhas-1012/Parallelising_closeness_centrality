#include <bits/stdc++.h>
using namespace std;

vector<pair<int,int>> edges;

void add_edge(int u, int v) {
    if (u != v)
        edges.push_back({u, v});
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: ./bi n m\n";
        return 0;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);

    if (n < 6) {
        cout << "Need n >= 6 for multiple BCCs\n";
        return 0;
    }

    edges.clear();

    int k = 4; // number of BCC blocks (you can tune this)
    int base_size = n / k;

    vector<vector<int>> blocks;
    int curr = 0;

    // Step 1: create blocks
    for (int i = 0; i < k; i++) {
        int sz = (i == k-1) ? (n - curr) : base_size;
        vector<int> block;

        for (int j = 0; j < sz; j++)
            block.push_back(curr++);

        blocks.push_back(block);
    }

    // Step 2: make each block a cycle (=> each is a BCC)
    for (auto &block : blocks) {
        int sz = block.size();
        for (int i = 0; i < sz; i++)
            add_edge(block[i], block[(i+1)%sz]);
    }

    // Step 3: connect blocks via articulation points
    for (int i = 0; i < k-1; i++) {
        int u = blocks[i].back();      // articulation
        int v = blocks[i+1].front();   // next block
        add_edge(u, v);
    }

    // Step 4: add extra edges ONLY inside blocks
    srand(time(0));

    while ((int)edges.size() < m) {
        int b = rand() % k; // pick a block
        auto &block = blocks[b];

        int u = block[rand() % block.size()];
        int v = block[rand() % block.size()];

        if (u != v)
            add_edge(u, v);
    }
    ofstream fout("biconnected_graph.csv");
    fout << n << " " << edges.size() << "\n";
    for (auto& e : edges)
        fout << e.first << " " << e.second << "\n";
    fout.close();

    return 0;
}