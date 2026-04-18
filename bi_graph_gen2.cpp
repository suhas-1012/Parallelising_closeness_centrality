#include <bits/stdc++.h>
using namespace std;

vector<pair<int,int>> edges;

void add_edge(int u, int v) {
    if (u != v)
        edges.emplace_back(u, v);
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

    edges.clear();

    int k = max(2, n / 2000); // number of BCC blocks (tunable)
    int curr = 0;

    vector<vector<int>> blocks;

    // Step 1: split nodes into blocks
    for (int i = 0; i < k; i++) {
        int remaining = n - curr;
        int sz = (i == k-1) ? remaining : max(3, remaining / (k - i));

        vector<int> block;
        for (int j = 0; j < sz; j++)
            block.push_back(curr++);

        blocks.push_back(block);
    }

    // Step 2: make each block a cycle (=> BCC)
    for (auto &block : blocks) {
        int sz = block.size();
        for (int i = 0; i < sz; i++)
            add_edge(block[i], block[(i+1)%sz]);
    }

    // Step 3: connect blocks (THIS ensures connectivity + articulation)
    for (int i = 0; i < (int)blocks.size() - 1; i++) {
        int u = blocks[i].back();      // articulation point
        int v = blocks[i+1].front();   // next block
        add_edge(u, v);
    }

    // Step 4: add extra edges ONLY inside blocks (preserve BCCs)
    srand(time(0));

    while ((int)edges.size() < m) {
        int b = rand() % blocks.size();
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