#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <omp.h>

using namespace std;

//graph
struct Graph {
    int n, m;
    vector<vector<int>> adj;
    Graph() : n(0), m(0) {}
    Graph(int n) : n(n), m(0), adj(n) {}
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
        m++;
    }
    bool hasEdge(int u, int v) const {
        for (int w : adj[u]) if (w == v) return true;
        return false;
    }
    static Graph generateRandom(int n, int seed) {
        Graph g(n);
        srand(seed);
        for (int u = 0; u < n; u++)
            for (int v = u + 1; v < n; v++)
                if (rand() % 100 < 3)
                    g.addEdge(u, v);
        // ensure connected
        vector<bool> visited(n, false);
        queue<int> q; q.push(0); visited[0] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : g.adj[u])
                if (!visited[v]) { visited[v] = true; q.push(v); }
        }
        for (int i = 1; i < n; i++) {
            if (!visited[i]) {
                g.addEdge(0, i);
                visited[i] = true;
            }
        }
        return g;
    }
};

//Bcc decomposition using tarzan
struct BCCDecomposition {
    int n, timer_cnt;
    vector<int>  disc, low, par;
    vector<bool> is_art;
    vector<vector<int>> bcc_nodes;
    vector<vector<int>> bcc_arts;

    BCCDecomposition(int n) : n(n), timer_cnt(0),
        disc(n, -1), low(n, -1), par(n, -1), is_art(n, false) {}

    void dfs(const Graph& g, int u, vector<pair<int,int>>& stk) {
        disc[u] = low[u] = timer_cnt++;
        int children = 0;
        for (int v : g.adj[u]) {
            if (disc[v] == -1) {
                children++;
                par[v] = u;
                stk.push_back({u, v});
                dfs(g, v, stk);
                low[u] = min(low[u], low[v]);
                if ((par[u] == -1 && children > 1) ||
                    (par[u] != -1 && low[v] >= disc[u])) {
                    is_art[u] = true;
                    extract_bcc(stk, u, v);
                }
            } else if (v != par[u] && disc[v] < disc[u]) {
                low[u] = min(low[u], disc[v]);
                stk.push_back({u, v});
            }
        }
    }

    void extract_bcc(vector<pair<int,int>>& stk, int u, int v) {
        set<int> nodes;
        while (!stk.empty()) {
            auto [a, b] = stk.back(); stk.pop_back();
            nodes.insert(a); nodes.insert(b);
            if (a == u && b == v) break;
        }
        if (!nodes.empty())
            bcc_nodes.push_back(vector<int>(nodes.begin(), nodes.end()));
    }

    void run(const Graph& g) {
        vector<pair<int,int>> stk;
        for (int i = 0; i < n; i++) {
            if (disc[i] == -1) {
                dfs(g, i, stk);
                if (!stk.empty()) {
                    set<int> nodes;
                    while (!stk.empty()) {
                        auto [a, b] = stk.back(); stk.pop_back();
                        nodes.insert(a); nodes.insert(b);
                    }
                    bcc_nodes.push_back(vector<int>(nodes.begin(), nodes.end()));
                }
            }
        }
        bcc_arts.resize(bcc_nodes.size());
        for (int b = 0; b < (int)bcc_nodes.size(); b++)
            for (int gv : bcc_nodes[b])
                if (is_art[gv]) bcc_arts[b].push_back(gv);
    }

    int num_bccs() const { return (int)bcc_nodes.size(); }
};

//block cut tree structure
struct BCT {
    vector<vector<int>> adj;
    vector<int> type;  //0=bcc 1=art
    vector<int> id;
    vector<int> bcc_to_bct;
    unordered_map<int,int> art_to_bct;

    int size() const { return (int)adj.size(); }

    void build(const BCCDecomposition& bcc) {
        int nb = bcc.num_bccs();
        bcc_to_bct.resize(nb);
        for (int b = 0; b < nb; b++) {
            bcc_to_bct[b] = (int)adj.size();
            adj.push_back({}); type.push_back(0); id.push_back(b);
        }
        for (int i = 0; i < bcc.n; i++) {
            if (bcc.is_art[i]) {
                art_to_bct[i] = (int)adj.size();
                adj.push_back({}); type.push_back(1); id.push_back(i);
            }
        }
        for (int b = 0; b < nb; b++) {
            int bct_b = bcc_to_bct[b];
            for (int ga : bcc.bcc_arts[b]) {
                int bct_a = art_to_bct.at(ga);
                adj[bct_b].push_back(bct_a);
                adj[bct_a].push_back(bct_b);
            }
        }
    }
};

//local bcc structure
struct LocalBCC {
    int K;
    vector<int> l2g;
    vector<int> g2l_vec;  // sized to graph n, -1 for non-members
    vector<vector<int>> ladj;

    //phase 2 results
    vector<long long> intra_sum;
    int num_arts;
    vector<int> art_local;
    vector<int> art_global;
    vector<vector<int>> dist_to_art;  // K × num_arts

    //for cross-BCC routing
    int home_count;                    //non-art nodes in this BCC
    vector<long long> art_home_sum;    //[a] = Σ dist(art_a, v) for HOME v
    set<int> art_local_set;            //for quick lookup
};

//bfs within a local bcc
static vector<int> bfs_local(const LocalBCC& lb, int src) {
    vector<int> dist(lb.K, -1);
    queue<int> q;
    dist[src] = 0; q.push(src);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : lb.ladj[u])
            if (dist[v] == -1) { dist[v] = dist[u] + 1; q.push(v); }
    }
    return dist;
}

//phase 2 bcc confined bitwise spmm
void phase2_bcc_spmm(const Graph& g, const BCCDecomposition& bcc,
                     vector<LocalBCC>& lccs) {

    int nb = bcc.num_bccs();
    lccs.resize(nb);

    #pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nb; b++) {

        LocalBCC& lb = lccs[b];
        lb.K = (int)bcc.bcc_nodes[b].size();
        lb.l2g = bcc.bcc_nodes[b];

        lb.g2l_vec.assign(g.n, -1);
        for (int i = 0; i < lb.K; i++)
            lb.g2l_vec[lb.l2g[i]] = i;

        lb.ladj.assign(lb.K, {});
        for (int i = 0; i < lb.K; i++) {
            int u = lb.l2g[i];
            for (int v : g.adj[u]) {
                int j = lb.g2l_vec[v];
                if (j >= 0) lb.ladj[i].push_back(j);
            }
        }

        lb.art_global = bcc.bcc_arts[b];
        lb.num_arts = lb.art_global.size();
        lb.art_local.resize(lb.num_arts);

        for (int i = 0; i < lb.num_arts; i++)
            lb.art_local[i] = lb.g2l_vec[lb.art_global[i]];

        lb.art_local_set = set<int>(lb.art_local.begin(), lb.art_local.end());
        lb.home_count = lb.K - lb.num_arts;

        lb.intra_sum.assign(lb.K, 0);

        for (int base = 0; base < lb.K; base += 64) {

            int batch = min(64, lb.K - base);

            vector<uint64_t> frontier(lb.K, 0);
            vector<uint64_t> visited(lb.K, 0);
            vector<uint64_t> nextF(lb.K, 0);

            for (int i = 0; i < batch; i++) {
                uint64_t bit = 1ULL << i;
                frontier[base + i] = bit;
                visited[base + i] = bit;
            }

            int level = 0;
            bool active = true;

            while (active) {

                level++;
                active = false;

                size_t frontier_count = 0;
                for (int i = 0; i < lb.K; i++)
                    if (frontier[i]) frontier_count++;

                bool use_pull = (frontier_count > lb.K / 10);

                fill(nextF.begin(), nextF.end(), 0ULL);

                if (!use_pull) {

                    //push parallel, no race via local buffers
                    vector<vector<uint64_t>> thread_local_next(omp_get_max_threads(),
                                                               vector<uint64_t>(lb.K, 0));

                    #pragma omp parallel
                    {
                        int tid = omp_get_thread_num();
                        auto& local = thread_local_next[tid];

                        #pragma omp for schedule(dynamic)
                        for (int u = 0; u < lb.K; u++) {
                            if (!frontier[u]) continue;
                            for (int v : lb.ladj[u])
                                local[v] |= frontier[u];
                        }
                    }

                    for (auto& local : thread_local_next)
                        for (int i = 0; i < lb.K; i++)
                            nextF[i] |= local[i];

                } else {

                    //pull safe parallel
                    #pragma omp parallel for schedule(dynamic)
                    for (int v = 0; v < lb.K; v++) {
                        uint64_t bits = 0;
                        for (int u : lb.ladj[v])
                            bits |= frontier[u];
                        nextF[v] = bits;
                    }
                }

                #pragma omp parallel for reduction(||:active)
                for (int i = 0; i < lb.K; i++) {

                    uint64_t nw = nextF[i] & ~visited[i];

                    if (nw) {
                        visited[i] |= nw;
                        frontier[i] = nw;
                        lb.intra_sum[i] += (long long)__builtin_popcountll(nw) * level;
                        active = true;
                    } else {
                        frontier[i] = 0;
                    }
                }
            }
        }

        lb.dist_to_art.assign(lb.K, vector<int>(lb.num_arts, -1));

        for (int a = 0; a < lb.num_arts; a++) {
            vector<int> dist(lb.K, -1);
            queue<int> q;

            int s = lb.art_local[a];
            dist[s] = 0;
            q.push(s);

            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int v : lb.ladj[u]) {
                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        q.push(v);
                    }
                }
            }

            for (int i = 0; i < lb.K; i++)
                lb.dist_to_art[i][a] = dist[i];
        }

        lb.art_home_sum.assign(lb.num_arts, 0);

        for (int a = 0; a < lb.num_arts; a++) {
            for (int v = 0; v < lb.K; v++) {
                if (lb.art_local_set.count(v)) continue;
                if (lb.dist_to_art[v][a] >= 0)
                    lb.art_home_sum[a] += lb.dist_to_art[v][a];
            }
        }
    }
}

//phase 3 cross bcc dfs on bct

struct CrossBCCResult {
    long long ext_cnt;
    long long ext_dist;
};

CrossBCCResult dfs_cross(int bct_art, int from_bct_bcc,
                         const BCT& bct, const vector<LocalBCC>& lccs) {
    long long cnt = 0, dist = 0;
    int art_global = bct.id[bct_art];

    for (int bct_bcc : bct.adj[bct_art]) {
        if (bct_bcc == from_bct_bcc) continue;  // don't go back

        int b_idx = bct.id[bct_bcc];
        const LocalBCC& lb = lccs[b_idx];
        int entry_local = lb.g2l_vec[art_global];

        int entry_art_idx = -1;
        for (int a = 0; a < lb.num_arts; a++)
            if (lb.art_global[a] == art_global) { entry_art_idx = a; break; }

        cnt += lb.home_count;
        if (entry_art_idx >= 0)
            dist += lb.art_home_sum[entry_art_idx];

        for (int a = 0; a < lb.num_arts; a++) {
            if (lb.art_global[a] == art_global) continue;  // skip entry art

            int bridge = lb.dist_to_art[entry_local][a];
            if (bridge < 0) continue;

            cnt += 1;
            dist += bridge;

            int exit_bct_art = bct.art_to_bct.at(lb.art_global[a]);
            auto [rc, rd] = dfs_cross(exit_bct_art, bct_bcc, bct, lccs);
            cnt += rc;
            dist += rd + (long long)bridge * rc;
        }
    }
    return {cnt, dist};
}

//phase 4 assemble cc
void phase4_assemble(const Graph& g, const BCCDecomposition& bcc,
                     const vector<LocalBCC>& lccs,
                     const BCT& bct, vector<double>& cc) {
    int n = g.n;
    cc.assign(n, 0.0);

    vector<int> home_bcc(n, -1);
    vector<int> home_local(n, -1);
    for (int b = 0; b < (int)lccs.size(); b++)
        for (int li = 0; li < lccs[b].K; li++) {
            int gv = lccs[b].l2g[li];
            if (home_bcc[gv] == -1) {
                home_bcc[gv] = b;
                home_local[gv] = li;
            }
        }

    int nb = (int)lccs.size();
    vector<vector<CrossBCCResult>> bcc_ext(nb);

    for (int b = 0; b < nb; b++) {
        const LocalBCC& lb = lccs[b];
        bcc_ext[b].resize(lb.num_arts);
        int bct_bcc = bct.bcc_to_bct[b];
        for (int a = 0; a < lb.num_arts; a++) {
            int bct_art = bct.art_to_bct.at(lb.art_global[a]);
            bcc_ext[b][a] = dfs_cross(bct_art, bct_bcc, bct, lccs);
        }
    }

    #pragma omp parallel for schedule(dynamic, 32)
    for (int v = 0; v < n; v++) {
        int b = home_bcc[v];
        if (b < 0) continue;
        int lv = home_local[v];
        const LocalBCC& lb = lccs[b];

        long long D = lb.intra_sum[lv];

        for (int a = 0; a < lb.num_arts; a++) {
            int d_v_a = lb.dist_to_art[lv][a];
            if (d_v_a < 0) continue;

            D += (long long)d_v_a * bcc_ext[b][a].ext_cnt + bcc_ext[b][a].ext_dist;
        }

        long long ext_total = 0;
        for (int a = 0; a < lb.num_arts; a++)
            ext_total += bcc_ext[b][a].ext_cnt;
        long long reachable = min((long long)(lb.K) + ext_total - 1, (long long)(n - 1));

        if (D > 0)
            cc[v] = (double)reachable / (double)D;
    }
}

//dynamic edge insertion support
struct DynamicResult {
    double timeMs;
    int edgesAdded;
};

DynamicResult dynamicInsert(Graph& g, const vector<pair<int,int>>& edges,
                            vector<double>& cc) {
    auto t0 = chrono::high_resolution_clock::now();

    for (auto [u, v] : edges)
        g.addEdge(u, v);

    // Full redecompose (correct for insertions — BCCs may merge)
    BCCDecomposition bcc(g.n);
    bcc.run(g);

    BCT bct;
    bct.build(bcc);

    vector<LocalBCC> lccs;
    phase2_bcc_spmm(g, bcc, lccs);

    phase4_assemble(g, bcc, lccs, bct, cc);

    auto t1 = chrono::high_resolution_clock::now();
    return {chrono::duration<double, milli>(t1 - t0).count(), (int)edges.size()};
}

//baseline sequential BFS CC
void baseline_cc(const Graph& g, vector<double>& cc) {
    int n = g.n;
    cc.assign(n, 0.0);
    for (int s = 0; s < n; s++) {
        vector<int> dist(n, -1);
        queue<int> q;
        dist[s] = 0; q.push(s);
        long long td = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int nb : g.adj[u])
                if (dist[nb] == -1) { dist[nb] = dist[u]+1; q.push(nb); td += dist[nb]; }
        }
        if (td > 0) cc[s] = (double)(n - 1) / td;
    }
}


int main(int argc, char* argv[]) {
    int n       = (argc > 1) ? atoi(argv[1]) : 500;
    int seed    = (argc > 2) ? atoi(argv[2]) : 42;
    int threads = (argc > 3) ? atoi(argv[3]) : 4;
    omp_set_num_threads(threads);

    Graph g = Graph::generateRandom(n, seed);

    
    printf("==============================================\n");
    printf("BCC-Confined Vectorized Closeness Centrality\n");
    printf("Sariyuce 2014 + Shukla 2020 -- Novel MIMD\n");
    printf("==============================================\n");
    printf("n=%d  m=%d  threads=%d\n\n", n, g.m, threads);

    //phase 1A: BCC decomposition
    auto t0 = chrono::high_resolution_clock::now();
    BCCDecomposition bcc(n);
    bcc.run(g);
    auto t1 = chrono::high_resolution_clock::now();
    int num_arts = 0;
    for (int i = 0; i < n; i++) if (bcc.is_art[i]) num_arts++;
    printf("[Phase 1A] BCC decomposition: %.2f ms  BCCs=%d  arts=%d\n",
           chrono::duration<double,milli>(t1-t0).count(), bcc.num_bccs(), num_arts);

    //phase 1B: BCT
    BCT bct;
    bct.build(bcc);
    auto t2 = chrono::high_resolution_clock::now();
    printf("[Phase 1B] Block-Cut Tree:     %.2f ms  BCT nodes=%d\n",
           chrono::duration<double,milli>(t2-t1).count(), bct.size());

    //phase 2: BCC-confined SpMM
    vector<LocalBCC> lccs;
    phase2_bcc_spmm(g, bcc, lccs);
    auto t3 = chrono::high_resolution_clock::now();
    printf("[Phase 2]  BCC-SpMM (parallel): %.2f ms\n",
           chrono::duration<double,milli>(t3-t2).count());

    //phase 3+4: Cross-BCC DFS + assemble
    vector<double> my_cc;
    phase4_assemble(g, bcc, lccs, bct, my_cc);
    auto t4 = chrono::high_resolution_clock::now();
    printf("[Phase 3+4] Cross-BCC + Assemble: %.2f ms\n",
           chrono::duration<double,milli>(t4-t3).count());

    double total_ms = chrono::duration<double,milli>(t4-t0).count();

    //baseline
    auto tb0 = chrono::high_resolution_clock::now();
    vector<double> ref_cc;
    baseline_cc(g, ref_cc);
    auto tb1 = chrono::high_resolution_clock::now();
    double base_ms = chrono::duration<double,milli>(tb1-tb0).count();

    double max_err = 0.0, avg_err = 0.0;
    for (int v = 0; v < n; v++) {
        double e = fabs(my_cc[v] - ref_cc[v]);
        max_err = max(max_err, e);
        avg_err += e;
    }
    avg_err /= n;

    printf("\n==============================================\n");
    printf("RESULTS\n");
    printf("==============================================\n");
    printf("Baseline (seq BFS):     %.2f ms\n", base_ms);
    printf("BCC-SpMM (our method):  %.2f ms\n", total_ms);
    printf("Speedup:                %.2fx\n", base_ms / total_ms);
    printf("Max CC error:           %.2e\n", max_err);
    printf("Avg CC error:           %.2e\n", avg_err);
    printf("Validation: %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");

    printf("\n-- BCC stats --\n");
    long long total_work = 0;
    for (auto& lb : lccs) total_work += (long long)lb.K * lb.K;
    printf("Sum K^2 across BCCs:    %lld (vs N^2=%lld, ratio=%.3f)\n",
           total_work, (long long)n*n, (double)total_work / (n*n));

    printf("\n=== Dynamic Edge Insertions ===\n");
    srand(123);
    for (int b = 0; b < 3; b++) {
        vector<pair<int,int>> edges;
        for (int i = 0; i < 10; i++) {
            int u = rand() % n, v = rand() % n;
            while (v == u || g.hasEdge(u, v)) v = rand() % n;
            edges.push_back({u, v});
        }

        auto res = dynamicInsert(g, edges, my_cc);

        vector<double> check_cc;
        baseline_cc(g, check_cc);
        double err = 0;
        for (int i = 0; i < n; i++) err = max(err, fabs(my_cc[i] - check_cc[i]));

        printf("Batch %d: %d insertions | %.1f ms | err=%.1e %s\n",
               b + 1, res.edgesAdded, res.timeMs, err,
               err < 1e-6 ? "OK" : "ERR");
    }

    return 0;
}
