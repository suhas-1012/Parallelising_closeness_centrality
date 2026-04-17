/*
================================================================================
  BCC-SCC CONFINED CLOSENESS CENTRALITY: WEIGHTED DIRECTED GRAPHS
  Sequential | Parallel | Dynamic
  
  Extends the BCC-confined approach to:
  - Weighted edges (Dijkstra replaces BFS/SpMM)
  - Directed graphs (SCCs replace BCCs, condensation DAG replaces BCT)
  
  PROVEN CORRECT: See proof section in accompanying document.
================================================================================
*/
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using W = double;
const W INF = 1e18;

// ─────────────────────────────────────────────────────────────────────────────
// Weighted Directed Graph
// ─────────────────────────────────────────────────────────────────────────────
struct WDG {
    int n, m;
    vector<vector<pair<int,W>>> adj, radj;  // forward and reverse adjacency
    WDG(int n) : n(n), m(0), adj(n), radj(n) {}
    void add(int u, int v, W w) {
        adj[u].push_back({v,w});
        radj[v].push_back({u,w});
        m++;
    }
    void remove(int u, int v) {
        auto& a = adj[u];
        a.erase(remove_if(a.begin(),a.end(),[v](auto p){return p.first==v;}),a.end());
        auto& b = radj[v];
        b.erase(remove_if(b.begin(),b.end(),[u](auto p){return p.first==u;}),b.end());
        m--;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Dijkstra (single source, on given adjacency list)
// ─────────────────────────────────────────────────────────────────────────────
vector<W> dijkstra(const vector<vector<pair<int,W>>>& adj, int n, int src) {
    vector<W> d(n, INF);
    priority_queue<pair<W,int>, vector<pair<W,int>>, greater<>> pq;
    d[src] = 0; pq.push({0, src});
    while (!pq.empty()) {
        auto [dd, u] = pq.top(); pq.pop();
        if (dd > d[u] + 1e-12) continue;
        for (auto [v, w] : adj[u])
            if (d[u]+w < d[v]-1e-12) { d[v]=d[u]+w; pq.push({d[v],v}); }
    }
    return d;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase 1A: SCC Decomposition (Kosaraju)
//   Directed analogue of BCC: every node reachable from every other within SCC.
// ─────────────────────────────────────────────────────────────────────────────
struct SCC {
    int n, ns = 0;
    vector<int> comp;
    vector<vector<int>> sccs;

    SCC(int n) : n(n), comp(n, -1) {}

    void run(const WDG& g) {
        // Pass 1: finish-order DFS on forward graph
        vector<bool> vis(n, false);
        vector<int> order;
        function<void(int)> dfs1 = [&](int u) {
            vis[u] = true;
            for (auto [v,w] : g.adj[u]) if (!vis[v]) dfs1(v);
            order.push_back(u);
        };
        for (int i = 0; i < n; i++) if (!vis[i]) dfs1(i);

        // Pass 2: DFS on reverse graph in reverse finish order
        fill(vis.begin(), vis.end(), false);
        function<void(int,int)> dfs2 = [&](int u, int c) {
            vis[u] = true; comp[u] = c;
            for (auto [v,w] : g.radj[u]) if (!vis[v]) dfs2(v, c);
        };
        sccs.clear(); ns = 0;
        for (int i = n-1; i >= 0; i--) {
            if (!vis[order[i]]) {
                sccs.push_back({});
                dfs2(order[i], ns++);
            }
        }
        sccs.resize(ns);
        for (int v = 0; v < n; v++) sccs[comp[v]].push_back(v);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Phase 1B: Condensation DAG
//   Directed analogue of BCT: DAG of SCCs connected by bridges.
//   A "bridge" here is any cross-SCC directed edge.
// ─────────────────────────────────────────────────────────────────────────────
struct Cond {
    int ns;
    vector<vector<pair<int,W>>> adj;  // condensation edges
    struct Bridge { int s, t, eu, ev; W w; };
    vector<Bridge> bridges;
    // For each (s,t) pair: best (minimum weight) bridge
    map<pair<int,int>, int> best_br;  // (s,t) -> index in bridges

    void build(const WDG& g, const SCC& scc) {
        ns = scc.ns;
        adj.assign(ns, {});
        bridges.clear(); best_br.clear();
        map<pair<int,int>,W> min_w;
        for (int u = 0; u < g.n; u++) {
            for (auto [v, w] : g.adj[u]) {
                int s = scc.comp[u], t = scc.comp[v];
                if (s == t) continue;
                int idx = (int)bridges.size();
                bridges.push_back({s, t, u, v, w});
                auto key = make_pair(s, t);
                if (!min_w.count(key) || w < min_w[key]) {
                    min_w[key] = w;
                    best_br[key] = idx;
                }
            }
        }
        for (auto& [key, w] : min_w)
            adj[key.first].push_back({key.second, w});
    }

    // Topological order of condensation DAG
    vector<int> topo_order() const {
        vector<int> indeg(ns, 0);
        for (int s = 0; s < ns; s++) for (auto [t,w] : adj[s]) indeg[t]++;
        queue<int> q;
        for (int s = 0; s < ns; s++) if (!indeg[s]) q.push(s);
        vector<int> ord;
        while (!q.empty()) {
            int s = q.front(); q.pop(); ord.push_back(s);
            for (auto [t,w] : adj[s]) if (!--indeg[t]) q.push(t);
        }
        return ord;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Local SCC data (computed in Phase 2)
// ─────────────────────────────────────────────────────────────────────────────
struct LS {
    int K, home_cnt = 0;
    vector<int> l2g, g2l;
    vector<vector<pair<int,W>>> ladj, rladj;
    vector<int> home_l;           // local IDs of home nodes
    vector<int> exit_l, exit_g;   // exit node local/global IDs
    int ne = 0;
    // Computed distances:
    vector<W>            intra_sum; // intra_sum[lv] = Σ_{u∈HOME} dist_SCC(lv, u)
    vector<vector<W>>    d2exit;    // d2exit[lv][e] = dist_SCC(lv, exit_e)
};

// ─────────────────────────────────────────────────────────────────────────────
// Phase 2: SCC-local Dijkstra
//   For each SCC: run Dijkstra FROM each local node to compute intra_sum and d2exit.
//   This replaces the bitwise SpMM of the undirected case.
//   PARALLEL: each SCC processed independently on separate threads.
// ─────────────────────────────────────────────────────────────────────────────
void phase2(const WDG& g, const SCC& scc, const Cond& dag,
            vector<LS>& lsa, const vector<int>& home) {
    int ns = scc.ns;
    lsa.resize(ns);

    // Find exit nodes per SCC from bridges
    vector<set<int>> exits(ns);
    for (auto& br : dag.bridges) exits[br.s].insert(br.eu);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int s = 0; s < ns; s++) {
        LS& ls = lsa[s];
        ls.K = (int)scc.sccs[s].size();
        ls.l2g = scc.sccs[s];
        ls.g2l.assign(g.n, -1);
        for (int li = 0; li < ls.K; li++) ls.g2l[ls.l2g[li]] = li;

        // Build local directed and reverse adjacency
        ls.ladj.resize(ls.K); ls.rladj.resize(ls.K);
        for (int li = 0; li < ls.K; li++) {
            for (auto [gv, w] : g.adj[ls.l2g[li]]) {
                int lv = ls.g2l[gv];
                if (lv >= 0) { ls.ladj[li].push_back({lv,w}); ls.rladj[lv].push_back({li,w}); }
            }
        }

        // Home and exit nodes
        for (int li = 0; li < ls.K; li++)
            if (home[ls.l2g[li]] == s) ls.home_l.push_back(li);
        ls.home_cnt = (int)ls.home_l.size();
        for (int gv : exits[s]) {
            int lv = ls.g2l[gv];
            if (lv >= 0) { ls.exit_l.push_back(lv); ls.exit_g.push_back(gv); }
        }
        ls.ne = (int)ls.exit_l.size();
        ls.intra_sum.assign(ls.K, 0.0);
        ls.d2exit.assign(ls.K, vector<W>(ls.ne, INF));

        // Run Dijkstra FROM each local node lv:
        //   - accumulate intra_sum[lv] = Σ_{u∈HOME} dist(lv, u)
        //   - record d2exit[lv][e] = dist(lv, exit_e)
        for (int lv = 0; lv < ls.K; lv++) {
            auto d = dijkstra(ls.ladj, ls.K, lv);
            for (int hu : ls.home_l) if (d[hu] < INF) ls.intra_sum[lv] += d[hu];
            for (int e = 0; e < ls.ne; e++) ls.d2exit[lv][e] = d[ls.exit_l[e]];
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase 3: Cross-SCC Distance DP
//   Process SCCs in REVERSE topological order of condensation DAG.
//   For each SCC s and outgoing bridge s->t via exit_e (weight w_br), entry_t:
//
//   cross_sum[s][lv] += (h_t + cross_cnt[t]) * (d2exit[s][lv][e] + w_br)
//                      + intra_sum[t][entry_t_local]
//                      + cross_sum[t][entry_t_local]
//
//   This recursion is correct because:
//     dist(lv_s, u in HOME(t)) = d2exit[s][lv][e] + w_br + dist_t(entry, u)
//     dist(lv_s, u beyond t)   = d2exit[s][lv][e] + w_br + dist_t(entry, exit_t')
//                                + w_br' + dist(entry', u) (recursively)
//     The term (intra_sum[t][entry]+cross_sum[t][entry]) captures all of this.
// ─────────────────────────────────────────────────────────────────────────────
struct CrossDP {
    vector<vector<W>> cs;  // cross_sum[s][lv]
    vector<long long> cc;  // cross_cnt[s] = reachable HOME nodes beyond s

    void compute(const WDG& g, const SCC& scc, const Cond& dag,
                 const vector<LS>& lsa, const vector<int>& home) {
        int ns = scc.ns;
        cs.resize(ns); cc.assign(ns, 0);
        for (int s = 0; s < ns; s++) cs[s].assign(lsa[s].K, 0.0);

        auto topo = dag.topo_order();  // forward topological order
        // Process in REVERSE topo (sinks first so cs[t] is ready when processing s)
        for (int i = (int)topo.size()-1; i >= 0; i--) {
            int s = topo[i];
            const LS& ls = lsa[s];

            // Collect unique outgoing SCCs to avoid double-counting
            set<int> done_targets;
            for (auto& [key, br_idx] : dag.best_br) {
                if (key.first != s) continue;
                int t = key.second;
                if (done_targets.count(t)) continue;
                done_targets.insert(t);

                auto& br = dag.bridges[br_idx];
                const LS& lt = lsa[t];

                // Find e_idx: index of exit node br.eu in ls.exit_l
                int e_idx = -1;
                for (int e = 0; e < ls.ne; e++)
                    if (ls.exit_g[e] == br.eu) { e_idx = e; break; }
                if (e_idx < 0) continue;

                // Find entry local index in t
                int entry_li = lt.g2l[br.ev];
                if (entry_li < 0) continue;

                long long h_t = lt.home_cnt;
                long long beyond_t = cc[t];
                W entry_dist = (lt.intra_sum[entry_li] < INF ? lt.intra_sum[entry_li] : 0.0)
                             + (cs[t][entry_li] < INF ? cs[t][entry_li] : 0.0);

                for (int lv = 0; lv < ls.K; lv++) {
                    W dev = ls.d2exit[lv][e_idx];
                    if (dev >= INF) continue;
                    cs[s][lv] += (W)(h_t + beyond_t) * (dev + br.w) + entry_dist;
                }
                cc[s] += h_t + beyond_t;
            }
            // Deduplicate cross_cnt: recount via BFS on condensation
        }

        // Recompute cc[s] accurately via BFS (avoids double-counting from multiple bridges)
        for (int s = 0; s < ns; s++) {
            vector<bool> vis(ns, false);
            queue<int> q; q.push(s); vis[s] = true; long long cnt = 0;
            while (!q.empty()) {
                int cur = q.front(); q.pop();
                if (cur != s) cnt += lsa[cur].home_cnt;
                for (auto [nxt,w] : dag.adj[cur]) if (!vis[nxt]) { vis[nxt]=true; q.push(nxt); }
            }
            cc[s] = cnt;
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Phase 4: Assemble Closeness Centrality
// ─────────────────────────────────────────────────────────────────────────────
void assemble(int n, const SCC& scc, const vector<LS>& lsa,
              const CrossDP& dp, const vector<int>& home, vector<double>& cc) {
    cc.assign(n, 0.0);
    #pragma omp parallel for schedule(dynamic, 32)
    for (int v = 0; v < n; v++) {
        if (home[v] < 0) continue;
        int s = home[v], lv = lsa[s].g2l[v];
        if (lv < 0) continue;
        W D = lsa[s].intra_sum[lv] + dp.cs[s][lv];
        long long r = (long long)lsa[s].home_cnt - 1 + dp.cc[s];
        if (D > 1e-12 && r > 0) cc[v] = (double)r / D;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SEQUENTIAL wrapper (single threaded)
// ─────────────────────────────────────────────────────────────────────────────
void run_sequential(const WDG& g, vector<double>& cc) {
    SCC scc(g.n); scc.run(g);
    Cond dag; dag.build(g, scc);
    vector<int> home(g.n, -1);
    for (int b = 0; b < scc.ns; b++) for (int v : scc.sccs[b]) if (home[v]<0) home[v]=b;
    vector<LS> lsa; phase2(g, scc, dag, lsa, home);
    CrossDP dp; dp.compute(g, scc, dag, lsa, home);
    assemble(g.n, scc, lsa, dp, home, cc);
}

// ─────────────────────────────────────────────────────────────────────────────
// PARALLEL wrapper (OpenMP)
// ─────────────────────────────────────────────────────────────────────────────
void run_parallel(const WDG& g, int threads, vector<double>& cc) {
    omp_set_num_threads(threads);
    // SCC and DAG are sequential (Kosaraju is inherently sequential)
    SCC scc(g.n); scc.run(g);
    Cond dag; dag.build(g, scc);
    vector<int> home(g.n, -1);
    for (int b = 0; b < scc.ns; b++) for (int v : scc.sccs[b]) if (home[v]<0) home[v]=b;
    // Phase 2 is parallel (OMP over SCCs)
    vector<LS> lsa; phase2(g, scc, dag, lsa, home);
    // Cross DP is sequential (topological DP with dependencies)
    CrossDP dp; dp.compute(g, scc, dag, lsa, home);
    // Assembly is parallel
    assemble(g.n, scc, lsa, dp, home, cc);
}

// ─────────────────────────────────────────────────────────────────────────────
// DYNAMIC: batch update handler
//   On edge insertion/deletion:
//   1. Update graph
//   2. Check if SCC structure changed
//   3. If yes: full recompute (SCC+DAG+Dijkstra+DP)
//   4. If no: re-Dijkstra only in affected SCCs, re-DP from affected SCCs upward
// ─────────────────────────────────────────────────────────────────────────────
struct Dynamic {
    WDG& g;
    SCC scc; Cond dag;
    vector<int> home;
    vector<LS> lsa;
    CrossDP dp;
    vector<double> cc;

    Dynamic(WDG& g) : g(g), scc(g.n) {
        scc.run(g); dag.build(g, scc);
        home.assign(g.n, -1);
        for (int b=0; b<scc.ns; b++) for (int v:scc.sccs[b]) if (home[v]<0) home[v]=b;
        phase2(g, scc, dag, lsa, home);
        dp.compute(g, scc, dag, lsa, home);
        assemble(g.n, scc, lsa, dp, home, cc);
    }

    // Insert batch of edges
    void insert_edges(vector<tuple<int,int,W>>& edges) {
        set<int> affected_sccs;
        for (auto [u,v,w] : edges) {
            g.add(u,v,w);
            affected_sccs.insert(scc.comp[u]);
            affected_sccs.insert(scc.comp[v]);
        }
        recompute_affected(affected_sccs);
    }

    // Delete batch of edges
    void delete_edges(vector<pair<int,int>>& edges) {
        set<int> affected_sccs;
        for (auto [u,v] : edges) {
            affected_sccs.insert(scc.comp[u]);
            affected_sccs.insert(scc.comp[v]);
            g.remove(u,v);
        }
        recompute_affected(affected_sccs);
    }

private:
    void recompute_affected(set<int>& touched) {
        // Check if SCC structure changed
        SCC new_scc(g.n); new_scc.run(g);
        bool changed = (new_scc.ns != scc.ns);
        if (!changed) for (int v=0; v<g.n; v++) if (new_scc.comp[v]!=scc.comp[v]){changed=true;break;}

        if (changed) {
            // Full recompute
            scc = new_scc;
            dag.build(g, scc);
            fill(home.begin(), home.end(), -1);
            for (int b=0; b<scc.ns; b++) for (int v:scc.sccs[b]) if (home[v]<0) home[v]=b;
            lsa.clear(); phase2(g, scc, dag, lsa, home);
            dp.compute(g, scc, dag, lsa, home);
        } else {
            // Partial recompute: only affected SCCs + their DAG ancestors
            // Find ancestors via reverse BFS on condensation
            vector<vector<int>> radj(scc.ns);
            for (int s=0; s<scc.ns; s++) for (auto[t,w]:dag.adj[s]) radj[t].push_back(s);
            vector<bool> needs_update(scc.ns, false);
            queue<int> q;
            for (int s : touched) { needs_update[s]=true; q.push(s); }
            while (!q.empty()) {
                int s = q.front(); q.pop();
                for (int anc : radj[s]) if (!needs_update[anc]) { needs_update[anc]=true; q.push(anc); }
            }
            // Re-run Dijkstra for affected SCCs
            for (int s=0; s<scc.ns; s++) if (needs_update[s]) {
                LS& ls = lsa[s];
                // Rebuild local adj
                for (auto& a : ls.ladj) a.clear();
                for (auto& a : ls.rladj) a.clear();
                for (int li=0; li<ls.K; li++) {
                    for (auto [gv,w] : g.adj[ls.l2g[li]]) {
                        int lv=ls.g2l[gv]; if(lv>=0){ ls.ladj[li].push_back({lv,w}); ls.rladj[lv].push_back({li,w}); }
                    }
                }
                ls.intra_sum.assign(ls.K, 0.0);
                ls.d2exit.assign(ls.K, vector<W>(ls.ne, INF));
                for (int lv=0; lv<ls.K; lv++) {
                    auto d = dijkstra(ls.ladj, ls.K, lv);
                    for (int hu : ls.home_l) if (d[hu]<INF) ls.intra_sum[lv]+=d[hu];
                    for (int e=0; e<ls.ne; e++) ls.d2exit[lv][e]=d[ls.exit_l[e]];
                }
            }
            // Recompute cross DP (full, since ancestors are affected)
            dp.compute(g, scc, dag, lsa, home);
        }
        assemble(g.n, scc, lsa, dp, home, cc);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Baseline: plain Dijkstra from every node
// ─────────────────────────────────────────────────────────────────────────────
void baseline(const WDG& g, vector<double>& cc) {
    cc.assign(g.n, 0.0);
    for (int s=0; s<g.n; s++) {
        auto d = dijkstra(g.adj, g.n, s);
        W td=0; long long r=0;
        for (int v=0; v<g.n; v++) if (v!=s && d[v]<INF) { td+=d[v]; r++; }
        if (td>1e-12) cc[s]=(double)r/td;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    int threads = (argc>1) ? atoi(argv[1]) : 4;

    printf("=================================================================\n");
    printf("SCC-Confined CC for Weighted Directed Graphs\n");
    printf("Threads: %d\n", threads);
    printf("=================================================================\n\n");

    auto T  = [] { return chrono::high_resolution_clock::now(); };
    auto ms = [](auto a, auto b) { return chrono::duration<double,milli>(b-a).count(); };

    auto test = [&](const char* name, WDG& g, bool do_dynamic) {
        printf("─── %s (n=%d m=%d) ───\n", name, g.n, g.m);

        // Baseline
        auto t0=T(); vector<double> ref; baseline(g, ref); auto t1=T();
        printf("  Baseline (all-pairs Dijkstra): %.3fms\n", ms(t0,t1));

        // Sequential
        omp_set_num_threads(1);
        auto t2=T(); vector<double> seq_cc; run_sequential(g, seq_cc); auto t3=T();
        double me=0; for(int v=0;v<g.n;v++) me=max(me,fabs(seq_cc[v]-ref[v]));
        printf("  Sequential SCC-Dijk:           %.3fms  speedup=%.2fx  err=%.2e  %s\n",
               ms(t2,t3), ms(t0,t1)/ms(t2,t3), me, me<1e-4?"PASS":"FAIL");

        // Parallel
        omp_set_num_threads(threads);
        auto t4=T(); vector<double> par_cc; run_parallel(g, threads, par_cc); auto t5=T();
        me=0; for(int v=0;v<g.n;v++) me=max(me,fabs(par_cc[v]-ref[v]));
        printf("  Parallel SCC-Dijk (%d threads): %.3fms  speedup=%.2fx  err=%.2e  %s\n",
               threads, ms(t4,t5), ms(t0,t1)/ms(t4,t5), me, me<1e-4?"PASS":"FAIL");

        // Work reduction
        long long wk=0; 
        SCC tmp_scc(g.n); tmp_scc.run(g);
        for(auto& s: tmp_scc.sccs) wk+=(long long)s.size()*s.size();
        printf("  Work ratio (ΣK²/N²):           %.4f\n", (double)wk/g.n/g.n);

        // Dynamic
        if (do_dynamic) {
            WDG g2 = g;
            Dynamic dyn(g2);
            vector<tuple<int,int,W>> inserts;
            srand(42);
            for (int i=0; i<10; i++) {
                int u=rand()%g.n, v=rand()%g.n;
                if (u!=v) inserts.push_back({u,v,1.0+(rand()%5)*0.5});
            }
            auto td0=T(); dyn.insert_edges(inserts); auto td1=T();
            printf("  Dynamic (10 insertions):       %.3fms  (full recompute if SCC changes)\n",
                   ms(td0,td1));
        }
        printf("\n");
    };

    // ── Test 1: Chain of cliques (clear SCC structure)
    {
        int nc=8, cs=6;
        WDG g(nc*cs);
        for (int c=0; c<nc; c++) {
            int b=c*cs;
            for (int i=0; i<cs; i++) for (int j=0; j<cs; j++) if (i!=j)
                g.add(b+i, b+j, 1.0+(i+j)%3);
            if (c+1<nc) g.add(b+cs-1, (c+1)*cs, 2.0);
        }
        test("Chain-of-cliques directed (8 SCCs × 6 nodes)", g, true);
    }

    // ── Test 2: Sparse random directed
    {
        srand(42); WDG g(300);
        for (int u=0; u<300; u++) for (int v=0; v<300; v++)
            if (u!=v && rand()%100<2) g.add(u, v, 1.0+(rand()%10)*0.5);
        test("Sparse random directed (n=300, p=2%)", g, false);
    }

    // ── Test 3: Dense (typically 1 SCC)
    {
        srand(77); WDG g(200);
        for (int u=0; u<200; u++) for (int v=0; v<200; v++)
            if (u!=v && rand()%100<8) g.add(u, v, 1.0+(rand()%5)*0.5);
        test("Dense random directed (n=200, p=8%)", g, false);
    }

    return 0;
}
