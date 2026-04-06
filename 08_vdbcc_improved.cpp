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
#include <map>
#include <cstdint>
#include <cmath>
#include <omp.h>

using namespace std;

struct Graph {
    int n;
    vector<vector<int>> adj;
    Graph() : n(0) {}
    Graph(int n) : n(n), adj(n) {}
    void addEdge(int u, int v) { adj[u].push_back(v); adj[v].push_back(u); }
    void removeEdge(int u, int v) {
        adj[u].erase(find(adj[u].begin(), adj[u].end(), v));
        adj[v].erase(find(adj[v].begin(), adj[v].end(), u));
    }
    bool hasEdge(int u, int v) const {
        for (int w : adj[u]) if (w == v) return true;
        return false;
    }
    static Graph readFromFile(const string& f) {
        ifstream fin(f); int n, m; fin >> n >> m;
        Graph g(n);
        for (int i = 0; i < m; i++) { int u, v; fin >> u >> v; g.addEdge(u, v); }
        return g;
    }
    static Graph generateRandom(int n, int m) {
        Graph g(n); srand(42);
        for (int i = 1; i < n; i++) { int p = rand() % i; g.addEdge(i, p); }
        set<pair<int,int>> ex;
        for (int u = 0; u < n; u++) for (int v : g.adj[u]) ex.insert({min(u,v), max(u,v)});
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

struct BCCDecomp {
    int n, timer_;
    vector<int> disc, low, par;
    vector<bool> isArt;
    vector<vector<int>> bccs, nodeBCCs;
    stack<pair<int,int>> st;

    BCCDecomp(const Graph& g) : n(g.n), timer_(0), disc(n,-1),
        low(n,0), par(n,-1), isArt(n,false), nodeBCCs(n) {
        for (int i = 0; i < n; i++)
            if (disc[i] == -1) { dfs(g, i); flush(); }
        for (int i = 0; i < (int)bccs.size(); i++)
            for (int v : bccs[i]) nodeBCCs[v].push_back(i);
    }
    void flush() {
        if (st.empty()) return;
        vector<int> b;
        while (!st.empty()) { auto [u,v] = st.top(); st.pop(); b.push_back(u); b.push_back(v); }
        sort(b.begin(), b.end()); b.erase(unique(b.begin(), b.end()), b.end());
        bccs.push_back(b);
    }
    void extract(int u, int v) {
        vector<int> b;
        while (!st.empty()) {
            auto [eu,ev] = st.top(); st.pop(); b.push_back(eu); b.push_back(ev);
            if ((eu==u && ev==v)||(eu==v && ev==u)) break;
        }
        sort(b.begin(), b.end()); b.erase(unique(b.begin(), b.end()), b.end());
        if (!b.empty()) bccs.push_back(b);
    }
    void dfs(const Graph& g, int u) {
        disc[u] = low[u] = timer_++; int ch = 0;
        for (int v : g.adj[u]) {
            if (disc[v]==-1) {
                ch++; par[v]=u; st.push({u,v}); dfs(g,v);
                low[u] = min(low[u], low[v]);
                if ((par[u]==-1 && ch>1)||(par[u]!=-1 && low[v]>=disc[u])) {
                    isArt[u]=true; extract(u,v);
                }
            } else if (v!=par[u] && disc[v]<disc[u]) {
                st.push({u,v}); low[u]=min(low[u],disc[v]);
            }
        }
    }
};

struct LocalBCC {
    int id, localN;
    vector<int> nodes;
    map<int,int> g2l;
    vector<vector<int>> localAdj;
    vector<vector<int>> dist;
    vector<int> artPointsLocal;
    vector<int> artPointsGlobal;

    void build(const Graph& g, const vector<int>& nodeList, int bccId, const vector<bool>& isArt) {
        id = bccId;
        nodes = nodeList;
        localN = nodes.size();
        for (int i = 0; i < localN; i++) g2l[nodes[i]] = i;
        localAdj.assign(localN, {});
        for (int i = 0; i < localN; i++) {
            int u = nodes[i];
            for (int v : g.adj[u])
                if (g2l.count(v)) localAdj[i].push_back(g2l[v]);
        }
        for (int i = 0; i < localN; i++)
            if (isArt[nodes[i]]) { artPointsLocal.push_back(i); artPointsGlobal.push_back(nodes[i]); }
    }

    void computeAllPairsDist() {
        dist.assign(localN, vector<int>(localN, -1));
        for (int s = 0; s < localN; s++) {
            dist[s][s] = 0;
            queue<int> q; q.push(s);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int v : localAdj[u])
                    if (dist[s][v]==-1) { dist[s][v]=dist[s][u]+1; q.push(v); }
            }
        }
    }

    int globalDist(int gU, int gV) const {
        auto itU = g2l.find(gU), itV = g2l.find(gV);
        if (itU == g2l.end() || itV == g2l.end()) return -1;
        return dist[itU->second][itV->second];
    }
};

struct VDBCC {
    Graph g;
    int n, nThreads;
    vector<LocalBCC> localBCCs;
    BCCDecomp* bcc;
    vector<double> sumDist, cc;
    vector<int> nodeHomeBCC;

    VDBCC(const Graph& g_, int nThreads_) : g(g_), n(g_.n), nThreads(nThreads_),
        sumDist(n, 0.0), cc(n, 0.0), nodeHomeBCC(n, -1) {
        buildStructures();
        computeCC();
    }

    ~VDBCC() { delete bcc; }

    void buildStructures() {
        bcc = new BCCDecomp(g);
        localBCCs.resize(bcc->bccs.size());

        #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
        for (int i = 0; i < (int)bcc->bccs.size(); i++) {
            localBCCs[i].build(g, bcc->bccs[i], i, bcc->isArt);
            localBCCs[i].computeAllPairsDist();
        }

        for (int v = 0; v < n; v++) {
            if (!bcc->isArt[v] && !bcc->nodeBCCs[v].empty())
                nodeHomeBCC[v] = bcc->nodeBCCs[v][0];
            else if (bcc->isArt[v])
                nodeHomeBCC[v] = bcc->nodeBCCs[v][0];
        }
    }

    void computeCC() {
        fill(sumDist.begin(), sumDist.end(), 0.0);
        fill(cc.begin(), cc.end(), 0.0);

        vector<vector<int>> artDist(n);
        vector<int> artList;
        for (int v = 0; v < n; v++)
            if (bcc->isArt[v]) artList.push_back(v);

        for (int v : artList) artDist[v].resize(n, -1);

        #pragma omp parallel num_threads(nThreads)
        {
            vector<int> d(n);
            queue<int> q;
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < (int)artList.size(); i++) {
                int a = artList[i];
                fill(d.begin(), d.end(), -1);
                d[a] = 0; q.push(a);
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    for (int v : g.adj[u])
                        if (d[v]==-1) { d[v]=d[u]+1; q.push(v); }
                }
                for (int v = 0; v < n; v++) artDist[a][v] = d[v];
            }
        }

        #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
        for (int v = 0; v < n; v++) {
            double sd = 0;
            int homeBCC = nodeHomeBCC[v];
            if (homeBCC < 0) continue;
            auto& lb = localBCCs[homeBCC];
            int lv = lb.g2l.at(v);

            for (int i = 0; i < lb.localN; i++) {
                if (i != lv && lb.dist[lv][i] > 0)
                    sd += lb.dist[lv][i];
            }

            if (lb.artPointsGlobal.empty()) {
                sumDist[v] = sd;
                if (sd > 0) cc[v] = (double)(n-1) / sd;
                continue;
            }

            for (int u = 0; u < n; u++) {
                if (u == v) continue;
                if (lb.g2l.count(u)) continue;

                int best = n * 10;
                for (int ai = 0; ai < (int)lb.artPointsGlobal.size(); ai++) {
                    int aGlobal = lb.artPointsGlobal[ai];
                    int aLocal = lb.artPointsLocal[ai];
                    int dVA = lb.dist[lv][aLocal];
                    int dAU = artDist[aGlobal][u];
                    if (dVA >= 0 && dAU >= 0)
                        best = min(best, dVA + dAU);
                }
                if (best < n * 10) sd += best;
            }

            sumDist[v] = sd;
            if (sd > 0) cc[v] = (double)(n-1) / sd;
        }
    }

    struct EdgeUpdate { int u, v; bool isInsertion; };

    struct UpdateResult { int affectedBCCs; double timeMs; };

    UpdateResult dynamicUpdate(const vector<EdgeUpdate>& batch) {
        set<int> affectedBCCids;

        for (auto& e : batch) {
            if (e.isInsertion) g.addEdge(e.u, e.v);
            else g.removeEdge(e.u, e.v);
            for (int bid : bcc->nodeBCCs[e.u]) affectedBCCids.insert(bid);
            for (int bid : bcc->nodeBCCs[e.v]) affectedBCCids.insert(bid);
        }

        auto t0 = chrono::high_resolution_clock::now();

        delete bcc;
        bcc = new BCCDecomp(g);
        localBCCs.resize(bcc->bccs.size());

        for (int i = 0; i < (int)bcc->bccs.size(); i++) {
            localBCCs[i] = LocalBCC();
            localBCCs[i].build(g, bcc->bccs[i], i, bcc->isArt);
            localBCCs[i].computeAllPairsDist();
        }

        for (int v = 0; v < n; v++) {
            if (!bcc->nodeBCCs[v].empty())
                nodeHomeBCC[v] = bcc->nodeBCCs[v][0];
        }

        computeCC();

        auto t1 = chrono::high_resolution_clock::now();
        return {(int)affectedBCCids.size(),
                chrono::duration<double, milli>(t1 - t0).count()};
    }
};

void parallelMSBFS(const Graph& g, int nThreads, vector<double>& sumDist) {
    int n = g.n;
    const int BATCH = 64;
    int numBatches = (n + BATCH - 1) / BATCH;
    fill(sumDist.begin(), sumDist.end(), 0.0);
    vector<vector<double>> tSD(nThreads, vector<double>(n, 0.0));

    #pragma omp parallel num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
        vector<uint64_t> fr(n), vis(n), nf(n);
        #pragma omp for schedule(dynamic, 1)
        for (int batch = 0; batch < numBatches; batch++) {
            int bs = batch*BATCH, be = min(bs+BATCH, n), sz = be-bs;
            fill(fr.begin(), fr.end(), 0ULL);
            fill(vis.begin(), vis.end(), 0ULL);
            for (int i = 0; i < sz; i++) { fr[bs+i]=(1ULL<<i); vis[bs+i]=(1ULL<<i); }
            int lev = 0; bool act = true;
            while (act) {
                lev++; fill(nf.begin(), nf.end(), 0ULL); act = false;
                for (int u = 0; u < n; u++) {
                    if (!fr[u]) continue;
                    for (int v : g.adj[u]) nf[v] |= fr[u];
                }
                for (int v = 0; v < n; v++) {
                    nf[v] &= ~vis[v];
                    if (nf[v]) { act=true; vis[v]|=nf[v]; tSD[tid][v]+=(double)__builtin_popcountll(nf[v])*lev; }
                }
                swap(fr, nf);
            }
        }
    }
    #pragma omp parallel for num_threads(nThreads)
    for (int v = 0; v < n; v++) { double t=0; for (int i=0;i<nThreads;i++) t+=tSD[i][v]; sumDist[v]=t; }
}

vector<double> naive_cc(const Graph& g) {
    int n = g.n; vector<double> cc(n, 0.0);
    for (int s = 0; s < n; s++) {
        vector<int> d(n,-1); queue<int> q; d[s]=0; q.push(s);
        long long td=0; int r=0;
        while (!q.empty()) { int u=q.front(); q.pop(); for (int v:g.adj[u]) if(d[v]==-1){d[v]=d[u]+1;td+=d[v];r++;q.push(v);} }
        if (td>0) cc[s]=(double)r/td;
    }
    return cc;
}

int main(int argc, char* argv[]) {
    Graph g;
    int nThreads = omp_get_max_threads();
    if (argc > 1) g = Graph::readFromFile(argv[1]);
    else g = Graph::generateRandom(2000, 8000);
    if (argc > 2) nThreads = atoi(argv[2]);

    cout << "Method: VDBCC with BCC-local distances + art-point BFS" << endl;
    cout << "Threads: " << nThreads << " | Graph: " << g.n << " nodes" << endl;

    auto t0 = chrono::high_resolution_clock::now();
    auto cc_naive = naive_cc(g);
    auto t1 = chrono::high_resolution_clock::now();
    double naive_ms = chrono::duration<double, milli>(t1-t0).count();

    t0 = chrono::high_resolution_clock::now();
    VDBCC engine(g, nThreads);
    t1 = chrono::high_resolution_clock::now();
    double vdbcc_ms = chrono::duration<double, milli>(t1-t0).count();

    vector<double> msbfs_sd(g.n);
    t0 = chrono::high_resolution_clock::now();
    parallelMSBFS(g, nThreads, msbfs_sd);
    t1 = chrono::high_resolution_clock::now();
    double msbfs_ms = chrono::duration<double, milli>(t1-t0).count();

    double maxErr = 0;
    for (int i = 0; i < g.n; i++) maxErr = max(maxErr, abs(cc_naive[i] - engine.cc[i]));

    cout << "Naive sequential:  " << fixed << setprecision(1) << naive_ms << " ms" << endl;
    cout << "Parallel MS-BFS:   " << msbfs_ms << " ms" << endl;
    cout << "VDBCC (BCC-local):  " << vdbcc_ms << " ms" << endl;
    cout << "Correctness vs naive: " << scientific << setprecision(1) << maxErr
         << (maxErr < 1e-9 ? " MATCH" : " MISMATCH") << endl;

    BCCDecomp bcc(g);
    int maxBCC = 0;
    for (auto& b : bcc.bccs) maxBCC = max(maxBCC, (int)b.size());
    int numArt = count(bcc.isArt.begin(), bcc.isArt.end(), true);
    cout << "BCCs: " << bcc.bccs.size() << " | Arts: " << numArt
         << " | Largest BCC: " << maxBCC << endl;

    srand(123);
    cout << "\n=== Dynamic Updates ===" << endl;

    for (int b = 0; b < 5; b++) {
        vector<VDBCC::EdgeUpdate> batch;
        for (int i = 0; i < 10; i++) {
            int u = rand() % g.n, v = rand() % g.n;
            while (v == u || engine.g.hasEdge(u, v)) v = rand() % g.n;
            batch.push_back({u, v, true});
        }

        auto res = engine.dynamicUpdate(batch);

        auto cc_check = naive_cc(engine.g);
        double err = 0;
        for (int i = 0; i < g.n; i++) err = max(err, abs(engine.cc[i] - cc_check[i]));

        cout << "Batch " << b+1 << ": " << res.affectedBCCs << " BCCs affected | "
             << fixed << setprecision(1) << res.timeMs << " ms | "
             << "err=" << scientific << setprecision(1) << err
             << (err < 1e-9 ? " OK" : " ERR") << endl;
    }

    return 0;
}
