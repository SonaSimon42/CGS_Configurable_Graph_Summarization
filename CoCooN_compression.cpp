#include<bits/stdc++.h>
using namespace std;

#define vvi vector< vector<int> >
#define vi vector< int >
#define vf vector< float >
#define vii vector< pair<int, int> >
#define REP(A, B) for(int A = 0; A < B; A++)
#define REPN(A, B, C) for(int A = B; A < C; A++)

vi remove(vi vec, int i, int j) {
    vi temp;
    REP(k, vec.size()) if(vec[k] != i && vec[k] != j) temp.push_back(vec[k]);
    vec.clear();
    return temp;
}
vvi intersection(int n, vi I, vi J) {
  vector< bool > v(n, 0);
  vi inter;
  vi ninter_i, ninter_j;

  REP(k, I.size()) v[I[k]] = 1;

  REP(k, J.size()) {
      int node = J[k];
      if(v[node]) inter.push_back(node);
      else ninter_j.push_back(node);
      v[node] = 0;
  }

  REP(k, I.size()) if(v[I[k]]) ninter_i.push_back(I[k]);

  vvi result;
  result.assign(3, vi ());
  result[0] = inter;
  result[1] = ninter_i;
  result[2] = ninter_j;

  return result;
}
vi Union(int n, vi I, vi J) {
  vector< bool > v(n, 0);
  vi result;

  REP(k, I.size()) {
    int node = I[k];
    if(node > n) cout << "Exception triggered at Union" << endl;
    v[node] = 1;
    result.push_back(node);
  }
  REP(k, J.size()) {
    int node = J[k];
    if(node > n) cout << "Exception triggered at Union" << endl;
    if(!v[node]) result.push_back(node);
  }

  return result;
}
class Heap {
  private:
  public:
    vector< pair<int, pair<int, int> > > tree;
    map< pair<int, int> , int > hash;
    vi index;

    Heap();
    void push(int value, pair<int, int> node);
    void insert(int value, pair<int, int> node);
    void make_heap();
    pair<int, pair<int, int> > top();
    void pop();
    void swap(int i, int j);
    void heapify_up(int i);
    void heapify_down(int i);
    bool empty();
    bool exists(pair<int, int> node);
    void decrease_key(pair<int, int> node, int value);
    void increase_key(pair<int, int> node, int value);

    void print(pair<int, int> node);
    void print_all();

};
Heap::Heap() {
    tree.assign(1, {0, {0, 0}});
    index.assign(1, -1);
}
bool Heap::empty() {
    if(tree.size() <= 1) return true;
    else return false;
}
void Heap::push(int value, pair<int, int> node) {

    if(hash[ node ] != 0) return;

    int key = index.size();
    hash[ node ] = key;

    tree.push_back({value, node});
    index.push_back(tree.size() - 1);

    heapify_up(tree.size() - 1);
}
void Heap::insert(int value, pair<int, int> node) {
    if(hash[ node ] != 0) return;

    int key = index.size();
    hash[ node ] = key;

    tree.push_back({value, node});
    index.push_back(tree.size() - 1);
}
void Heap::make_heap() {
    REPN(i, 1, tree.size()) heapify_down(i);
}
void Heap::heapify_up(int i) {
    if(i == 1) return;

    int p = i / 2;
    pair<int, pair<int, int> > cur = tree[i];
    pair<int, pair<int, int> > par = tree[p];

    if ( cur.first > par.first ) {
        swap(i, p);
        heapify_up(p);
    }
}
void Heap::heapify_down(int i) {
    if(i >= tree.size()) return;

    int cl = 2*i;
    int cr = 2*i + 1;

    int swap_with = -1;

    if(cl < tree.size() && cr < tree.size()) {
        if( tree[cl].first >= tree[cr].first && tree[cl].first > tree[i].first ) swap_with = cl;
        else if( tree[cr].first > tree[cl].first && tree[cr].first > tree[i].first ) swap_with = cr;
    }
    else if( cl < tree.size() && tree[cl].first > tree[i].first ) swap_with = cl;
    else if( cr < tree.size() && tree[cr].first > tree[i].first ) swap_with = cr;

    if(swap_with != -1) {
        swap(i, swap_with);
        heapify_down(swap_with);
    }
}
void Heap::swap(int i, int j) {
    pair<int, pair<int, int> > temp_1 = tree[i];
    tree[i] = tree[j];
    tree[j] = temp_1;

    int idx_1 = hash[ tree[i].second ];
    int idx_2 = hash[ tree[j].second ];

    int temp_2 = index[idx_1];
    index[idx_1] = index[idx_2];
    index[idx_2] = temp_2;
}
pair<int, pair<int, int> > Heap::top() {
    if(empty()) return {0, {0, 0}};
    return tree[1];
}
void Heap::pop() {
    if(tree.size() > 2) {
        swap(1, tree.size() - 1);
        hash.erase( hash.find( tree[tree.size() - 1].second ) );
        tree.pop_back();
        heapify_down(1);
    }

    else if(tree.size() == 2) {
        hash.erase( hash.find( tree[1].second ) );
        tree.pop_back();
    }
}
void Heap::decrease_key(pair<int, int> node, int value) {
    int key = hash[ node ];
    if(key == 0) return;

    int idx = index[key];
    pair<int, pair<int, int> > old = tree[idx];
    tree[idx] = {old.first - value, old.second};

    heapify_down(idx);
}
bool Heap::exists(pair<int, int> node) {
    if(hash[ node ] == 0) return false;
    return true;
}
void Heap::increase_key(pair<int, int> node, int value) {
    int key = hash[ node ];
    if(key == 0) return;

    int idx = index[key];
    pair<int, pair<int, int> > old = tree[idx];
    tree[idx] = {old.first + value, old.second};

    heapify_up(idx);
}
void Heap::print(pair<int, int> node) {
    int key = hash[ node ];
    if(key == 0) return;
    
    pair< int, pair<int,int> > item = tree[index[key]];
    cout << item.second.first << ", " << item.second.second << " : " << item.first << "\n";
}
void Heap::print_all() {

    REP(i, index.size()) {
    pair< int, pair<int,int> > item = tree[index[i]];
    }
}

//-----------------------------------------------------------------------
//-------------------------- Graph Structure ----------------------------
//-----------------------------------------------------------------------


class Graph {
  private:
  public:
    int n, m, size;
    vi d;
    vvi adj;
    bool is_complement;

    void input(char* filename);
    void input(Graph G);
    void complement(Graph G);
};

void Graph::input(char* filename) {
    ifstream in(filename);
    streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    is_complement = false;

    cin >> n >> m;
    adj.assign(n, vi ());
    d.assign(n, 0);

    int a, b;
    REP(i, m) {
        cin >> a >> b;
        adj[a].push_back(b); 
        adj[b].push_back(a); 
    }
    REP(i, n) d[i] = adj[i].size();

    size = n + 2*m;
    std::cin.rdbuf(cinbuf);   //reset to standard input again
}
void Graph::complement(Graph G) {
    is_complement = true;
    n = G.n;
    m = 0;
    size = G.size;

    adj.assign(n, vi ());
    d.assign(n, 0);

    vector< bool > v(n, 0);
    REP(i, n) {
        REP(k, G.adj[i].size()) v[G.adj[i][k]] = 1;
        REP(k, n) {
            if(!v[k] && k != i) adj[i].push_back(k);
            v[k] = 0;
        }
    }

    REP(i, n) d[i] = adj[i].size();
    REP(i, n) m += (int) adj[i].size();
    m /= 2;
}
void Graph::input(Graph G) {
  n = G.n;
  m = G.m;
  is_complement = G.is_complement;
  size = G.size;

  adj.assign(n, vi ());
  d.assign(n, 0);

  REP(i, n) REP(k, G.adj[i].size()) adj[i].push_back(G.adj[i][k]);
  REP(i, n) d[i] = adj[i].size();
}

//-----------------------------------------------------------------------
//---------------------------- Compression ------------------------------
//-----------------------------------------------------------------------

class CompressedGraph {
  private:
  public:
    int n_, m, n;
    vi d;
    vvi adj, org, par, cached_sN;
    bool is_complement;
    vector< bool > exists;
    vector< bool > removed;

    void input(Graph G);

    vi N2(int i);
    vi SN(int i);
    vi sN_cached(int i);
    int size();
};

void CompressedGraph::input(Graph G) {
  n_ = n = G.n;
  m = G.m;
  is_complement = G.is_complement;

  adj.assign(n, vi ());
  org.assign(n, vi ());
  cached_sN.assign(n, vi());
  d.assign(n, 0);
  exists.assign(n, true);
  removed.assign(n, false);

  REP(i, n) REP(k, G.adj[i].size()) adj[i].push_back(G.adj[i][k]);
  REP(i, n) d[i] = adj[i].size();
}
vi CompressedGraph::N2(int i) { // Get two hop neighborhood of node i
    vector<bool> v(n_, 0);
    vi two_hop;

    v[i] = 1;
    REP(k, adj[i].size()) {
        int neigh = adj[i][k];
        REP(l, adj[neigh].size()) {
            int node = adj[neigh][l];
            if(!v[node]) two_hop.push_back(node);
            v[node] = 1;
        }
    }

    return two_hop;
}

vi CompressedGraph::SN(int i) {
    // Get all the supernodes containing the node i from the "par" table
    // If no supernodes, then the node itself is returned.
    vi supernodes;
    vector<bool> v(n_, 0);
    stack<int> st;
    st.push(i);

    while(!st.empty()) {
        int top = st.top();
        st.pop();
        if(!v[top]) supernodes.push_back(top), v[top] = 1;
        REP(k, par[top].size()) st.push(par[top][k]);
    }

    return supernodes;
}

vi CompressedGraph::sN_cached(int i) {
  if(i < n) return vi(1, i); // If not a supernode, then return the node itself.
  else return cached_sN[i]; // If a supernode, then return its subnodes.
}
int CompressedGraph::size() {
    int size = 0;

    REP(i, n_) size += adj[i].size();                                       // 2*|E|
    REPN(i, n, n_) size += 2;                                               // Mapping Cost
    REP(i, n_) if(adj[i].size() != 0) size += 1;                            // |V|

    return size; 
}

//-----------------------------------------------------------------------
//---------------------------- CommonHood -------------------------------
//-----------------------------------------------------------------------

class CommonHood {
  private:
    Heap heap;
    Graph G;
    CompressedGraph G_;

    int fp;
    vf det;

    void greedy_cn();
    void build_heap();
    int cg(int i, int j);
    bool is_safe_merge(int i, int j);
    void merge(int i, int j);
    void update_heap(int i, int j, int idx);
    void update_degree(int i, int j);

  public:

    double compression_time;
    double heap_construction_time;
    double iteration_time;

    CommonHood(Graph target, int type, vf thres);
    CompressedGraph call();
};

CommonHood::CommonHood(Graph target, int type, vf thres) {

  // initializing the time variable to zero
  compression_time = 0;
  heap_construction_time = 0;
  iteration_time = 0;

  int m = target.m;
  int n = target.n;

  G.input(target), fp = type, det = thres;
  // checking the complement graph
  if( 1ll * m <= (1ll * n * (n-1)) / 4 ) G.input(target), fp = type, det = thres;
  else {
    G.complement(target), fp = -type;
    det.assign(n, 0);
    REP(i, n) det[i] = thres[i] * G.d[i] / (n - (1 - thres[i]) * G.d[i]);
  }
}
CompressedGraph CommonHood::call() {
  clock_t start, end;
  start = clock();
  greedy_cn();
  end = clock();
  compression_time = double(end - start) / double(CLOCKS_PER_SEC);
  return G_;
}
void CommonHood::greedy_cn() {
  G_.input(G);
  clock_t start, end;
  start = clock();
  build_heap();
  heap.print_all();
  end = clock();
  heap_construction_time = double(end - start) / double(CLOCKS_PER_SEC);

  start = clock();
  while(!heap.empty()) {
    pair< int, pair<int,int> > top = heap.top();
    heap.pop();
    int gain = top.first;
    int i = top.second.first;
    int j = top.second.second;
    if(gain <= 0) break;
    bool flag = G_.exists[i] && G_.exists[j] && is_safe_merge(i, j);
    if(flag) {
      if(fp == -1) update_degree(i, j);
      merge(i, j);   
    } 
  }
  end = clock();
  iteration_time = double(end - start) / double(CLOCKS_PER_SEC);
}
void CommonHood::build_heap() {
  REP(i, G_.n_) {
    vi hop2 = G_.N2(i);

    REP(k, hop2.size()) {
      if( heap.exists( {hop2[k], i} )) continue;
      int gain = cg(i, hop2[k]);
      if(gain > 0) heap.push(gain, {i, hop2[k]});
    }
  }
}
int CommonHood::cg(int i, int j) {
  vvi sets = intersection(G_.n_, G_.adj[i], G_.adj[j]);
  vi inter = sets[0];
  vi ninter_i = sets[1];
  vi ninter_j = sets[2];

  bool ij_edge = false;
  REP(k, ninter_i.size()) if(ninter_i[k] == j) ij_edge = true; 
  REP(k, ninter_j.size()) if(ninter_j[k] == i) ij_edge = true; 

  int gain = 0;

  if(fp == 0)
    gain = 2 * inter.size() - 3;

  else if(fp == -1) {
    gain = 2 * (inter.size() + ninter_i.size() + ninter_j.size()) - 3;
    if(ij_edge) gain -= 2;
  }
  
  else if(fp == 1) 
    gain = 2 * inter.size() - 3;

  if(i < G_.n) {
    if(fp == -1) gain++;
    else if(fp == 1 && !ij_edge) gain++;
    else if(fp == 0 && ninter_i.size() == 0) gain++;
  }

  if(j < G_.n) {
    if(fp == -1) gain++;
    else if(fp == 1 && !ij_edge) gain++;
    else if(fp == 0 && ninter_j.size() == 0) gain++;
  }

  return gain;
}
bool CommonHood::is_safe_merge(int i, int j) {
  vvi results = intersection(G_.n_, G_.adj[i], G_.adj[j]);
  vi inter = results[0];
  vi ninter_i = results[1];
  vi ninter_j = results[2];

  vi nI, nJ;

  REP(k, ninter_i.size()) {
    int node = ninter_i[k];
    if(node < G_.n) nI.push_back(node);
    else {
      vi subnodes = G_.sN_cached(node);
      REP(l, subnodes.size()) nI.push_back(subnodes[l]);
    }
  }

  REP(k, ninter_j.size()) {
    int node = ninter_j[k];
    if(fp == -1 && node == i)  continue;
    if(node < G_.n) nJ.push_back(node);
    else {
      vi subnodes = G_.sN_cached(node);
      REP(l, subnodes.size()) nJ.push_back(subnodes[l]);
    }
  }

  vi I = G_.sN_cached(i);
  vi J = G_.sN_cached(j);

  if(fp == -1) {

    vi loss(G.n, 0);
    REP(k, I.size()) loss[I[k]] += nI.size();
    REP(k, nI.size()) loss[nI[k]] += I.size();
    REP(k, J.size()) loss[J[k]] += nJ.size();
    REP(k, nJ.size()) loss[nJ[k]] += J.size(); 

    REP(k, I.size()) {
      int node = I[k];
      if( (G_.d[node] - loss[node]) / (1.0 * G.d[node]) < 1 - det[node] ) return false;
    }
    REP(k, nI.size()) {
      int node = nI[k];
      if( (G_.d[node] - loss[node]) / (1.0 * G.d[node]) < 1 - det[node] ) return false;
    }

    REP(k, J.size()) {
      int node = J[k];
      if( (G_.d[node] - loss[node]) / (1.0 * G.d[node]) < 1 - det[node] ) return false; 
    }
    REP(k, nJ.size()) {
      int node = nJ[k];
      if( (G_.d[node] - loss[node]) / (1.0 * G.d[node]) < 1 - det[node] ) return false;
    }
  }

  else if(fp == 1) {

    nI = remove(nI, j, -1);
    nJ = remove(nJ, i, -1);

    REP(k, I.size()) {
      int node = I[k];
      if( G.d[node] / (1.0 * G_.d[node] + nJ.size()) < 1 - det[node] ) return false;
    }
    REP(k, nI.size()) {
      int node = nI[k];
      if( G.d[node] / (1.0 * G_.d[node] + J.size()) < 1 - det[node] ) return false;
    }
    
    REP(k, J.size()) {
      int node = J[k];
      if( G.d[node] / (1.0 * G_.d[node] + nI.size()) < 1 - det[node] ) return false;
    }
    REP(k, nJ.size()) {
      int node = nJ[k];
      if( G.d[node] / (1.0 * G_.d[node] + I.size()) < 1 - det[node] ) return false;
    }

  }
  return true;
}
void CommonHood::merge(int i, int j) {
  vvi results = intersection(G_.n_, G_.adj[i], G_.adj[j]);
  vi inter = results[0];
  vi ninter_i = results[1];
  vi ninter_j = results[2];

  if(inter.size() == 0) return;

  bool ij_edge = false;
  REP(k, ninter_i.size()) if(ninter_i[k] == j) ij_edge = true; 
  REP(k, ninter_j.size()) if(ninter_j[k] == i) ij_edge = true; 

  int idx = G_.adj.size();
  G_.adj.push_back(vi ());
  G_.org.push_back(vi ());
  G_.cached_sN.push_back(vi ());
  G_.exists.push_back(true);
  G_.n_++;

  // Adding neighbourhood of s.
  if(fp == 0 || fp == -1) REP(k, inter.size()) G_.adj[idx].push_back(inter[k]);
  else if(fp == 1) {
    REP(k, inter.size()) G_.adj[idx].push_back(inter[k]);
    REP(k, ninter_i.size()) G_.adj[idx].push_back(ninter_i[k]);
    REP(k, ninter_j.size()) G_.adj[idx].push_back(ninter_j[k]);
  }

  // Updating neighbourhoods of nehoughbours.
  if(fp == 0) {
    REP(k, inter.size()) {
      int node = inter[k];
      G_.adj[node] = remove(G_.adj[node], i, j);
      G_.adj[node].push_back(idx);
    }
  } else {
    REP(k, inter.size()) {
      int node = inter[k];
      G_.adj[node] = remove(G_.adj[node], i, j);
      G_.adj[node].push_back(idx);
    }

    REP(k, ninter_i.size()) {
      int node = ninter_i[k];
      G_.adj[node] = remove(G_.adj[node], i, j);
      if(fp == 1) G_.adj[node].push_back(idx);
    }

    REP(k, ninter_j.size()) {
      int node = ninter_j[k];
      G_.adj[node] = remove(G_.adj[node], i, j);
      if(fp == 1) G_.adj[node].push_back(idx);
    }

    G_.adj[idx] = remove(G_.adj[idx], i, j);
  }
  
  // Updating: neighbourhoods of i and j; heap; degrees
  if(fp == -1 || fp == 1) {
    update_heap(i, j, idx);
    if(fp == 1) update_degree(i, j);
    G_.adj[i].clear(), G_.exists[i] = false;
    G_.adj[j].clear(), G_.exists[j] = false;
    if(fp == 1 && ij_edge) G_.adj[i].push_back(j), G_.adj[j].push_back(i);
  }
  else if(fp == 0) {
    update_heap(i, j, idx);
    G_.adj[i].clear(), G_.adj[j].clear();
    G_.adj[i] = ninter_i;
    G_.adj[j] = ninter_j;
  }

  // Making the parent and child mappings
  G_.org[idx].push_back(i);
  G_.org[idx].push_back(j);

  // Combine the subnodes
  // Will combine the subodes present in node i and j
  // The vector at any index of cached_sN will not contain any supernodes
  G_.cached_sN[idx] = Union(G_.n, G_.cached_sN[i], G_.cached_sN[j]);
  if(i < G_.n) G_.cached_sN[idx].push_back(i);
  if(j < G_.n) G_.cached_sN[idx].push_back(j);

  // Removing isolated nodes.
  if(G_.adj[i].size() == 0 && i < G_.n) G_.removed[i] = true;
  if(G_.adj[j].size() == 0 && j < G_.n) G_.removed[j] = true;

  // Updating the number of edges
  if(fp == 0 || fp == 1) G_.m -= inter.size();
  else if(fp == -1) {
    G_.m -= (inter.size() + ninter_i.size() + ninter_j.size());
    if(ij_edge) G_.m +=1;
  }
}
void CommonHood::update_heap(int a, int b, int idx) {
  vi hop2 = G_.N2(idx);
  REP(k, hop2.size()) {
    int gain = cg(idx, hop2[k]);
    if(gain > 0) heap.push(gain, {idx, hop2[k]});
  }

  vvi results = intersection(G_.n_, G_.adj[a], G_.adj[b]);
  vi inter = results[0];
  vi ninter_i = results[1];
  vi ninter_j = results[2];

  if(fp == 0) {

    // for each pair (𝑢, 𝑣) such that either 𝑢 or 𝑣 or both lie in 𝑁𝐺𝑐 (𝑢∗ ) ∩ 𝑁𝐺𝑐 (𝑣∗ ) do
      //Decrease-Key (𝐻,𝑢,𝑣,2)

    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, inter.size()) {
        int node_j = inter[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
      }
    }
    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, ninter_i.size()) {
        int node_j = ninter_i[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
        if(heap.exists({node_j, node_i})) heap.decrease_key({node_j, node_i}, 2); 
      }
    }
    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, ninter_j.size()) {
        int node_j = ninter_j[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2);
        if(heap.exists({node_j, node_i})) heap.decrease_key({node_j, node_i}, 2);
      }
    }

    // for each pair (𝑢, 𝑣∗ ) ∈ 𝐻 do
      //𝑣𝑎𝑙𝑢𝑒 ← 2|𝑁𝐺𝑐 (𝑢∗) ∩ 𝑁𝐺𝑐 (𝑣∗) ∩ 𝑁𝐺𝑐 (𝑢)| Decrease-Key (𝐻,𝑢, 𝑣∗, 𝑣𝑎𝑙𝑢𝑒)
    // for each pair (𝑢∗, 𝑣) ∈ 𝐻 do
      // 𝑣𝑎𝑙𝑢𝑒 ← 2|𝑁𝐺𝑐 (𝑢∗) ∩ 𝑁𝐺𝑐 (𝑣∗) ∩ 𝑁𝐺𝑐 (𝑣)| Decrease-Key (𝐻,𝑢∗, 𝑣, 𝑣𝑎𝑙𝑢𝑒)


    vi neighbourhood;
    vector<bool> v(G_.n_, 0);
    REP(i, inter.size()) {
      int node = inter[i];
      REP(j, G_.adj[node].size()) {
        int neigh = G_.adj[node][j];
        if(v[neigh]) continue;
        v[neigh] = 1;
        neighbourhood.push_back(neigh);
      }
    }
    v.clear();
    v.assign(G_.n_, 0);
    REP(i, inter.size()) v[inter[i]] = 1;
    REP(i, neighbourhood.size()) {
      int cnt = 0;
      int node = neighbourhood[i];
      REP(j, G_.adj[node].size()) if(v[ G_.adj[node][j] ]) cnt++;
      if( heap.exists({a, node}) ) heap.decrease_key({a, node}, 2*cnt);
      if( heap.exists({b, node}) ) heap.decrease_key({b, node}, 2*cnt);
      if( heap.exists({node, a}) ) heap.decrease_key({node, a}, 2*cnt);
      if( heap.exists({node, b}) ) heap.decrease_key({node, b}, 2*cnt);
    }
  }

  else if(fp == -1) {
    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, inter.size()) {
        int node_j = inter[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
      }
    }
    REP(i, ninter_i.size()) {
      int node_i = ninter_i[i];
      REP(j, ninter_i.size()) {
        int node_j = ninter_i[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
      }
    }
    REP(i, ninter_j.size()) {
      int node_i = ninter_j[i];
      REP(j, ninter_j.size()) {
        int node_j = ninter_j[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
      }
    }
    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, ninter_i.size()) {
        int node_j = ninter_i[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
        if(heap.exists({node_j, node_i})) heap.decrease_key({node_j, node_i}, 2); 
      }
    }
    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, ninter_j.size()) {
        int node_j = ninter_j[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2);
        if(heap.exists({node_j, node_i})) heap.decrease_key({node_j, node_i}, 2);
      }
    }
  }

  else if(fp == 1) {

    REP(i, ninter_i.size()) {
      int node_i = ninter_i[i];
      REP(j, ninter_j.size()) {
        int node_j = ninter_j[j];
        if(heap.exists({node_i, node_j})) heap.increase_key({node_i, node_j}, 2); 
      }
    }

    REP(i, inter.size()) {
      int node_i = inter[i];
      REP(j, inter.size()) {
        int node_j = inter[j];
        if(heap.exists({node_i, node_j})) heap.decrease_key({node_i, node_j}, 2); 
      }
    }
  }
}
void CommonHood::update_degree(int i, int j) {
  vvi results = intersection(G_.n_, G_.adj[i], G_.adj[j]);
  vi inter = results[0];
  vi ninter_i = results[1];
  vi ninter_j = results[2];

  vi nI, nJ;
  REP(k, ninter_i.size()) {
    int node = ninter_i[k];
    if(node < G_.n) nI.push_back(node);
    else {
      vi subnodes = G_.sN_cached(node);
      REP(l, subnodes.size()) nI.push_back(subnodes[l]);
    }
  }
  REP(k, ninter_j.size()) {
    int node = ninter_j[k];
    if(fp == -1 && node == i) continue;
    if(node < G_.n) nJ.push_back(node);
    else {
      vi subnodes = G_.sN_cached(node);
      REP(l, subnodes.size()) nJ.push_back(subnodes[l]);
    }
  }
  vi I = G_.sN_cached(i);
  vi J = G_.sN_cached(j);

  if (fp == -1) {
    REP(k, I.size()) G_.d[I[k]] -= nI.size();
    REP(k, nI.size()) G_.d[nI[k]] -= I.size();
    REP(k, J.size()) G_.d[J[k]] -= nJ.size();
    REP(k, nJ.size()) G_.d[nJ[k]] -= J.size();
  }

  else if(fp == 1) {
    nI = remove(nI, j, -1);
    nJ = remove(nJ, i, -1);
  
    REP(k, I.size()) G_.d[I[k]] += nJ.size();
    REP(k, nI.size()) G_.d[nI[k]] += J.size();
    REP(k, J.size()) G_.d[J[k]] += nI.size();
    REP(k, nJ.size()) G_.d[nJ[k]] += I.size();
  }

}

//-----------------------------------------------------------------------
//------------------------- Output Controller ---------------------------
//-----------------------------------------------------------------------

class OutputController {
  private:
    int fp;
    float tol;
    char *algo, *dataset, *filename, *output_graph_EL, *output_graph_SN;
    vf error;
    
    Graph G;
    CompressedGraph G_;
    CommonHood *Compressor;

    void compression();
    void compression_times();
    void save_compressed_graph();
    void save_supernode_graph();

  public:
    OutputController(Graph original, CompressedGraph compressed, CommonHood *compressionClient, char* outfile, char* dataset, char* _algo, int _fp, float _tol, char* _output_graph_EL, char* _output_graph_SN);
    void print_analysis();
};  
OutputController::OutputController(Graph original, CompressedGraph compressed, CommonHood *compressionClient, char* outfile, char* _dataset, char* _algo, int _fp, float _tol, char* _output_graph_EL, char* _output_graph_SN){
  G = original;
  G_ = compressed;
  Compressor = compressionClient;
  filename = outfile;
  dataset = _dataset;
  algo = _algo;
  fp = _fp;
  tol = _tol;
  output_graph_EL = _output_graph_EL;
  output_graph_SN = _output_graph_SN;
}
void OutputController::compression() {

  float cf = (G.size - G_.size()) / (1.0 * G.size);
  float cr = 1.0 - cf;
  float degree = G.m / (1.0 * G.n);
  cout<<"\nNodes in input graph (G): "<<G.n;
  cout<<"\nEdges in input graph (G): "<<G.m;
  cout<<"\nNodes in compressed graph (G_c): "<<G_.n_;
  cout<<"\nEdges in compressed graph (G_c): "<<G_.m;
  cout<<"\nSize of original graph (G): "<<G.size;
  cout<<"\nSize of compressed graph (G_c): "<<G_.size();
  cout<<"\nCompression Ratio: "<<cr;
}

void OutputController::compression_times() {
  cout<< "\nCompression time (in seconds): "<<Compressor->compression_time;
  cout<< "\nHeap construction time (in seconds): "<< Compressor->heap_construction_time;
  cout<< "\nIteration time (in seconds): "<< Compressor->iteration_time;

}

void OutputController::save_compressed_graph() {
    cout<< G_.n_ << "\t" << G_.m<< "\t" << G_.n <<"\n";
    REP(i, G_.n_) {
      REP(j, G_.adj[i].size()) {
        int node_a = i;
        int node_b = G_.adj[i][j];
        if(node_a > node_b)
          cout<< node_a << "\t" << node_b<<"\n";
      }
    }                                    
}

void OutputController::save_supernode_graph() {
  REPN(i, G_.n, G_.org.size()) {
    vi temp = G_.org[i];
    cout<< i;
    REP(j, temp.size()) {
      int node = temp[j];
      cout<<"\t" << node;
    }
    cout<<"\n";
  }     
}

void OutputController::print_analysis() {

  std::ofstream out(filename, ios::app);
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
  cout <<"Dataset: " << dataset << "\nAlgorithm: " << algo << "\nThreshold: " << tol;
  compression();
  compression_times();
  cout << "\n\n";
  std::cout.rdbuf(coutbuf); // reset to standard output again  
  
  std::ofstream out_file1(output_graph_EL, ios::trunc);
  std::streambuf *coutbuf1 = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out_file1.rdbuf()); //redirect std::cout to out.txt!
  save_compressed_graph();
  std::cout.rdbuf(coutbuf1); //reset to standard output again 

  std::ofstream out_file2(output_graph_SN, ios::trunc);
  std::streambuf *coutbuf2 = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out_file2.rdbuf()); //redirect std::cout to out.txt!
  save_supernode_graph();
  std::cout.rdbuf(coutbuf2); //reset to standard output again 
}

int main(int argc, char* argv[]) {

  char* filename = argv[1]; // input file path
  char* outfile = argv[2]; // output file path
  char* dataset = argv[3]; // input file name 
  char* algo = argv[4]; // compression type
  int fp;
  if(*algo == 'I')
    fp = -1;
  else if(*algo == 'E')
    fp = 0;
  else if(*algo == 'U')
    fp = 1;
  else
    cout<<"Invalid CoCooN Algorithm";
  float tol = atof(argv[5]); // tolerance

  string EL = "Compressed/EL_"+string(argv[4])+"_"+string(argv[5])+"_"+string(dataset); // compressed graph edge list
  char* output_graph_EL = &EL[0];
  string SN = "Compressed/SN_"+string(argv[4])+"_"+string(argv[5])+"_"+string(dataset); // compressed graph supernode list
  char* output_graph_SN = &SN[0];
  
  Graph G; // creating graph data structure
  G.input(filename); // storing the input graph as adjacency list
  vf det(G.n, tol); // setting the tolerance for each node (same tolerance)
  
  CommonHood CompressionClient(G, fp, det); // initializing the time variables and checking complement graph
  CompressedGraph G_ = CompressionClient.call(); // calling greedy_cn algorithm
  
  OutputController OutClient(G, G_, &CompressionClient, outfile, dataset, algo, fp, tol, output_graph_EL, output_graph_SN);
  OutClient.print_analysis();
  
  return 0;
}
