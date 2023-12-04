#include<bits/stdc++.h>
using namespace std;

#define vvi vector< vector<int> >
#define vi vector< int >
#define vf vector< float >
#define vii vector< pair<int, int> >
#define REP(A, B) for(int A = 0; A < B; A++)
#define REPN(A, B, C) for(int A = B; A < C; A++)

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

    vi NQ(int i);
    int SP(int i, int j); // Shortest path
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

vi Graph::NQ(int i) {
  if(!is_complement) return adj[i];
  
  vi result;
  vector< bool > v(n, 0);
  REP(k, adj[i].size()) v[adj[i][k]] = 1;
  REP(k, n) if(!v[k] && k != i) result.push_back(k);

  return result;
}

int Graph::SP(int i, int j) {
    vector<bool> visited(n, false);
    vi parent(n, -1);
    queue<int> q;
    int path = 0;

    visited[i] = true;
    q.push(i);

    while(!q.empty()) {
      int node = q.front();
      q.pop();
      if(node == j) break;

      vi N = NQ(node); 
      REP(k, N.size()) {
        int neigh = N[k];
        if(!visited[neigh]) {
            visited[neigh] = true;
            parent[neigh] = node;
            if(neigh == j) {
              q.empty();
              break;
            }
            q.push(neigh);
        }
      }
    }

    // If path does not exist, return -1
    if(!visited[j]) return -1;

    //  Calculate path length
    int node = j;
    while( parent[node] !=-1) {
      path += 1;
      node = parent[node];
    }
    return path; 
}

//-----------------------------------------------------------------------
//------------------------- Compressed Graph ----------------------------
//-----------------------------------------------------------------------

class CompressedGraph {
  private:
  public:
    int n_, m, n;
    vvi adj, org, par, cached_sN;
    bool is_complement;
    vi d;

    void input(char* EL_filename, char* SN_filename);
    void get_parent();
    int size();
    vi NQ(int i);
    vi SN(int i);
    vi sN_cached(int i);
    int SP(int i, int j); // Shortest path
};

void CompressedGraph::input(char* EL_filename, char* SN_filename) {

  // Storing the edge list of the compressed graph
  ifstream in_1(EL_filename);
  streambuf *cinbuf_1 = std::cin.rdbuf(); //save old buf
  cin.rdbuf(in_1.rdbuf()); //redirect std::cin to in.txt!

  cin >> n_ >> m >> n;
  adj.assign(n_, vi ());

  int a, b;
  REP(i, m) { 
      cin >> a >> b;
      adj[a].push_back(b); 
      adj[b].push_back(a);
  }

  d.assign(n_, 0);

  std::cin.rdbuf(cinbuf_1);   //reset to standard input again

  // Storing supernodes in org variable
  org.assign(n_, vi ());
  ifstream in_2(SN_filename);
  streambuf *cinbuf_2 = std::cin.rdbuf(); //save old buf
  cin.rdbuf(in_2.rdbuf()); //redirect std::cin to in.txt!

  is_complement = false;

  int x, y, z;
  REPN(i, n, n_) {
      cin >> x >> y >> z;
      org[x].push_back(y);
      org[x].push_back(z);
  }

  std::cin.rdbuf(cinbuf_2);   //reset to standard input again
}

vi CompressedGraph::NQ(int i) { // Get the neighborhood of node i
    vi temp;
    vi neighbourhood;    
    vi supernodes = SN(i); // Get the supernodes containing the node i, if nothing, it will give the node itself.

    REP(k, supernodes.size()) { // Add the neighborhood of the supernodes, to retreive the neighborhood of the node i
        int s = supernodes[k];
        REP(j, adj[s].size()) temp.push_back(adj[s][j]);
    }

    vector<bool> v(n, 0);

    REP(k, temp.size()) {
        int node = temp[k];
        // If the node is not a supernode, then add it to the neighborhood
        if(node < n && !v[node]) neighbourhood.push_back(node), v[node] = 1;
        else {
             // If the node is a supernode, then add its subnodes to the neighborhood
            vi subnodes = sN_cached(node);
            REP(j, subnodes.size()) {
                int s = subnodes[j];
                if(!v[s]) neighbourhood.push_back(s), v[s] = 1;
            }
        }
    }

    // If not complemented, return the neighborhood as it is.
    if(!is_complement) return neighbourhood; 

    // If complemented, return the complement of the neighborhood as it is.
    vi result;
    REP(k, n) if(!v[k] && k != i) result.push_back(k);
    return result;
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

void CompressedGraph::get_parent() {
  // Reverse the parent->child retaion stored in the variable "org" to obtain child->parent relation "par"
  par.assign(n_, vi ());
  cached_sN.assign(n_, vi ()); // to store subnodes
  REP(i, org.size()) {
    vi temp = org[i];
    REP(j, temp.size()) {
      int node = temp[j];
      par[node].push_back(i);
      if(node<n) cached_sN[i].push_back(node); 
      else {
        vi sub_node = cached_sN[node];
        REP(k, sub_node.size()) cached_sN[i].push_back(sub_node[k]);
      }
    }
  }
}

vi CompressedGraph::sN_cached(int i) {
  if(i < n) return vi(1, i); // If not a supernode, then return the node itself.
  else return cached_sN[i]; // If a supernode, then return its subnodes.
}

int CompressedGraph::size() {
    int size = 0;
  
    REP(i, n_) size += adj[i].size();                                       // 2*|E|
    REPN(i, n, n_) size += 2;                                               // Mapping Cost
    REP(i, n_) if(adj[i].size() != 0) size += 1;                                   // |V|

    return size; 
}

int CompressedGraph::SP(int i, int j) {
    vector<bool> visited(n_, false);
    vi parent(n_, -1);
    queue<int> q;
    int path = 0;
    bool flag = 0;

    visited[i] = true;
    q.push(i);

    while(!q.empty()) {
      int node = q.front();
      q.pop();
      
      if(node == j) flag = true;
      else {
        vi sub_nodes = sN_cached(node);
        REP(l, sub_nodes.size()) if(sub_nodes[l] == j) flag = true;
      }
      if(flag == true) break;
      
      vi N = NQ(node);
      REP(k, N.size()) {
        int neigh = N[k];
        if(!visited[neigh]) {
            visited[neigh] = true;
            parent[neigh] = node;
            //  Checking if the node == j or if the node contains j
            if(neigh == j) flag = true;
            else {
              vi sub_nodes = sN_cached(neigh);
              REP(l, sub_nodes.size()) {
                if(!visited[sub_nodes[l]]) {
                  visited[sub_nodes[l]] = true;
                  parent[sub_nodes[l]] = node;
                  if(sub_nodes[l] == j) flag = true;
                }
              }
            }

            if(flag == true) {
              q.empty();
              break;
            }
            q.push(neigh);
        }
      }
    }

    // If path does not exist, return 0
    if(!flag) return -1;

    //  Calculate path length
    int node = j;
    while( parent[node] !=-1) {
      path += 1;
      node = parent[node];
    }

    return path; 
}

//-----------------------------------------------------------------------
//------------------------- Output Controller ---------------------------
//-----------------------------------------------------------------------

class OutputController {
  private:
    int fp;
    char *dataset, *compressed_EL_filename, *compressed_SN_filename, *algo, *filename;
    
    Graph G;
    CompressedGraph G_;

    double RQ_time;
    double compressed_RQ_time;
    
    void reachability_query(int k);
    void RQ_times();
    void find_path();
    void DFS(int v, vi &visited, bool compressed);
    vii get_node_pairs(int n, int k);

  public:
    OutputController(char* dataset_, char* compressed_EL_filename_, char* compressed_SN_filename_, Graph original, CompressedGraph compressed, char* algo_, char* outfile);
    void print_analysis();
};  
OutputController::OutputController(char* dataset_, char* compressed_EL_filename_, char* compressed_SN_filename_, Graph original, CompressedGraph compressed, char* algo_, char* outfile) {
  G = original;
  G_ = compressed;
  dataset = dataset_;
  compressed_EL_filename = compressed_EL_filename_;
  compressed_SN_filename = compressed_SN_filename_;
  algo = algo_;

  filename = outfile;
  RQ_time = 0;
  compressed_RQ_time = 0;

  if(*algo == 'I')
    fp = -1;
  else if(*algo == 'E')
    fp = 0;
  else if(*algo == 'U')
    fp = 1;
  
}

void OutputController::DFS(int v, vi &visited, bool compressed) {
  if(!compressed) {

    visited[v] = 1;
    vi GN = G.NQ(v);
    REP(i, GN.size()) {
      int node = GN[i];
      if(!visited[node]) 
        DFS(node, visited, compressed);
    }
  }
  else {
    visited[v] = 1;
    vi G_N = G_.NQ(v);
    REP(i, G_N.size()) {
      int node = G_N[i];
      if(!visited[node]) 
        DFS(node, visited, compressed);
    }
  }
}

void OutputController::reachability_query(int k) {
  clock_t start, end;
  double avg = 0.0, maxx = 0.0, minn = 0.0, std = 0.0;
  int count =0, count_G = 0, count_G_ = 0;

  // --------- Decompress the graph to G_d ------------------
  
  // Initializing graph data structure
  Graph G_d;
  start = clock();
  G_d.n = G.n;
  G_d.m = 0;
  G_d.adj.assign(G_d.n, vi ());
  G_d.is_complement = G.is_complement;

  // ------------------------------------------------------------

  // Reachability querying using shortest path querying 

  vii node_pairs = get_node_pairs(G.n, k);

  REP(i, k){
    
    int node1 = node_pairs[i].first;
    int node2 = node_pairs[i].second;

    // On original graph
    start = clock();
    int path_G = G.SP(node1, node2);
    end = clock();
    RQ_time += double(end - start) / (double(CLOCKS_PER_SEC));

    // On compressed graph (by local decompression)
    start = clock();
    int path_G_ = G_.SP(node1, node2);
    end = clock();
    compressed_RQ_time += double(end - start) / (double(CLOCKS_PER_SEC));

    // Computing reachability query
    if(path_G == -1 && path_G_ == -1) count +=1;
    else {
      if(path_G == -1) count_G +=1;
      if(path_G_ == -1) count_G_ +=1;
    }
  }

  RQ_time = RQ_time /(1.0 * k);
  compressed_RQ_time = compressed_RQ_time /(1.0 * k);
  
  double reachability_accuracy = ((k - count) - count_G_)/(k - count);
  
  float cf = (G.size - G_.size()) / (1.0 * G.size);
  float cr = 1.0 - cf;
  
  cout<<"\n\nOriginal graph size (G): "<<G.size;
  cout<<"\nCompressed graph size (G_c): "<<G_.size();
  cout<<"\nCompression Ratio: "<<cr;

  cout<<"\n\nReachability query results on "<< k <<" random node pairs:";
  cout<<"\n----------------------------------------------------";
  
  cout<<"\nReachability Accuracy: "<< reachability_accuracy;
  cout<<"\nNo. of paths not existing in both G & G_c: "<<count;
  cout<<"\nNo. of paths not existing in G: "<<count_G;
  cout<<"\nNo. of paths not existing in G_c: "<<count_G_;
  
}

vii OutputController::get_node_pairs(int n, int k) {
  vii pairs;

  while(pairs.size() != k) {
    pair<int, int> node_pair;
    int i = rand() % n;
    int j = rand() % n;
    while (i == j) i = rand() % n;

    if(i<j) {
      node_pair.first = i;
      node_pair.second = j;
    }
    else {
      node_pair.first = j;
      node_pair.second = i;
    }

    bool exists = false;
    REP(m, pairs.size()) {
      if(pairs[m].first == node_pair.first)
        if(pairs[m].second == node_pair.second)
          exists = true;
    }
    if(!exists) pairs.push_back(node_pair);
  }
  return pairs;
}

void OutputController::RQ_times() {
  cout<< "\nReachability query time on G (in seconds): "<< RQ_time;
  cout<< "\nReachability query time on G_c (in seconds): "<< compressed_RQ_time;
}

void OutputController::print_analysis() {
  std::ofstream out(filename, ios::app);
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
  G_.get_parent();
  reachability_query(100000); // computes shortest path query for 100000 node pairs
  RQ_times();
  cout << "\n\n"; 
  std::cout.rdbuf(coutbuf); // reset to standard output again  

}

int main(int argc, char* argv[]) {

  char* filename = argv[1]; // input file path
  char* compressed_EL_filename = argv[2]; // compressed edge list file path
  char* compressed_SN_filename = argv[3]; // compressed supernode file path
  char* algo = argv[4]; // compression type
  char* outfile = argv[5]; // output file path

  Graph G; // creating graph data structure
  G.input(filename); // storing the input graph as adjacency list

  CompressedGraph G_; // creating compressed graph data structure
  G_.input(compressed_EL_filename, compressed_SN_filename); // storing the compressed graph

  OutputController OutClient(filename, compressed_EL_filename, compressed_SN_filename, G, G_, algo, outfile);
  OutClient.print_analysis();
  
  return 0;
}