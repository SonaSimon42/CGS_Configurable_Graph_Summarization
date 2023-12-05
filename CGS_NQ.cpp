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
    REP(i, n_) if(adj[i].size() != 0) size += 1;                            // |V|

    return size; 
}


//-----------------------------------------------------------------------
//------------------------- Output Controller ---------------------------
//-----------------------------------------------------------------------

class OutputController {
  private:
    int fp;
    char *dataset, *compressed_EL_filename, *compressed_SN_filename, *algo, *filename;
    vf error;
    
    Graph G;
    CompressedGraph G_;

    double neighbourhood_query_time;
    double compressed_neighbourhood_query_time;

    void neighbourhood_loss();
    void per_node_loss();
    void neighbourhood_times();

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
  neighbourhood_query_time = 0;
  compressed_neighbourhood_query_time = 0;

  if(*algo == 'I')
    fp = -1;
  else if(*algo == 'E')
    fp = 0;
  else if(*algo == 'U')
    fp = 1;
  
}
void OutputController::neighbourhood_loss() {
  clock_t start, end;
  int fp = 0, fn = 0;

  REP(i, G.n) {
    start = clock();
    vi GN = G.NQ(i);
    end = clock();
    neighbourhood_query_time += double(end - start) / (double(CLOCKS_PER_SEC) * G.n);

    start = clock();
    vi G_N = G_.NQ(i);
    end = clock();
    compressed_neighbourhood_query_time += double(end - start) / (double(CLOCKS_PER_SEC) * G.n);

    vvi results = intersection(G.n, GN, G_N);
    fp += results[2].size();
    fn += results[1].size();
  }
  float cf = (G.size - G_.size()) / (1.0 * G.size);
  float cr = 1.0 - cf;
  
  cout<<"\n\nOriginal graph size: "<<G.size;
  cout<<"\nCompressed graph size: "<<G_.size();
  cout<<"\nCompression Ratio: "<<cr;
  cout<<"\n\nNeighbourhood Query:";
  cout<<"\n---------------------";
  cout<<"\nNo. of extra edges: "<<fp/2;
  cout<<"\nNo. of missing edges: "<<fn/2;
}

void OutputController::per_node_loss() {
  double minn = 10000, maxx = 0;
  double avg = 0, loss, std =0;
  int sum;

  error.empty();
  REP(i, G.n) {
    G_.d[i] = G_.NQ(i).size();
    if(fp == -1) loss = ( (G.d[i] == 0) ? 0 : 1 - G_.d[i] / (1.0 * G.d[i]) );
    else         loss = ( (G_.d[i] == 0) ? 0 : 1 - G.d[i] / (1.0 * G_.d[i]) );
    avg += loss;
    error.push_back(loss);
    minn = min(minn, loss); 
    maxx = max(maxx, loss);
  }
  avg = avg / (1.0 * G.n);

  // Computing standard deviation
  REP(i, G.n) { std += pow((error[i] - avg),2);}
  double var = std/(1.0 * G.n);
  std = sqrt(std/(1.0 * G.n));

  cout<<"\nAverage loss per node: "<< avg;
  cout<<"\nMinimum loss per node: "<< minn;
  cout<<"\nMaximum loss per node: "<< maxx;
  cout<<"\nStandard deviation of loss per node: "<< std;
}

void OutputController::neighbourhood_times() {
  cout<< "\nNeighbourhood query time in G (in seconds): "<<neighbourhood_query_time;
  cout<< "\nNeighbourhood query time in G_c (in seconds): "<< compressed_neighbourhood_query_time;
}

void OutputController::print_analysis() {
  std::ofstream out(filename, ios::app);
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
  cout <<"Dataset: " << dataset << "\nCGS Variant: " << algo;
  G_.get_parent();
  neighbourhood_loss();
  per_node_loss();
  neighbourhood_times();
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