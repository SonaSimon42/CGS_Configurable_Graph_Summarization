# CoCooN: A Graph Compression Framework with Bounded Neighborhood Loss

## Running CoCooN
<ol>
<li> <strong>Input file format</strong>:
 The input file is a text file. The first line contains the number of nodes(V) and the number of edges(E) separated by space. Then, the next E lines represent the undirected edges, which is denoted by two integers separated by space. The nodes must be named from 0 to V-1. Save the graphs to be tested as text files in the <em>Data</em> folder.
</li>

<li><strong>Output file format</strong>:
 The output of the compression algorithm is written to a text file in the following order. <br>
 <em>Dataset name, var, Threshold,	No. of nodes G,	No. of edges G,	no_of_nodes_G',	No. of edges G', G size, G' size,	Compression ratio, False positive,	false_negative, avg_loss_per_node,	min_loss_per_node,	max_loss_per_node,	Compression time,	Heap construction time,	Iteration time,	Neighbourhood query time,	Compressed neighborhood query time, No. of connected components in G,	No. of connected components in G_, Avg error in shortest path,	Max error in shortest path,	Min error in shortest path,	Standard deviation of error in shortest path,	No. of paths not existing in G & G_,	No. of paths not existing in G,	No. of paths not existing in G_,	Reachability accuracy,	No. of pairs with zero error,	Shortest path query time,	Compressed shortest path query time,	Decompression time,	Decompressed shortest path query time</em></li>

<li><strong>Running the algorithm:</strong>

`./exec_CoCooN.sh input_file var δ `

* input_file: the input graph file in .txt format.
* var: ls will be -1 for CCN-I, 0 for CCN-E, and 1 for CCN-U.
* δ: threshold parameter for controlling the loss.
* Output file will be stored in <em>Results</em> folder.

</li>
</ol>
