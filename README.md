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

## Querying on Compressed Graph
<ol>

<li><strong>Running Neighbourhood Query:</strong>

`./exec_CoCooN_NQ.sh dataset_name edge_list.txt supernode_list.txt var`

* input_file: the input graph file in .txt format.
* var: I for CCN-I, E for CCN-E, and U for CCN-U.
* Results will be will be stored in <em>NQ_Results</em> folder as text file.

Example:

To run the algorithm:

`./exec_CoCooN_NQ.sh AP EL_I_0.5_AP.txt SN_I_0.5_AP.txt I`

Output saved at NQ_Results/AP.txt:

<em>Original graph size: 414872<br>
Compressed graph size: 196800<br>
Compression Ratio: 0.474363<br>

Neighbourhood Query:<br>
---------------------<br>
No. of extra edges: 0<br>
No. of missing edges: 99262<br>
Average loss per node: 0.251521<br>
Minimum loss per node: 0<br>
Maximum loss per node: 0.5<br>
Standard deviation of loss per node: 0.176506<br>
Neighbourhood query time in G (in seconds): 1.06542e-06<br>
Neighbourhood query time in G_c (in seconds): 3.19625e-06</em><br>
</li>

<li><strong>Running Shortest Path Query:</strong>

`./exec_CoCooN_SP.sh dataset_name edge_list.txt supernode_list.txt var`

* input_file: the input graph file in .txt format.
* var: I for CCN-I, E for CCN-E, and U for CCN-U.
* Results will be will be stored in <em>SP_Results</em> folder as text file.

Example:

To run the algorithm:

`./exec_CoCooN_SP.sh AP EL_I_0.5_AP.txt SN_I_0.5_AP.txt I`

Output saved at SP_Results/AP.txt:

<em>Original graph size (G): 414872<br>
Compressed graph size (G_c): 196800<br>
Compression Ratio: 0.474363<br>

Shortest path query results on 100000 random node pairs:<br>
----------------------------------------------------<br>
Average error in shortest path: 0.147192<br>
Minimum error in shortest path: 0<br>
Maximum error in shortest path: 6<br>
Standard deviation of error in shortest path: 0.301521<br>
Shortest path query time on G (in seconds): 0.0021431<br>
Shortest path query time on G_c (in seconds): 0.0103938<br>
Decompression time of G_c (in seconds): 0.04<br>
Shortest path query time on globally decompressed G_c (in seconds): 0.0018237</em><br>
</li>

<li><strong>Running Reachability Query:</strong>

`./exec_CoCooN_RQ.sh dataset_name edge_list.txt supernode_list.txt var`

* input_file: the input graph file in .txt format.
* var: I for CCN-I, E for CCN-E, and U for CCN-U.
* Results will be will be stored in <em>RQ_Results</em> folder as text file.

Example:

To run the algorithm:

`./exec_CoCooN_RQ.sh AP EL_I_0.5_AP.txt SN_I_0.5_AP.txt I`

Output saved at RQ_Results/AP.txt:

<em>Original graph size (G): 414872<br>
Compressed graph size (G_c): 196800<br>
Compression Ratio: 0.474363<br>

Reachability query results on 100000 random node pairs:<br>
----------------------------------------------------<br>
Reachability Accuracy: 0<br>
No. of paths not existing in both G & G_c: 9086<br>
No. of paths not existing in G: 0<br>
No. of paths not existing in G_c: 5237<br>
Reachability query time on G (in seconds): 0.0021387<br>
Reachability query time on G_c (in seconds): 0.0106238<br>
</em><br>
</li>

</ol>
