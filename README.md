# CGS - Configurable Graph Summarization with Bounded Neighborhood Loss and Query Support
## Running CGS
<ol>
<li> <strong>Input file format</strong>:
 The input file is a text file. The first line contains the number of nodes(V) and the number of edges(E) separated by space. Then, the next E lines represent the undirected edges, denoted by two integers separated by space. The nodes must be named from 0 to V-1. Save the graphs to be tested as text files in the <em>Dataset</em> folder.
</li>

<li><strong>Output file format</strong>:</li>
The edge list of the compressed graph is saved to the text file named 'EL_var_δ_dataset_name.txt' where the first line contains the number of nodes(V) and the number of edges(E) separated by a tab, and the next E lines represent the undirected edges, which is denoted by two integers separated by a tab.<br>
The supernode mapping of the compressed graph is saved to the text file named 'SL_var_δ_dataset_name.txt' where each line contains the supernode and its two subnodes separated by a tab.<br>The files are saved at <em>Compressed</em> folder.

<li><strong>Running the algorithm:</strong>

`./exec_CGS_Compression.sh input_file var δ `

* input_file: the input graph file in .txt format.
* var: I for CGS-I, E for CGS-E, and U for CGS-U.
* δ: threshold parameter for controlling the loss, a value between 0 and 1.
* Output files will be stored in <em>Compressed</em> folder.
* Output results will be stored in <em>Compressed_Results</em> folder.

Example:

To run the algorithm:

`./exec_CGS_Compression.sh AP.txt I 0.5`

Output saved at Compression_Results/AP.txt:

 <em>Dataset: AP.txt<br>
CGS Variant: I<br>
Threshold: 0.5<br>
Nodes in input graph (G): 18772<br>
Edges in input graph (G): 198050<br>
Nodes in compressed graph (G_c): 24351<br>
Edges in compressed graph (G_c): 86225<br>
Size of original graph (G): 414872<br>
Size of compressed graph (G_c): 196800<br>
Compression Ratio: 0.474363<br>
Compression time (in seconds): 131.31<br>
Heap construction time (in seconds): 24.56<br>
Iteration time (in seconds): 106.74</em><br>
</li>
</ol>

## Querying on Compressed Graph
<ol>

<li><strong>Running Neighbourhood Query:</strong>

`./exec_CGS_NQ.sh dataset_name edge_list.txt supernode_list.txt var`

* input_file: the input graph file in .txt format.
* var: I for CGS-I, E for CGS-E, and U for CGS-U.
* Results will be will be stored in <em>NQ_Results</em> folder as text file.

Example:

To run the algorithm:

`./exec_CGS_NQ.sh AP EL_I_0.5_AP.txt SN_I_0.5_AP.txt I`

Output saved at NQ_Results/AP.txt:

<em>
CGS Variant: I<br>
Original graph size: 414872<br>
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
Neighbourhood query time in G (in seconds): 0<br>
Neighbourhood query time in G_c (in seconds): 4.26167e-06
</em><br>
</li>

<li><strong>Running Shortest Path Query:</strong>

`./exec_CGS_SP.sh dataset_name edge_list.txt supernode_list.txt var`

* input_file: the input graph file in .txt format.
* var: I for CGS-I, E for CGS-E, and U for CGS-U.
* Results will be will be stored in <em>SP_Results</em> folder as text file.

Example:

To run the algorithm:

`./exec_CGS_SP.sh AP EL_I_0.5_AP.txt SN_I_0.5_AP.txt I`

Output saved at SP_Results/AP.txt:

<em>
CGS Variant: I<br>
Original graph size (G): 414872<br>
Compressed graph size (G_c): 196800<br>
Compression Ratio: 0.474363<br>

Shortest path query results on 100000 random node pairs:<br>
----------------------------------------------------<br>
Average error in shortest path: 0.147192<br>
Minimum error in shortest path: 0<br>
Maximum error in shortest path: 6<br>
Standard deviation of error in shortest path: 0.301521<br>
Shortest path query time on G (in seconds): 0.0020145<br>
Shortest path query time on G_c (in seconds): 0.0096335<br>
Decompression time of G_c (in seconds): 0.04<br>
Shortest path query time on globally decompressed G_c (in seconds): 0.0016715
</em><br>
</li>

<li><strong>Running Reachability Query:</strong>

`./exec_CGS_RQ.sh dataset_name edge_list.txt supernode_list.txt var`

* input_file: the input graph file in .txt format.
* var: I for CGS-I, E for CGS-E, and U for CGS-U.
* Results will be will be stored in <em>RQ_Results</em> folder as text file.

Example:

To run the algorithm:

`./exec_CGS_RQ.sh AP EL_I_0.5_AP.txt SN_I_0.5_AP.txt I`

Output saved at RQ_Results/AP.txt:

<em>CGS Variant: I<br>
Original graph size (G): 414872<br>
Compressed graph size (G_c): 196800<br>
Compression Ratio: 0.474363<br>

Reachability query results on 100000 random node pairs:<br>
----------------------------------------------------<br>
Reachability Accuracy: 0.942396<br>
No. of paths not existing in both G & G_c: 9086<br>
No. of paths not existing in G: 0<br>
No. of paths not existing in G_c: 5237<br>
Reachability query time on G (in seconds): 0.0020053<br>
Reachability query time on G_c (in seconds): 0.0098913
</em><br>
</li>

</ol>
