#!/bin/bash
rm a.out
g++ -O3 CoCooN_RQ.cpp -w

# $1 represents dataset name
# $2 represents filename of edge list of compressed graph
# $3 represents filename of supernode list of compressed graph
# $4 represents the type of compression done

./a.out "./Testing/$1" "./Compressed/$2" "./Compressed/$3" "$4" "./RQ_Results/$1"&

echo "Code running in the backgroud. Results will be saved in the 'RQ_Results' folder."