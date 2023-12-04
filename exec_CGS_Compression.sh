#!/bin/bash
rm a.out
g++ -O3 CGS_compression.cpp -w

# $1 represents dataset name
# $2 represents the CGS Algorithm i.e, {I,E,U} for {CCN-I, CCN-E, CCN-U} respectively.
# $3 represents the threshold value [0,1]

./a.out "./Testing/$1" "./Compression_Results/$1" $1 $2 $3 &

echo "Code running in the backgroud. Results will be saved in the 'Compression_Results' folder."