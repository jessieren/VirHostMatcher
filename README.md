Sequence comparison oligonucleotide frequency (ONF)
===========

Basic tools for computing various oligonucleotide frequency (ONF) based distance/dissimilarity measures


Requirements
---------------

The source code is written by C++. Thus it requires a C++ compiler. It has been tested on Mac and Linux, but it should work under Windows, Linux or Mac environment. 


Usage
---------------

This program is used to compute various oligonucleotide frequency (ONF) based distance/dissimilarity measures between a pair of DNA sequences. These measures include Euclidean distance (Eu), Manhattan distance (Ma), Chebyshev distance (Ch), Jensen-Shannon divergence (JS), d2, d2*, d2S, Hao, Teeling, EuF and Willner. See paper "Alignment-free d2* oligonucleotide frequency dissimilarity measure improves accuracy of predicting virus-host interactions" for the definitions. 

To use the tool, please simply follow the steps and copy and paste the following commands to the terminal command line. You can find an folder named "test" containing 2 phage sequences and 3 host sequences in fasta format. Here we use this test data to show how to use the tool.

* Step 0: parameter settings

We set the kmer length k=6, and the Markov chain order=2.

k=6
order=2


* Step 1: path settings

Let codeDIR be the directory where the two c++ scripts locate.

codeDIR=/Users/jessie/Desktop/alignment-free/script/

Users need to put the fasta sequences under a directory. The phage sequences and the host sequences can be in different folders. Let phageFaDIR be the path to the phage fasta files and hostFaDIR be the path to the host fasta files respectively.

phageFaDIR=/Users/jessie/Desktop/alignment-free/script/test/phage
hostFaDIR=/Users/jessie/Desktop/alignment-free/script/test/host

An output directory need to be created and set, 

outDIR=/Users/jessie/Desktop/alignment-free/script/test/output
mkdir $outDIR

* Step 2: compile the two c++ scripts

Use a C++ compiler to compile the two c++ script countKmer.cpp and computeMeasure.cpp.

g++ $codeDIR/countKmer.cpp -o $codeDIR/countKmer.out
g++ $codeDIR/computeMeasure.cpp -o $codeDIR/computeMeasure.out

* Step 3: count kmer frequency and compute the various measures

Two files containing the list of phages and the list of hosts are going to be generated. These files will be used in the final step: computing the various measures. 

> $outDIR/hostList
> $outDIR/phageList

Now let count the 1-6-mers for each of the phage and host fasta sequences.

for fa in $phageFaDIR/*
do
name=`basename $fa`
kmerDIR=$outDIR/kmerCount/$name
for klength in $(seq 1 $k)
do
$codeDIR/countKmer.out -l -k $klength -i $fa -o $kmerDIR
done
echo $name $kmerDIR $order >> $outDIR/phageList
done

for fa in $hostFaDIR/*
do
name=`basename $fa`
kmerDIR=$outDIR/kmerCount/$name
for klength in $(seq 1 $k)
do
$codeDIR/countKmer.out -l -k $klength -i $fa -o $kmerDIR
done
echo $name $kmerDIR $order >> $outDIR/hostList
done

Now we can finally compute the various distance/dissimialrity measures

$codeDIR/computeMeasure.out -k $k -i $outDIR/phageList -j $outDIR/hostList > $outDIR/results.csv


* Congratulations! The results can be find in $outDIR/results.csv. 



Contacts and bug reports
------------------------
Jie Ren
renj@usc.edu

Fengzhu Sun
fsun@usc.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
fixed.
2. Check that your input is in the correct format and you have selected the
correct options.
3. Please reduce your input to the smallest possible size that still produces
the bug; we will need your input data to reproduce the problem, and the
smaller you can make it, the easier it will be.


Copyright and License Information
---------------------------------
Copyright (C) 2016 University of Southern California, Jie Ren

Authors: Jie Ren

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

