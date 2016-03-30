Sequence comparison oligonucleotide frequency (ONF)
===========

Basic tools for computing various oligonucleotide frequency (ONF) based distance/dissimialrity measures


Requirements
---------------

The source code is written by C++. Thus it requires a C++ compiler. It has been tested on Mac and Linux, but it should work under Windows, Linux or Mac environment. 


Usage
---------------

This program is used to compute various oligonucleotide frequency (ONF) based distance/dissimialrity measures between a pair of DNA sequences. These measures include Euclidian distance (Eu), Manhattan distance (Ma), Chebyshev distance (Ch), Jensen-Shannon divergence (JS), d2, d2*, d2S, Hao, Teeling, EuF and Willner. See paper "Alignment-free d2* oligonucleotide frequency dissimilarity measure improves accuracy of predicting virus-host interactions" for the definitions. 


1. parameter setting
> k=6
> order=2

2. path setting
2.1 directory to the c++ code
> codeDIR=/Users/jessie/Desktop/alignment-free/script/
2.2 directories to the phage and host fasta files ##
> phageFaDIR=/Users/jessie/Desktop/alignment-free/script/test/phage
> hostFaDIR=/Users/jessie/Desktop/alignment-free/script/test/host
## directory to the output ##
> outDIR=/Users/jessie/Desktop/alignment-free/script/test/output
> mkdir $outDIR

#### compile the code ####
> g++ $codeDIR/countKmer.cpp -o $codeDIR/countKmer.out
> g++ $codeDIR/computeMeasure.cpp -o $codeDIR/computeMeasure.out

#### copy the following to the command line ####
## generate the file containing the phage and host file paths ##
> > $outDIR/hostList
> > $outDIR/phageList

## count the 1-6-mers for phage and host fasta sequences ##
> for fa in $phageFaDIR/*
do
name=`basename $fa`
kmerDIR=$outDIR/kmerCount/$name
for klength in $(seq 1 $k)
do
$codeDIR/countKmer.out -l -k $klength -i $fa -o $kmerDIR
done
echo $name $kmerDIR $order >> $outDIR/phageList
done

>for fa in $hostFaDIR/*
do
name=`basename $fa`
kmerDIR=$outDIR/kmerCount/$name
for klength in $(seq 1 $k)
do
$codeDIR/countKmer.out -l -k $klength -i $fa -o $kmerDIR
done
echo $name $kmerDIR $order >> $outDIR/hostList
done

## compute distance/dissimialrity using various measures ##
> $codeDIR/computeMeasure.out -k $k -i $outDIR/phageList -j $outDIR/hostList > $outDIR/results.csv

#### the results can be find in $outDIR/results.csv done! ####



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