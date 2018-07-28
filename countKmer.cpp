#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include <iostream>
#include <fstream>
#include <cstring>
////#include<cmath>
#include <vector>
//
#include <algorithm>
#include <unordered_map>
//using namespace std::tr1;
//#include<ext/hash_map>
//using __gnu_cxx::hash_map;


#include <libgen.h>
using namespace std;



//unsigned long power;//parameter, power = 4^(k-1)

//receptacle of kmer counts
unordered_map<unsigned long,unsigned long> HashTable;
//
int seqlength=999999999;
// seq: save read/genome(A C G T) in each line from file
char seq[999999999];
// inverse of seq
char seq_inverse[999999999];




/* Define scientific numbers and operations */


// In order to give more precise calculation and avoid zero in dominator,
// we introduce the SCIENTIFIC_NUMBER calculation.
// We will re-define addtion, multiplication and power in the following sub rountine functions.
// It is contributed by Prof. Minghua Deng, Peking University.

struct SCIENTIFIC_NUMBER
{
	int factor;
	double value;
};



// Scientific Number calculation funtion:

// TransToReal : trans a SCIENTIFIC NUMBER to a read number
double TransToReal(SCIENTIFIC_NUMBER dSci)
{
	double dReal=0;
	dReal = dSci.value * pow(10.0,dSci.factor);
	return dReal;
	
}


// TransToScientific : trans a read number to a SCIENTIFIC NUMBER
SCIENTIFIC_NUMBER TransToScientific(double dReal)
{
	SCIENTIFIC_NUMBER sciTemp;
	int count;
	
	if( dReal==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
	}
	else if(dReal>10.0 || dReal<-10.0)
	{
		count=0;
		while(dReal>10.0 || dReal<-10.0)
		{
			dReal /=10.0;
			count++;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else if( dReal<1.0 && dReal>-1.0)
	{
		count=0;
		while( dReal<1.0 && dReal>-1.0 )
		{
			dReal *=10.0;
			count--;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else
	{
		sciTemp.value=dReal;
		sciTemp.factor=0;
	}
	
	return sciTemp;
}


// SciMultiple : Multiplication of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
	//    cout << "SciMultiple " << endl;
	
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int count;
	
	if( left.value==0.0 || right.value==0.0 )
	{
		//        cout << "Both 0 " << endl;
		
		sciTemp.value=0.0;
		sciTemp.factor=0;
		
		return sciTemp;
	}
	
	// now both left and right element are nonzero
	dTemp=left.value * right.value;
	
	//    cout << "left.value " << left.value << endl;
	//    cout << "right.value " << right.value << endl;
	//    cout << "dTemp " << dTemp << endl;
	
	
	if( dTemp>10.0 || dTemp<-10.0 )
	{
		
		//        cout << "10 < dTemp or dTemp < -10 " << endl;
		
		count=0;
		while(dTemp>10.0 || dTemp<-10.0 )
		{
			dTemp /=10.0;
			count++;
		}
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
	}
	else if( dTemp<1.0 && dTemp>-1.0)
	{
		//        cout << "dTemp < 1 or dTemp > -1 " << dTemp << endl;
		
		count=0;
		while( dTemp<1.0 && dTemp>-1.0 )
		{
			dTemp *=10.0;
			count--;
		}
		
		//        cout << "count " << count << endl;
		
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
		
		//        cout << "sciTemp " << sciTemp.value << " " << sciTemp.factor << endl;
	}
	else
	{
		
		//        cout << "dTemp normal " << dTemp << endl;
		
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor;
	}
	
	return sciTemp;
}


// SciMultiple : Multiplication between a SCIENTIFIC NUMBERS and a read number
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
	
	sciTemp=TransToScientific(right);
	sciTemp=SciMultiple(left,sciTemp);
	
	return sciTemp;
}





// SciAddition : addition of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int i,count;
	
	if( left.value==0.0 || right.value==0.0 )
	{
		if( left.value==0.0 )
			return right;
		else
			return left;
	}
	
	// now the two element are both non zero
	if( left.factor>=right.factor)
	{
		// left element is larger than right element
		dTemp=right.value;
		for(i=0;i<(left.factor-right.factor);i++)
			dTemp /=10.0;
		dTemp +=left.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}
		
		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while(dTemp>10.0 || dTemp<-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor;
		}
		return sciTemp;
	}
	else
	{
		// right element  is larger than left element
		dTemp=left.value;
		for(i=0;i<(right.factor-left.factor);i++)
			dTemp /=10.0;
		dTemp +=right.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}
		
		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while( dTemp>10.0 || dTemp <-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor;
		}
		return sciTemp;
	}
}


// SciAddition : addition between a SCIENTIFIC NUMBERS and a real number
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
	
	sciTemp=TransToScientific(right);
	sciTemp=SciAddition(left,sciTemp);
	
	return sciTemp;
}

// SciPow : give the power of a scientific number
SCIENTIFIC_NUMBER SciPow(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
	double dTemp;
	int iTemp;
	/*
	 if(left.value==0.0 )
	 {
		printf("the base of the power is nagative\n");
		exit(1);
	 }
	 */
	if(left.value==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
		return sciTemp;
	}
	
	dTemp=(log10(fabs(left.value))+left.factor)*right;
	
	if( dTemp>0.0 )
		iTemp=int(ceil(dTemp));    //ceil(a)是求不小于a的最小整数。floor(a)表示求不大于a的最大整数
	else
		iTemp=int(floor(dTemp));
	sciTemp.factor=iTemp;
	sciTemp.value=pow(10.0,dTemp-iTemp);
	
	return sciTemp;
}



/////////
// The algorithm in SeqKmerCount: to count the frequency of each kmer in a seqence of read/genome
//	k=6
//	ACGTCACGTACGT...
//	ACGTCA index1
//	 CGTCAC index2 = floor(index1/4^5)*4 + 1(C)
//	  GTCACG index3 = floor(index2/4^5)*4 + 2(G)
//	   TCACGT index4 = floor(index3/4^5)*4 + 3(T)
//	    CACGTA index5 = floor(index4/4^5)*4 + 0(A)
//	     ACGTAC index6 = ...

// NO inv_seq 20130124

unsigned long SeqKmerCountSingle(char* seq, int k, unsigned long power)
{
	//..........................
	// int count: The length of char seq
	int count = 0;
	//..........................
	int i=0,j=0;
	unsigned long index = 0;
	unsigned long total = 0;//total number of the kmers counted in seq
	while(seq[i])
	{
		//kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
		//current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
		if(seq[i]=='A'|| seq[i] == 'a') {j++; }
		else if(seq[i]=='C'|| seq[i] == 'c') { j++; index++; }
		else if(seq[i]=='G'|| seq[i] == 'g') { j++; index+=2; }
		else if(seq[i]=='T'|| seq[i] == 't') { j++; index+=3; }
		else { j=0; index=0; }//If seq[i] is ambiguous, reset j and index
		
		if( j == k )
		{
			HashTable[index]++;
			total++;
			index %= power;// current index = floor(previous index/4^(k-1))
			j--;//the lengh of seq[i+1,i+2,...i+k-1]
		}
		index*=4;//current index = floor(previous index/4^(k-1))*4
		i++;
		count++;
	}
	
	
	return total;
}








unsigned long SeqKmerCount(char* seq, int k, unsigned long power)
{
	//..........................
	// int count: The length of char seq
	int count = 0;
	//..........................
	int i=0,j=0;
	unsigned long index = 0;
	unsigned long total = 0;//total number of the kmers counted in seq
	while(seq[i])
	{
		//kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
		//current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
		if(seq[i]=='A'|| seq[i] == 'a') {j++;}
		else if(seq[i]=='C'|| seq[i] == 'c') { j++; index++;}
		else if(seq[i]=='G'|| seq[i] == 'g') { j++; index+=2; }
		else if(seq[i]=='T'|| seq[i] == 't') { j++; index+=3; }
		else { j=0; index=0; }//If seq[i] is ambiguous, reset j and index
		
		if( j == k )
		{
			HashTable[index]++;
			total++;
			index %= power;// current index = floor(previous index/4^(k-1))
			j--;//the lengh of seq[i+1,i+2,...i+k-1]
		}
		index*=4;//current index = floor(previous index/4^(k-1))*4
		i++;
		count++;
	}
	
	//...................
	
	// Inversed seq from right to left, but not inverse A-T, G-C
	for (int m = 1; m != count + 1; m++)
		seq_inverse[m - 1] = *(seq + count - m);
	seq_inverse[count] = '\0';
	i=0,j=0;
	index = 0;
	//...................
	while(seq_inverse[i])
	{
		// Kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
		// Current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
		if(seq_inverse[i]=='T'|| seq_inverse[i] == 't') {j++;}
		else if(seq_inverse[i]=='G'|| seq_inverse[i] == 'g') { j++; index++;}
		else if(seq_inverse[i]=='C'|| seq_inverse[i] == 'c') { j++; index+=2;}
		else if(seq_inverse[i]=='A'|| seq_inverse[i] == 'a') { j++; index+=3;}
		else { j=0; index=0; }// If seq[i] is ambiguous, reset j and index
		
		if( j == k )
		{
			HashTable[index]++;
			index %= power;// Current index = floor(previous index/4^(k-1))
			j--;// The lengh of seq[i+1,i+2,...i+k-1]
		}
		index*=4;// Current index = floor(previous index/4^(k-1))*4
		i++;
	}
	
	return total;
}






// The algorithm to count Kmer for a dataset

unsigned long DataKmerCount(char* combinefileName, char* outputFileDir, char* shortName, int k, unsigned long power, bool doubleStrand, bool zeroCountOut)
{
	
	// Initialize total_DataNum. total_DataNum records the total number of kmer in rth dataset
	unsigned long total = 0;
	
	ifstream fin(combinefileName);   //EDIT fin(argv[1])
	
	//  unsigned long copynumber;
	//  cout << "begin scan" << endl;
	
	// Begin to scan lines in dataset from top to the bottom
	// countKmer ONLY for single strand

	while(fin.getline(seq,seqlength))
	{
		//cout << seq << endl;
		if (doubleStrand == false) {
			// sinlge strand
			total += SeqKmerCountSingle(seq, k, power);
			//     cout << "seq" << endl;
		}else{
			// double strand
			total += SeqKmerCount(seq, k, power);
		}

	}
	cout << "total characters: " << total << endl;
	
	fin.close();


	// Output kmer count-pw files
	char kstr[5]; sprintf(kstr, "%d", k);
	char kmerCountfile[1000];
	if (outputFileDir == NULL) {
		strcpy(kmerCountfile,shortName);
		strcat(kmerCountfile,"_k");
		strcat(kmerCountfile,kstr);
	}else{
		strcpy(kmerCountfile,outputFileDir);
		strcat(kmerCountfile,shortName);
		strcat(kmerCountfile,"_k");
		strcat(kmerCountfile,kstr);
	}
	if (doubleStrand == false) {
		// single strand
		strcat(kmerCountfile,"_ss_wc");
	}else{
		strcat(kmerCountfile,"_ds_wc");
	}
	cout << "kmerCountFileName: " << kmerCountfile << endl;
	ofstream foutKmerCount(kmerCountfile);
	
	if (zeroCountOut == false) {
		// only output those occurring
		// Sort the kmer
		vector<unsigned long> temp;
		for( unordered_map<unsigned long,unsigned long>::iterator i = HashTable.begin(); i!= HashTable.end(); i++) temp.push_back(i->first);
		sort(temp.begin(), temp.end());
		// print
		for(vector<unsigned long>::iterator j = temp.begin(); j!=temp.end(); j++)
		{
			unsigned long key = *j;
			foutKmerCount << key << "," << HashTable[key] << endl;
		}
	}else {
		// output all kmers
		for(unsigned long key = 0; key < pow(4,k); key++)
		{
			foutKmerCount << key << "," << HashTable[key] << endl;
		}

	}

	foutKmerCount.close();


	return total;
	
}








// Data Pre-process: Transform .fasta(/.fastq) to read sequence only data.
char* dataPreProcess(char *inputFileName, char *outputFileDir, char *shortName, bool fastq, bool longseq)
{
	static char combinefileName[1000];
	if (outputFileDir == NULL) {
		strcpy(combinefileName,shortName);
		strcat(combinefileName,"_combine");
	}else{
		strcpy(combinefileName,outputFileDir);
		strcat(combinefileName,shortName);
		strcat(combinefileName,"_combine");
	}
	cout << "combinefileName: " << combinefileName << endl;
	ofstream fout(combinefileName);
	
	//char seq[99999999];
	//int seqlength=99999999;
	ifstream infile(inputFileName);
	
	if(fastq == false)
	{
	  // fasta format
		int flagA = 0;
		while(infile.getline(seq,seqlength))
		{
			if (longseq == true) {
				// connect seq in diff lines
				if(seq[0]=='>') {
					if(flagA == 1){fout << endl;}
					continue;
				}else {
					flagA = 1;
					fout << seq;
				}
			}else {
				// no connect
				if(seq[0]=='>'){
					continue;
				}else{
					fout << seq << endl;
				}
			}
		}
		
	}else{
		//cout << fastq << endl;
		// fastq format
		int flagQ = 0;
		unsigned long lineNum = 1;
		while(infile.getline(seq,seqlength))
		{
			if (longseq == true) {
				//cout << "lineNum: " << lineNum << endl;
				if (lineNum % 4 == 2) {
					flagQ = 1;
					fout << seq;
					//cout << seq;
				}else {
					if(flagQ == 1){fout << endl;}
					flagQ = 0;
				}
			}else {
				if (lineNum % 4 == 2) {
					fout << seq << endl;
				}

			}

			lineNum = lineNum + 1;

		}
		
		
	}
	
	fout.close();
	
	return combinefileName;
	
}




///////////////////////////////////
/* Define the global variables */

bool zeroCountOut = false;
bool doubleStrand = false;
bool fastq = false;
bool longseq = false;
bool preprocess = true;
int k = 0;
char *inputFileName = NULL;
char *outputFileDir = NULL;
char *shortName = NULL;






int main (int argc, char **argv)
{
	/* getOptions from command line */
	int c;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options don’t set a flag.
			 We distinguish them by their indices. */
			{"zeroCountOut",     no_argument,       0, 'z'},
			{"doubleStrand",     no_argument,       0, 'd'},
			{"fastq",  no_argument,       0, 'q'},
			{"longseq",  no_argument,       0, 'l'},
			{"noPreprocess",  no_argument,       0, 'p'},
			// necessary arguments: kvalue and inputFileName !!
			{"kvalue",  required_argument, 0, 'k'},
			{"inputFileName",  required_argument, 0, 'i'},
			{"outputFileDir",  required_argument, 0, 'o'},
			{"shortName", required_argument, 0, 's'},
			{0,         0,                 0,  0 }
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		c = getopt_long (argc, argv, "zdqlpk:i:o:s:",
										 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1)
			break;
		
		switch (c)
		{
			case 'z':
				zeroCountOut = true;
				puts ("option -z: output zeroCount k-tuples. ");
				break;
				
			case 'd':
				doubleStrand = true;
				puts ("option -d: count the k-tuples for double strands. ");
				break;
				
			case 'q':
				fastq = true;
				puts ("option -q: the input file is in fastq (not fasta). ");
				break;
				
			case 'l':
				longseq = true;
				puts ("option -l: the input file is longseq (need to concatenate lines). ");
				break;
				
			case 'p':
				preprocess = false;
				puts ("option -p: no need for pre-process (combine file exists). ");
				break;

			case 'k':
				k = atoi(optarg) ;
				printf ("option -k, the value of k, with value `%s'\n", optarg);
				break;
				
			case 'i':
				inputFileName = optarg;
				printf ("option -i, the input filename, with value `%s'\n", optarg);
				break;
				
			case 'o':
				outputFileDir = optarg;
				printf ("option -o, the output directory, with value `%s'\n", optarg);
				break;
				
			case 's':
				shortName = optarg;
				printf ("option -s, the short name, with value `%s'\n", optarg);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				abort ();
		}
	}
	
	
	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
		return 0;
	}
	

	/* Set the missing options to be default. */
	if ( shortName == NULL )
	{
		std::string fullpath = inputFileName;
		int beginIdx = fullpath.rfind('/');
		std::string filename = fullpath.substr(beginIdx + 1);
		int beginIdx2 = filename.rfind('.');
		std::string filename2 = filename.substr(0, beginIdx2);
		shortName = new char[filename.length()+1];
		//strcpy(shortName, filename2.c_str()) ;
		strcpy(shortName, filename.c_str()) ;
		cout << "shortName: " << shortName << endl;

	}
	
	if (outputFileDir != NULL) {
		std::string outputPath = outputFileDir;
		if (outputPath[outputPath.length()-1] != '/') {
			strcat(outputFileDir,"/");
		}
		cout << "outputFileDir: " << outputFileDir << endl;
	
		// 0520 if outputFileDir not exist, mkdir
		//int status = mkdir (outputFileDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		char mkcmd[100];
		sprintf(mkcmd, "mkdir -p %s", outputFileDir);
		system(mkcmd);
		//if( status == -1 )
		//{
			//cout << "WARNING: outputFileDir already exists. " << endl;
		//}
	}
	
	
	
	// Parameters: k and power = 4^(k-1)
	unsigned long power = 1; for( int i = 0; i < k-1; i++) power *= 4;
	HashTable.clear();
	
	if (preprocess == true) {
		/* Preprocess the data */
		char* combinefileName = dataPreProcess(inputFileName, outputFileDir, shortName, fastq, longseq);
		
		/* Count the k-tuple */
		DataKmerCount(combinefileName, outputFileDir, shortName, k, power, doubleStrand, zeroCountOut);
		
	}else {
		/* Count the k-tuple */
		DataKmerCount(inputFileName, outputFileDir, shortName, k, power, doubleStrand, zeroCountOut);
		
	}


	

	

	
	
	
	
	
	
	
	
	
	
	
	
	return 0;
}
