
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <cstring>
//#include<cmath>
#include <vector>
#include <sstream>
#include <string>
#include <stdlib.h> 

#include <iomanip>

//strtol(s.c_str(),0,10);
using namespace std;

# include <stdio.h>
# include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>

#include<algorithm>
#include <unordered_map>
//using namespace std::tr1;

//#include<ext/hash_map>
//using __gnu_cxx::hash_map;

using namespace std;



struct SCIENTIFIC_NUMBER
{
	int factor;
	double value;
};


int ZI = 4;
//int k = 3;
//int order = 1;



struct SPECIESINFO
{
	std::string name;
	std::string dir;
	int order;
	
};






struct KMERINFO
{
	//receptacle of kmer counts
	unordered_map<unsigned long, unsigned long> HashTable;
	
	unordered_map<unsigned long, unsigned long> HashTableK_1;
	// kmer length = k - 2
	unordered_map<unsigned long, unsigned long> HashTableK_2;
	
	unordered_map<unsigned long, unsigned long> HashTableOrder;
	// kmer length = order + 1
	unordered_map<unsigned long, unsigned long> HashTableOrder_1;
	
	//receptacle of Pw(probability of a kmer word)
	unordered_map<unsigned long, SCIENTIFIC_NUMBER > HashPw;
	
	unsigned long totalKmer;
	unsigned long totalK_1;
	unsigned long totalK_2;
	unsigned long totalOrder;
	unsigned long totalOrder_1;
	
	vector<SCIENTIFIC_NUMBER> *probIIDPointer;
	
	vector< vector<SCIENTIFIC_NUMBER> > *transMatrixPointer;
	vector<SCIENTIFIC_NUMBER> *iniProbPointer;
};


//KMERINFO.p = new vector<vector> transMatrixNew(col, row)
//transMatrixNew = transmatrix
//*KMERINFO.p[colNum][rowNum]

//vector< vector<SCIENTIFIC_NUMBER> > transMatrix(transRowSize, vector<SCIENTIFIC_NUMBER>(transColSize));
//vector<SCIENTIFIC_NUMBER> iniProb;


// Scientific Number calculation funtion:

// TransToReal : trans a SCIENTIFIC NUMBER to a read number
double TransToReal(SCIENTIFIC_NUMBER dSci)
{
  double dReal=0;
  dReal = dSci.value * pow(10,dSci.factor);
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


SCIENTIFIC_NUMBER SciNegative(SCIENTIFIC_NUMBER Sci)
{
  SCIENTIFIC_NUMBER sciNeg;
  sciNeg.value = - Sci.value; sciNeg.factor = Sci.factor;
  return sciNeg;
}

SCIENTIFIC_NUMBER SciInverse(SCIENTIFIC_NUMBER Sci)
{
  SCIENTIFIC_NUMBER sciInv;
  sciInv.value = 1/Sci.value; sciInv.factor = -Sci.factor;
  return sciInv;
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


// SciMultiple : Multiplication between a SCIENTIFIC NUMBERS and a real number
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
//    cout << "this step"  << endl;
		if( left.value==0.0 )
    {
			return right;
    }else{
			return left;
    }
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





// SciPow : give the power of a scientific number
SCIENTIFIC_NUMBER SciLn(SCIENTIFIC_NUMBER sci)
{
  SCIENTIFIC_NUMBER sciTemp = TransToScientific( log( sci.value ) + sci.factor * log(10) ) ;
  
  return sciTemp;
}





void printFour(vector<int> four)
{
  cout << "print ";
  for(int it=0; it<four.size(); it++)
  {
    cout << four[it] << "," ;
  }
  cout << endl;

}



vector<int> ten2four(unsigned long ten, int k)
{
  vector<int> four (k,0);
  unsigned long tmp = ten; 
  //int currentPos = k-1;
  for(int currentPos = k-1; currentPos >=0; --currentPos)
  {
    four[currentPos]=tmp%ZI;
    tmp/=ZI; 
  }
  //while(tmp>=(ZI-1)) {four[currentPos]=tmp%ZI;tmp/=ZI;currentPos--; }
  //four[currentPos] = tmp;
  return four;
}


int four2ten(vector<int> four, int k)
{
  unsigned long ten = 0;
  for(int currentPos=(k-1); currentPos >= 0; --currentPos)
  { 
    int tmp = four[currentPos] * pow(ZI,(k-1 - currentPos));
    ten = ten + tmp;
    //cout << currentPos << " " << ten << endl;
  }
  return ten;
  
}




void loadKmerCountHash(string currentKmerFilePathName, KMERINFO* speciesKmerInfo, string kmerUsage)
{
	ifstream fin(currentKmerFilePathName.c_str());
	string currentKmerLine;
	while(getline(fin, currentKmerLine, '\n'))
	{
		unsigned long currentKmerID; unsigned long currentKmerCount;
		std::istringstream ss(currentKmerLine);
		std::string token;
		int colCount = 0;
		while(std::getline(ss, token, ',')) {
			colCount++;
			if(colCount == 1){
				currentKmerID = atoi(token.c_str());
			}else{
				currentKmerCount = atoi(token.c_str());
			}//std::cout << token << '\n';
		}
		//cout << "kmerID:" << currentKmerID << " kmerCount:" << currentKmerCount << endl;
		if(kmerUsage == "kmerCount"){
			speciesKmerInfo->HashTable[currentKmerID] = currentKmerCount;
			speciesKmerInfo->totalKmer += currentKmerCount;
		}else if(kmerUsage == "kmerCount_1"){
			speciesKmerInfo->HashTableK_1[currentKmerID] = currentKmerCount;
			speciesKmerInfo->totalK_1 += currentKmerCount;
		}else if(kmerUsage == "kmerCount_2"){
			speciesKmerInfo->HashTableK_2[currentKmerID] = currentKmerCount;
			speciesKmerInfo->totalK_2 += currentKmerCount;
		}else if(kmerUsage == "kmerOrder"){
			speciesKmerInfo->HashTableOrder[currentKmerID] = currentKmerCount;
			speciesKmerInfo->totalOrder += currentKmerCount;
		}else if(kmerUsage == "kmerOrder+1"){
			speciesKmerInfo->HashTableOrder_1[currentKmerID] = currentKmerCount;
			speciesKmerInfo->totalOrder_1 += currentKmerCount;
		}else{cout << "ERROR: wrong kmerUsage" << endl;}
		
	}
	return;
	
}





void loadSpeciesInfo(string speciesInfoFilePathName, vector<SPECIESINFO>& speciesInfoList)
{
	ifstream fin(speciesInfoFilePathName.c_str());
	string currentSpeciesLine;
	while(getline(fin, currentSpeciesLine, '\n'))
	{
		string currentName; string currentDir; int currentOrder;
		std::istringstream ss(currentSpeciesLine);
		std::string token;
		int colCount = 0;
		while(std::getline(ss, token, ' ')) {
			colCount++;
			if(colCount == 1){
				currentName = token;
			}else if(colCount == 2){
				currentDir = token;
				if (currentDir[currentDir.length()-1] != '/')
				{
					currentDir = currentDir + "/";
				}
			}else{
				currentOrder = atoi(token.c_str());
			}//std::cout << token << '\n';
		}
		SPECIESINFO currentSpeciesInfo;
		currentSpeciesInfo.name = currentName;
		currentSpeciesInfo.dir = currentDir;
		currentSpeciesInfo.order = currentOrder;
		//cout << "speciesInfo.name:" << currentSpeciesInfo.name << ", dir:" << currentSpeciesInfo.dir << ", oroder:" << currentSpeciesInfo.order << endl;
		speciesInfoList.push_back(currentSpeciesInfo);
	}
	
	return;
}




int loadTaxaInfo (string taxaFile, vector<string>& hostNCBIName, vector<string>& hostName, vector<string>& hostSuperkingdom, vector<string>& hostPhylum, vector<string>& hostClass, vector<string>& hostOrder, vector<string>& hostFamily, vector<string>& hostGenus, vector<string>& hostSpecies)
{
	ifstream taxaFileIn(taxaFile.c_str());
	string currentTaxaLine;
	string token;
	int lineNum = 0;
	int taxaColCount = 0;
	int errorTaxa = 0;
	
	while( getline(taxaFileIn, currentTaxaLine) )
	{
		//cout << currentTaxaLine << endl;
		lineNum++;
		
		if( lineNum == 1 )
		{
			// first line is the header
			std::istringstream ss(currentTaxaLine);
			taxaColCount = 0;
			while(std::getline(ss, token, '\t')) {
				taxaColCount++;
				//cout << taxaColCount << endl;
				//cout << token << endl;
				if(taxaColCount == 1){
					if(token != "hostNCBIName"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 2){
					if(token != "hostName"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 3){
					if(token != "hostSuperkingdom"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 4){
					if(token != "hostPhylum"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 5){
					if(token != "hostClass"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 6){
					if(token != "hostOrder"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 7){
					if(token != "hostFamily"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 8){
					if(token != "hostGenus"){
						errorTaxa = 1;
					}
				}else if(taxaColCount == 9){
				if(token != "hostSpecies"){
					errorTaxa = 1;
				}
			}
			}
			if(errorTaxa == 1)
			{
				cerr << "ERROR: the format of taxaFile is not correct!" << endl;
				return 0;
			}
			
		}else{
			
			std::istringstream ss(currentTaxaLine);
			taxaColCount = 0;
			while(std::getline(ss, token, '\t')) {
				taxaColCount++;
				//cout << taxaColCount << endl;
				//cout << token << endl;
				if(taxaColCount == 1){
					hostNCBIName.push_back(token);
					//cout << "check token " << hostNCBIName[0] << endl;
				}else if(taxaColCount == 2){
					hostName.push_back(token);
					//cout << "check token " << token << endl;
				}else if(taxaColCount == 3){
					hostSuperkingdom.push_back(token);
					//cout << "check token " << token << endl;
				}else if(taxaColCount == 4){
					hostPhylum.push_back(token);
				}else if(taxaColCount == 5){
					hostClass.push_back(token);
				}else if(taxaColCount == 6){
					hostOrder.push_back(token);
				}else if(taxaColCount == 7){
					hostFamily.push_back(token);
				}else if(taxaColCount == 8){
					hostGenus.push_back(token);
					//cout << "check token " << token << endl;
				}else if(taxaColCount == 9){
					hostSpecies.push_back(token);
					//cout << "check token " << token << endl;
				}
			}
		}
		
	}
	
	return lineNum-1;
}




vector<int> reverseFour(vector<int> Four)
{
  vector<int> reverseFour(Four.size(), 4);
  for(int revPos = 0; revPos < Four.size(); revPos++)
  { 
    reverseFour[revPos] = 3 - Four[Four.size()- 1 - revPos];
  }
  return reverseFour;
  
}
  


void pwIID(int ZI, int k, KMERINFO* speciesKmerInfo, vector<SCIENTIFIC_NUMBER>& probIID )
{
  // compute pw for each kmer word
  //unordered_map<unsigned long,SCIENTIFIC_NUMBER > probIID;
  //cout << HashTableOrder_1[speciesID][0]/double(totalOrder_1[speciesID]) << endl;
  for(int index = 0; index < ZI; index++)
  { 
    //cout << HashTableOrder_1[speciesID][index] << endl;
    //probIID[index] = TransToScientific(speciesKmerInfo->HashTableOrder_1[index]/double(speciesKmerInfo->totalOrder_1));
		probIID.push_back( TransToScientific(speciesKmerInfo->HashTableOrder_1[index]/double(speciesKmerInfo->totalOrder_1)) );
  }
  
  for(unsigned long ten = 0; ten < pow(ZI, k); ten++ )
  {
    SCIENTIFIC_NUMBER pw;
    pw.value = 1; pw.factor = 0;
    //vector<int> four (k, 0);
    vector<int> currentKmerFour = ten2four(ten, k);
    for(int pos = 0; pos < k; pos ++)
    {
      pw = SciMultiple(pw, probIID[currentKmerFour[pos]]);
    }
    
    speciesKmerInfo->HashPw[ten] = pw;
  }
}



void pwMC(int ZI, int k, int order, KMERINFO* speciesKmerInfo, vector<SCIENTIFIC_NUMBER>& iniProb, vector< vector<SCIENTIFIC_NUMBER> > &transMatrix)
{
	
  // compute the transition probability matrix
  int transColSize = ZI;
  int transRowSize = pow(ZI, order);
  //vector< vector<SCIENTIFIC_NUMBER> > transMatrix(transRowSize, vector<SCIENTIFIC_NUMBER>(transColSize));
  //vector<SCIENTIFIC_NUMBER> iniProb;
  for(int currentRow = 0; currentRow < transRowSize; currentRow++)
  {
    vector<int> currentRowFour = ten2four(currentRow, order);
    unsigned long countBelow = speciesKmerInfo->HashTableOrder[currentRow];
    double probBelow = countBelow/double(speciesKmerInfo->totalOrder);
    SCIENTIFIC_NUMBER probBelowSci = TransToScientific(probBelow);
    iniProb.push_back(probBelowSci);
    // %% it could be 0 ! %%
    // 20140910: the denominator cannot be 0
    SCIENTIFIC_NUMBER inv_probBelowSci;
    inv_probBelowSci.value = 0; inv_probBelowSci.factor = 0;
    if(probBelowSci.value != 0)
    {
      inv_probBelowSci = SciInverse(probBelowSci);
    }
    
    SCIENTIFIC_NUMBER rowSum;
    rowSum.value = 0; rowSum.factor = 0;
  
    for(int currentCol = 0; currentCol < transColSize; currentCol++)
    {
      //cout << "row:" << currentRow << ", col:" << currentCol << endl;
      vector<int> currentRowColFour(currentRowFour);
      currentRowColFour.push_back(currentCol); 
      for(int it=0; it< currentRowColFour.size(); it++)
      {
        //cout << currentRowColFour[it];
      }
      //cout << endl;
      int currentRowColTen = four2ten(currentRowColFour, (order+1));
      unsigned long countAbove = speciesKmerInfo->HashTableOrder_1[currentRowColTen];
      double probAbove = countAbove/double(speciesKmerInfo->totalOrder_1);
      SCIENTIFIC_NUMBER probAboveSci = TransToScientific(probAbove);

      if( probBelow != 0 )
      {
        transMatrix[currentRow][currentCol] = SciMultiple(probAboveSci,inv_probBelowSci);
      }else{
        
        transMatrix[currentRow][currentCol] = TransToScientific(0);
      }
      
      rowSum = SciAddition(rowSum, transMatrix[currentRow][currentCol]);
    
    //cout << currentRow << ", " << currentCol << ", " << countAbove << ", " << countBelow << ", " << totalOrder_1[speciesID] << ", " << totalOrder[speciesID] << ", " << countAbove/double(countBelow) << ", " << probAbove << ", " << probBelow << ", " << TransToReal(transMatrix[currentRow][currentCol]) << endl;
    //cout << totalKmer[speciesID] << ", " << totalOrder[speciesID] << ", "<< totalOrder_1[speciesID] << endl;
    }
  
    
    
    // normalize trans matrix
    for(int currentCol = 0; currentCol < transColSize; currentCol++)
    {
      transMatrix[currentRow][currentCol] = SciMultiple(transMatrix[currentRow][currentCol], SciInverse(rowSum));
      //cout << TransToReal(transMatrix[currentRow][currentCol]) << endl;
      
    }
  }
	
	
	//vector<SCIENTIFIC_NUMBER> *iniProbPointer = &iniProb;
	//cout << "iniProb[0]" << TransToReal(iniProb[0]) << endl;
	//cout << TransToReal((*iniProbPointer)[0]) << endl;
	
	//speciesKmerInfo->iniProbPointer = &iniProb;
	//cout << TransToReal((*(speciesKmerInfo->iniProbPointer))[0]) << endl;
	//speciesKmerInfo->transMatrixPointer = &transMatrix;
  
  
  // compute pw for each kmer word
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++ )
  {
    SCIENTIFIC_NUMBER pw;
    //vector<int> four (k, 0);
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    //cout << endl << endl;
    //printFour(currentKmerFour);
    vector<int> iniWordFour;
    for(int iniWordPos = 0; iniWordPos < order; iniWordPos ++)
    {
      iniWordFour.push_back(currentKmerFour[iniWordPos]);
    }
    unsigned long iniWordTen = four2ten(iniWordFour, iniWordFour.size());
    pw = iniProb[iniWordTen];
    //cout << "wordten: " << currentKmerTen << endl << "wordfour: ";
    //cout << "iniProb " << TransToReal(pw);
    
    for(int currentTransToPos = order; currentTransToPos < k; currentTransToPos++)
    {
      vector<int> currentTransFromWordFour;
      for(int transWordPos = currentTransToPos - order; transWordPos < currentTransToPos; transWordPos ++)
      {
        //cout << currentTransToPos << ", " << transWordPos << endl;
        //cout << currentKmerFour[transWordPos] << ",";
        currentTransFromWordFour.push_back(currentKmerFour[transWordPos]);
      }
      //printFour(currentTransFromWordFour);
      unsigned long currentTransFromWordTen = four2ten(currentTransFromWordFour, currentTransFromWordFour.size());
      pw = SciMultiple(transMatrix[currentTransFromWordTen][currentKmerFour[currentTransToPos]], pw);
      //cout << " trans from ";
      //printFour(currentTransFromWordFour);
      //cout << "to " << currentKmerFour[currentTransToPos] << " with prob " << TransToReal(transMatrix[currentTransFromWordTen][currentKmerFour[currentTransToPos]]) << endl;
      //cout << pw.value << ", " << pw.factor << endl;
    }
    
    speciesKmerInfo->HashPw[currentKmerTen] = pw;
    //printFour(currentKmerFour);
    //cout << " pw:" << TransToReal(pw) << endl;
  }
	
}





vector<double> D2C2computeNGS(int ZI, int k, KMERINFO* speciesKmerInfoA, KMERINFO* speciesKmerInfoB)
{
  // compute D2 statistics
  SCIENTIFIC_NUMBER D2, D2star, D2shepp;
  D2.value = 0; D2.factor = 0;  
  D2star.value = 0; D2star.factor = 0;
  D2shepp.value = 0; D2shepp.factor = 0;
  SCIENTIFIC_NUMBER C2_below[2];
  C2_below[0].value = 0; C2_below[0].factor = 0;
  C2_below[1].value = 0; C2_below[1].factor = 0;
  SCIENTIFIC_NUMBER C2star_below[2];
  C2star_below[0].value = 0; C2star_below[0].factor = 0;
  C2star_below[1].value = 0; C2star_below[1].factor = 0;
  SCIENTIFIC_NUMBER C2shepp_below[2];
  C2shepp_below[0].value = 0; C2shepp_below[0].factor = 0;
  C2shepp_below[1].value = 0; C2shepp_below[1].factor = 0;
  
  //cout << TransToReal(D2) << endl;
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
  {
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
    
    unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
    //cout << currentKmerRevTen << endl;
    //printFour(currentKmerFour);
    //printFour(currentKmerFour);
    //printFour(currentKmerRevFour);
    
    SCIENTIFIC_NUMBER p_w[2], EX_w[2], X_w[2], X_w_tilde[2], X_w_tilde_var[2];
    
    for(int speciesID=0; speciesID<2; speciesID++)
    {
			KMERINFO* speciesKmerInfo;
			if (speciesID == 0)
				speciesKmerInfo = speciesKmerInfoA;
			else
				speciesKmerInfo = speciesKmerInfoB;
			
      p_w[speciesID] = SciAddition(speciesKmerInfo->HashPw[currentKmerTen], speciesKmerInfo->HashPw[currentKmerRevTen]);
      EX_w[speciesID] = SciMultiple(p_w[speciesID], speciesKmerInfo->totalKmer);
      X_w[speciesID] = TransToScientific(speciesKmerInfo->HashTable[currentKmerTen] + speciesKmerInfo->HashTable[currentKmerRevTen]);
      //cout << "Species" << speciesID << " Ten:" << currentKmerTen << " Four: " ;
      //printFour(currentKmerFour);
      //cout << "forward:" << HashTable[speciesID][currentKmerTen] << " reverse:" << HashTable[speciesID][currentKmerRevTen] << TransToReal(X_w[speciesID]) << endl;
      X_w_tilde[speciesID] = SciAddition(X_w[speciesID], SciNegative(EX_w[speciesID]));
      
      //cout << TransToReal(p_w[speciesID]) << ", " << TransToReal(EX_w[speciesID]) << ", " << TransToReal(X_w[speciesID]) << endl;
      
      // prepare for C2 compute
      C2_below[speciesID] = SciAddition(C2_below[speciesID], SciPow(X_w[speciesID], 2));
      // prepare for C2star compute
      // 20140910: the denominator cannot be 0
      if(EX_w[speciesID].value != 0)
      {
        X_w_tilde_var[speciesID] = SciMultiple(X_w_tilde[speciesID], SciInverse(SciPow(EX_w[speciesID],0.5)));
        C2star_below[speciesID] = SciAddition(C2star_below[speciesID], SciPow(X_w_tilde_var[speciesID], 2));
      }
      
      //cout << TransToReal(X_w[speciesID]) << endl;
    }
    
    //cout << "D2:" << D2.value << ", " << D2.factor << endl;
    //if(D2.value == 0.0){cout << "D2.value = 0.0" << endl;}
    D2 = SciAddition(SciMultiple(X_w[0], X_w[1]), D2);
    //cout << "X_w[0]:" << TransToReal(X_w[1]) << " X_w[1]:" << TransToReal(X_w[1]) << " D2:" << TransToReal(D2) << endl;
    
    SCIENTIFIC_NUMBER D2star_above = SciMultiple(X_w_tilde[0], X_w_tilde[1]);
    SCIENTIFIC_NUMBER D2star_below = SciPow(SciMultiple(EX_w[0], EX_w[1]), 0.5);
    if(D2star_below.value != 0)
    {
      D2star = SciAddition(D2star, SciMultiple(D2star_above, SciInverse(D2star_below)));
    }
    
    SCIENTIFIC_NUMBER D2shepp_below = SciPow(SciAddition(SciPow(X_w_tilde[0],2), SciPow(X_w_tilde[1],2)),0.5);
    // 20140910: the denominator cannot be 0
    if(D2shepp_below.value != 0)
    {
      D2shepp = SciAddition(D2shepp, SciMultiple(D2star_above, SciInverse(D2shepp_below)));
      // prepare for c2shepp compute
      C2shepp_below[0] = SciAddition(C2shepp_below[0], SciMultiple(SciPow(X_w_tilde[0],2), SciInverse(D2shepp_below)));
      C2shepp_below[1] = SciAddition(C2shepp_below[1], SciMultiple(SciPow(X_w_tilde[1],2), SciInverse(D2shepp_below)));
    }
    
  }
  
  SCIENTIFIC_NUMBER C2 = SciMultiple(D2, SciInverse(SciMultiple(SciPow(C2_below[0],0.5),SciPow(C2_below[1],0.5))));
  
  SCIENTIFIC_NUMBER C2star = SciMultiple(D2star, SciInverse( SciMultiple(SciPow(C2star_below[0], 0.5), SciPow(C2star_below[1], 0.5)) ));
  
  SCIENTIFIC_NUMBER C2shepp = SciMultiple(D2shepp, SciInverse( SciMultiple(SciPow(C2shepp_below[0],0.5),SciPow(C2shepp_below[1],0.5) )) );
  //cout << "d2shepp_below " << TransToReal(SciMultiple(SciPow(C2shepp_below[0],0.5),SciPow(C2shepp_below[1],0.5) )) << endl;
  

  vector<double> C2values;
  C2values.push_back(0.5*(1-TransToReal(C2)));
  C2values.push_back(0.5*(1-TransToReal(C2star)));
  C2values.push_back(0.5*(1-TransToReal(C2shepp)));
  //cout << D2C2values[0] << endl;

  return C2values;
}




vector<double> EuMaChCombineDistNGS(int ZI, int k, KMERINFO* speciesKmerInfoA, KMERINFO* speciesKmerInfoB)
{
  //cout << TransToReal(D2) << endl;
  SCIENTIFIC_NUMBER EuDist ; EuDist.value=0; EuDist.factor=0;
  SCIENTIFIC_NUMBER MaDist ; MaDist.value=0; MaDist.factor=0;
  SCIENTIFIC_NUMBER ChDist ; ChDist.value=0; ChDist.factor=0;
	
	SCIENTIFIC_NUMBER EuUDist ; EuUDist.value=0; EuUDist.factor=0;
	SCIENTIFIC_NUMBER MaUDist ; MaUDist.value=0; MaUDist.factor=0;
	SCIENTIFIC_NUMBER ChUDist ; ChUDist.value=0; ChUDist.factor=0;
	
	SCIENTIFIC_NUMBER EuFDist ; EuFDist.value=0; EuFDist.factor=0;
  
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
  {
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
    unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
    
    // frequency
    SCIENTIFIC_NUMBER X1 = TransToScientific( ( speciesKmerInfoA->HashTable[currentKmerTen] + speciesKmerInfoA->HashTable[currentKmerRevTen] ) / double ( 2 * speciesKmerInfoA->totalKmer ) );
    SCIENTIFIC_NUMBER X2 = TransToScientific( ( speciesKmerInfoB->HashTable[currentKmerTen] + speciesKmerInfoB->HashTable[currentKmerRevTen] ) / double ( 2 * speciesKmerInfoB->totalKmer ) );
    
    SCIENTIFIC_NUMBER diff = SciAddition(X1, SciNegative(X2)) ;
    diff.value = fabs(diff.value);
    
    EuDist = SciAddition( EuDist, SciPow( diff , 2 ) );
    
    MaDist = SciAddition( MaDist, diff ) ;
    
    if( TransToReal(diff) > TransToReal(ChDist) )
    {
      ChDist = diff;
    }
		
		// EuF
		SCIENTIFIC_NUMBER pw1 = SciAddition( speciesKmerInfoA->HashPw[currentKmerTen], speciesKmerInfoA->HashPw[currentKmerRevTen]);
		SCIENTIFIC_NUMBER Fw1;
		Fw1.value = 0; Fw1.factor = 0;
		if (pw1.value != 0)
		{
			Fw1 = SciMultiple( X1, SciInverse( SciMultiple(pw1, 0.5)) );
		}
		
		SCIENTIFIC_NUMBER pw2 = SciAddition( speciesKmerInfoB->HashPw[currentKmerTen], speciesKmerInfoB->HashPw[currentKmerRevTen] );
		SCIENTIFIC_NUMBER Fw2;
		Fw2.value = 0; Fw2.factor = 0;
		if (pw2.value != 0)
		{
			Fw2 = SciMultiple( X2, SciInverse( SciMultiple(pw2, 0.5) ));
		}
		
		SCIENTIFIC_NUMBER diffF = SciAddition(Fw1, SciNegative(Fw2)) ;
		
		EuFDist = SciAddition( EuFDist, SciPow( diffF , 2 ) );
		
  }
	//cout << "count: " << count << endl;

  vector<double> EuMaChCombineDistvalues;
  EuMaChCombineDistvalues.push_back(TransToReal(SciPow(EuDist, 0.5)));
  EuMaChCombineDistvalues.push_back(TransToReal(MaDist));
  EuMaChCombineDistvalues.push_back(TransToReal(ChDist));
	
	EuMaChCombineDistvalues.push_back(TransToReal(SciPow(EuFDist, 0.5)));

  return EuMaChCombineDistvalues;
  
}






// Willner et al. Di, Tri, Tetra
double WillnerDiNGS(int ZI, int k, KMERINFO* speciesKmerInfoA, KMERINFO* speciesKmerInfoB )
{
  
  double deltaDi = 0;
  //cout << "test" << ", k" << k << endl;

  // Di
  if( k == 2 )
  {

    for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
    {
      //cout << "test" << endl;
      vector<int> currentKmerFour = ten2four(currentKmerTen, k);
      vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
      unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
      
  //    cout << "test" << endl;
      
      // frequency
      double fab1 = ( speciesKmerInfoA->HashTable[currentKmerTen] + speciesKmerInfoA->HashTable[currentKmerRevTen] ) / double( 2 * speciesKmerInfoA->totalKmer );
      double fab2 = ( speciesKmerInfoB->HashTable[currentKmerTen] + speciesKmerInfoB->HashTable[currentKmerRevTen] ) / double( 2 * speciesKmerInfoB->totalKmer );
      
      //cout << "fab1" << fab1 << endl;
      
      double fa1 = ( speciesKmerInfoA->HashTableOrder_1[currentKmerFour[0]] + speciesKmerInfoA->HashTableOrder_1[3-currentKmerFour[0]] ) / double( 2 * speciesKmerInfoA->totalOrder_1 );
      double fb1 = ( speciesKmerInfoA->HashTableOrder_1[currentKmerFour[1]] + speciesKmerInfoA->HashTableOrder_1[3-currentKmerFour[1]] ) / double( 2 * speciesKmerInfoA->totalOrder_1 );
      
      double fa2 = ( speciesKmerInfoB->HashTableOrder_1[currentKmerFour[0]] + speciesKmerInfoB->HashTableOrder_1[3-currentKmerFour[0]] ) / double( 2 * speciesKmerInfoB->totalOrder_1 );
      double fb2 = ( speciesKmerInfoB->HashTableOrder_1[currentKmerFour[1]] + speciesKmerInfoB->HashTableOrder_1[3-currentKmerFour[1]] ) / double( 2 * speciesKmerInfoB->totalOrder_1 );
      
      deltaDi = deltaDi + fabs( fab1 / ( fa1 * fb1) - fab2 / ( fa2 * fb2) );
      
    } 
    return deltaDi;
  }
  
  double gammaTri = 0;
  if( k == 3 )
  { 
    for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
    {
      vector<int> currentKmerFour = ten2four(currentKmerTen, k);
      vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
      unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
      
      // frequency
      double fabc1 = ( speciesKmerInfoA->HashTable[currentKmerTen] + speciesKmerInfoA->HashTable[currentKmerRevTen] ) / double( 2 * speciesKmerInfoA->totalKmer ) ;
      double fabc2 = ( speciesKmerInfoB->HashTable[currentKmerTen] + speciesKmerInfoB->HashTable[currentKmerRevTen] ) / double( 2 * speciesKmerInfoB->totalKmer ) ;
      
      vector<double> fTwoLetters1; // 3 in total
      vector<double> fTwoLetters2; // 3 in total
      for( int missPos = 0; missPos < 4; missPos++)
      {
        double wordFreq1=0;
        double wordFreq2=0;
        for( int missLetter = 0; missLetter < 4; missLetter++)
        {
          vector<int> currentWordFour = currentKmerFour;
          currentWordFour[missPos] = missLetter;
          unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
          unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
          wordFreq1 = wordFreq1 + ( speciesKmerInfoA->HashTable[currentWordTen] + speciesKmerInfoA->HashTable[currentWordRevTen] ) ;
          wordFreq2 = wordFreq2 + ( speciesKmerInfoB->HashTable[currentWordTen] + speciesKmerInfoB->HashTable[currentWordRevTen] ) ;
        }
        fTwoLetters1.push_back( wordFreq1 / double (2 * speciesKmerInfoA->totalKmer ) ) ;
        fTwoLetters2.push_back( wordFreq2 / double (2 * speciesKmerInfoB->totalKmer ) ) ;

      }
      
      vector<double> fOneLetters1; // 3 in total
      vector<double> fOneLetters2; // 3 in total
      for( int missPosFirst = 0; missPosFirst < 4; missPosFirst++)
      {
        for( int missPosSecond = missPosFirst + 1; missPosSecond < 4; missPosSecond++)
        {
          double wordFreq1=0;
          double wordFreq2=0;
          for( int missLetterFirst = 0; missLetterFirst < 4; missLetterFirst++)
          {
            for( int missLetterSecond = 0; missLetterSecond < 4; missLetterSecond++)
            {
              //cout << "missPosFirst," << missPosFirst << ",missPosSecond," << missPosSecond << endl;
              //cout << "missLetterFirst," << missLetterFirst << ",missLetterSecond," << missLetterSecond << endl; 
              vector<int> currentWordFour = currentKmerFour;
              currentWordFour[missPosFirst] = missLetterFirst;
              currentWordFour[missPosSecond] = missLetterSecond;
              unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
              unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
              wordFreq1 = wordFreq1 + ( speciesKmerInfoA->HashTable[currentWordTen] + speciesKmerInfoA->HashTable[currentWordRevTen] ) ;
              wordFreq2 = wordFreq2 + ( speciesKmerInfoB->HashTable[currentWordTen] + speciesKmerInfoB->HashTable[currentWordRevTen] ) ;
              //cout << "missPosFirst," << missPosFirst << ",missPosSecond," << missPosSecond << endl; 
              //printFour(currentKmerFour);
              //cout <<  "wordFreq," << wordFreq1 << endl;
            }
          }
          fOneLetters1.push_back( wordFreq1 / double (2 * speciesKmerInfoA->totalKmer ) ) ;
          fOneLetters2.push_back( wordFreq2 / double (2 * speciesKmerInfoB->totalKmer ) ) ;
        }
      }
      //cout << fOneLetters1[0] << "," << fOneLetters1[1] << "," << fOneLetters1[2] << endl; 
      //cout << fOneLetters1[0] << "," << fOneLetters1[1] << "," << fOneLetters1[2] << endl;
			double gamma1 = 0;
			double gamma2 = 0;
			if( fTwoLetters1[0] != 0 && fTwoLetters1[1] != 0 && fTwoLetters1[2] != 0 )
			{
				gamma1 = fabc1 * fOneLetters1[0] * fOneLetters1[1] * fOneLetters1[2] / ( fTwoLetters1[0] * fTwoLetters1[1] * fTwoLetters1[2] ) ;
			}
			
			if( fTwoLetters2[0] != 0 && fTwoLetters2[1] != 0 && fTwoLetters2[2] != 0 )
			{
				gamma2 = fabc2 * fOneLetters2[0] * fOneLetters2[1] * fOneLetters2[2] / ( fTwoLetters2[0] * fTwoLetters2[1] * fTwoLetters2[2] ) ;
			}
      gammaTri = gammaTri + fabs( gamma1 - gamma2) ;
      //cout << gamma1 << ","  << gamma2 << endl;
      //cout << fabs(gamma1-gamma2) << ".." << gammaTri << endl;
    }
    //cout << a << "," << gammaTri << endl;
    return gammaTri;
  }
  
  double tauTetra = 0;
  if( k == 4)
  {
    for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
    { 
      vector<int> currentKmerFour = ten2four(currentKmerTen, k);
      vector<int> currentKmerRevFour = reverseFour(currentKmerFour);
      unsigned long currentKmerRevTen = four2ten(currentKmerRevFour, currentKmerFour.size());
      
      double fabcd1 = ( speciesKmerInfoA->HashTable[currentKmerTen] + speciesKmerInfoA->HashTable[currentKmerRevTen] ) / double( 2 * speciesKmerInfoA->totalKmer ) ;
      double fabcd2 = ( speciesKmerInfoB->HashTable[currentKmerTen] + speciesKmerInfoB->HashTable[currentKmerRevTen] ) / double( 2 * speciesKmerInfoB->totalKmer ) ;
      
      vector<double> fThreeLetters1; // 4 in total
      vector<double> fThreeLetters2; // 4 in total
      for( int missPos = 0; missPos < 4; missPos++)
      {
        double wordFreq1=0;
        double wordFreq2=0;
        for( int missLetter = 0; missLetter < 4; missLetter++)
        {
          vector<int> currentWordFour = currentKmerFour;
          currentWordFour[missPos] = missLetter;
          unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
          unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
          wordFreq1 = wordFreq1 + ( speciesKmerInfoA->HashTable[currentWordTen] + speciesKmerInfoA->HashTable[currentWordRevTen] ) ;
          wordFreq2 = wordFreq2 + ( speciesKmerInfoB->HashTable[currentWordTen] + speciesKmerInfoB->HashTable[currentWordRevTen] ) ;
        }
        fThreeLetters1.push_back( wordFreq1 / double (2 * speciesKmerInfoA->totalKmer ) ) ;
        fThreeLetters2.push_back( wordFreq2 / double (2 * speciesKmerInfoB->totalKmer ) ) ;
        
      }
      
      vector<double> fTwoLetters1; // 6 in total
      vector<double> fTwoLetters2; // 6 in total
      for( int missPosFirst = 0; missPosFirst < 4; missPosFirst++)
      {
        for( int missPosSecond = missPosFirst + 1; missPosSecond < 4; missPosSecond++)
        {
          double wordFreq1=0;
          double wordFreq2=0;
          for( int missLetterFirst = 0; missLetterFirst < 4; missLetterFirst++)
          {
            for( int missLetterSecond = 0; missLetterSecond < 4; missLetterSecond++)
            {
              vector<int> currentWordFour = currentKmerFour;
              currentWordFour[missPosFirst] = missLetterFirst;
              currentWordFour[missPosSecond] = missLetterSecond;
              unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
              unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
              wordFreq1 = wordFreq1 + ( speciesKmerInfoA->HashTable[currentWordTen] + speciesKmerInfoA->HashTable[currentWordRevTen] ) ;
              wordFreq2 = wordFreq2 + ( speciesKmerInfoB->HashTable[currentWordTen] + speciesKmerInfoB->HashTable[currentWordRevTen] ) ;
            }
          }
          fTwoLetters1.push_back( wordFreq1 / double (2 * speciesKmerInfoA->totalKmer ) ) ;
          fTwoLetters2.push_back( wordFreq2 / double (2 * speciesKmerInfoB->totalKmer ) ) ;
        }
      }

      vector<double> fOneLetters1; // 4 in total
      vector<double> fOneLetters2; // 4 in total
      for( int missPosFirst = 0; missPosFirst < 4; missPosFirst++)
      {
        for( int missPosSecond = missPosFirst + 1; missPosSecond < 4; missPosSecond++)
        {
          for( int missPosThird = missPosSecond + 1; missPosThird < 4; missPosThird++)
          {
            double wordFreq1=0;
            double wordFreq2=0;
            for( int missLetterFirst = 0; missLetterFirst < 4; missLetterFirst++)
            {
              for( int missLetterSecond = 0; missLetterSecond < 4; missLetterSecond++)
              {
                for( int missLetterThird = 0; missLetterThird < 4; missLetterThird++)
                {
                  vector<int> currentWordFour = currentKmerFour;
                  currentWordFour[missPosFirst] = missLetterFirst;
                  currentWordFour[missPosSecond] = missLetterSecond;
                  currentWordFour[missPosThird] = missLetterThird;
                  unsigned long currentWordTen = four2ten(currentWordFour, currentWordFour.size());
                  unsigned long currentWordRevTen = four2ten(reverseFour(currentWordFour), currentWordFour.size());
                  wordFreq1 = wordFreq1 + ( speciesKmerInfoA->HashTable[currentWordTen] + speciesKmerInfoA->HashTable[currentWordRevTen] ) ;
                  wordFreq2 = wordFreq2 + ( speciesKmerInfoB->HashTable[currentWordTen] + speciesKmerInfoB->HashTable[currentWordRevTen] ) ;
                }
              }
            }
            
            fOneLetters1.push_back( wordFreq1 / double (2 * speciesKmerInfoA->totalKmer ) ) ;
            fOneLetters2.push_back( wordFreq2 / double (2 * speciesKmerInfoB->totalKmer ) ) ;
            
            //cout << missPosFirst << "," << missPosSecond << "," << missPosThird << endl;
            //cout << fOneLetters1.size() << "," << fOneLetters1[fOneLetters1.size()-1] << endl;
          }
        }
      }
      //cout << "fTwoLetters1_0," << fTwoLetters1[0] << endl;
      //cout << fOneLetters1[1] << endl;
      //cout << fTwoLetters1[1] << endl;
      //cout << fThreeLetters1[1] << endl;
			double tau1 = 0;
			double tau2 = 0;
			if( fThreeLetters1[0] != 0 && fThreeLetters1[1] != 0 && fThreeLetters1[2] != 0 && fThreeLetters1[3] != 0 )
			{
				tau1 = fabcd1 * fTwoLetters1[0] * fTwoLetters1[1] * fTwoLetters1[2] * fTwoLetters1[3] * fTwoLetters1[4] * fTwoLetters1[5] / ( fThreeLetters1[0] * fThreeLetters1[1] * fThreeLetters1[2] * fThreeLetters1[3] * fOneLetters1[0] * fOneLetters1[1] * fOneLetters1[2] * fOneLetters1[3]) ;
			}
			
			if( fThreeLetters2[0] != 0 && fThreeLetters2[1] != 0 && fThreeLetters2[2] != 0 && fThreeLetters2[3] != 0 )
			{
				tau2 = fabcd2 * fTwoLetters2[0] * fTwoLetters2[1] * fTwoLetters2[2] * fTwoLetters2[3] * fTwoLetters2[4] * fTwoLetters2[5] / ( fThreeLetters2[0] * fThreeLetters2[1] * fThreeLetters2[2] * fThreeLetters2[3] * fOneLetters2[0] * fOneLetters2[1] * fOneLetters2[2] * fOneLetters2[3]) ;
			}
      
      //cout << "tau1," << tau1 << endl;
      
      tauTetra = tauTetra + fabs(tau1 - tau2) ;
    }

    return tauTetra ;
  }
  
  // If k is not 2, 3, or 4
  return 0.0;
}



/// need to write it as complimentary chains

// not consider complimentary chains
// the version of Teeling on complimentary chains is not exactly right
/// need to write it as complimentary chains
vector<double> HAOTeelingcompute(int ZI, int k, KMERINFO* speciesKmerInfoA, KMERINFO* speciesKmerInfoB)
{
  SCIENTIFIC_NUMBER HAO_above;
  HAO_above.value = 0; HAO_above.factor = 0;
  SCIENTIFIC_NUMBER HAO_below[2]; 
  HAO_below[0].value = 0; HAO_below[0].factor = 0; 
  HAO_below[1].value = 0; HAO_below[1].factor = 0;
	
	SCIENTIFIC_NUMBER Teeling_above;
	Teeling_above.value = 0; Teeling_above.factor = 0;
	SCIENTIFIC_NUMBER Teeling_below[2];
	Teeling_below[0].value = 0; Teeling_below[0].factor = 0;
	Teeling_below[1].value = 0; Teeling_below[1].factor = 0;
	
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++ )
  {
    //vector<int> four (k, 0);
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    //cout << endl << endl;
    //printFour(currentKmerFour);
    vector<int> currentKmerRemoveLastFour;
    vector<int> currentKmerRemoveFirstFour;
    vector<int> currentKmerRemoveTwoFour;
    for(int position = 0; position < k; position++)
    {
      if( position != k)
      {
        currentKmerRemoveLastFour.push_back(currentKmerFour[position]);
      }
      if( position != 0)
      {
        currentKmerRemoveFirstFour.push_back(currentKmerFour[position]);
      }
      if( position != 0 && position != k)
      {
        currentKmerRemoveTwoFour.push_back(currentKmerFour[position]);
      }
    }
    unsigned long currentKmerRemoveLastTen = four2ten(currentKmerRemoveLastFour, (k-1));
    unsigned long currentKmerRemoveFirstTen = four2ten(currentKmerRemoveFirstFour, (k-1));
    unsigned long currentKmerRemoveTwoTen = four2ten(currentKmerRemoveTwoFour, (k-2));
    
    SCIENTIFIC_NUMBER avalue[2]; 
    avalue[0].value = 0; avalue[0].factor = 0; 
    avalue[1].value = 0; avalue[1].factor = 0;
		
		SCIENTIFIC_NUMBER Zvalue[2];
		Zvalue[0].value = 0; Zvalue[0].factor = 0;
		Zvalue[1].value = 0; Zvalue[1].factor = 0;
		
    for(int speciesID=0; speciesID<2; speciesID++)
    {
			KMERINFO* speciesKmerInfo;
			if (speciesID == 0)
				speciesKmerInfo = speciesKmerInfoA;
			else
				speciesKmerInfo = speciesKmerInfoB;
			
      SCIENTIFIC_NUMBER pw0; pw0.value = 0; pw0.factor = 0;
			SCIENTIFIC_NUMBER Ew; Ew.value = 0; Ew.factor = 0;
			SCIENTIFIC_NUMBER Vw; Vw.value = 0; Vw.factor = 0;
      if(speciesKmerInfo->HashTableK_2[currentKmerRemoveTwoTen] != 0)
      {
				// hao
        //cout << "HashTableK_1 " << "species " << speciesID << " word " << currentKmerRemoveLastTen << " kmercount " << speciesKmerInfo->HashTableK_1[currentKmerRemoveLastTen] << " totalK_1 " << speciesKmerInfo->totalK_1 << "; " << speciesKmerInfo->HashTableK_1[currentKmerRemoveFirstTen] << " totalK_1 " << speciesKmerInfo->totalK_1 << "; " << speciesKmerInfo->HashTableK_2[currentKmerRemoveTwoTen] << " totalK_2 " << speciesKmerInfo->totalK_2 << endl;
        SCIENTIFIC_NUMBER pw0_up1 = TransToScientific(speciesKmerInfo->HashTableK_1[currentKmerRemoveLastTen]/double(speciesKmerInfo->totalK_1));
        SCIENTIFIC_NUMBER pw0_up2 = TransToScientific(speciesKmerInfo->HashTableK_1[currentKmerRemoveFirstTen]/double(speciesKmerInfo->totalK_1));
        SCIENTIFIC_NUMBER pw0_down = TransToScientific(speciesKmerInfo->HashTableK_2[currentKmerRemoveTwoTen]/double(speciesKmerInfo->totalK_2));
        
        pw0 = SciMultiple(SciMultiple(pw0_up1, pw0_up2),SciInverse(pw0_down));
        //cout << HashTableK_1[speciesID][currentKmerRemoveLastTen]/totalK_1[speciesID] << " " << TransToReal(pw0_up1) << " " << TransToReal(pw0_up2) << " " << TransToReal(pw0_down) << " " << TransToReal(pw0) << endl;
        //fw0 = totalKmer[speciesID] * pw0;
				
				// teeling
				SCIENTIFIC_NUMBER Ew_up1 = TransToScientific(speciesKmerInfo->HashTableK_1[currentKmerRemoveLastTen]);
				SCIENTIFIC_NUMBER Ew_up2 = TransToScientific(speciesKmerInfo->HashTableK_1[currentKmerRemoveFirstTen]);
				SCIENTIFIC_NUMBER Ew_down = TransToScientific(speciesKmerInfo->HashTableK_2[currentKmerRemoveTwoTen]);
				Ew = SciMultiple(SciMultiple(Ew_up1, Ew_up2),SciInverse(Ew_down));
				
				SCIENTIFIC_NUMBER one; one.value = 1; one.factor = 0;
				SCIENTIFIC_NUMBER Vw_1 = SciAddition(one, SciNegative(SciMultiple(Ew_up1, SciInverse(Ew_down))) );
				SCIENTIFIC_NUMBER Vw_2 = SciAddition(one, SciNegative(SciMultiple(Ew_up2, SciInverse(Ew_down))) );
				Vw = SciMultiple( SciMultiple(Ew, Vw_1), Vw_2 );
				
      }
      // hao
      SCIENTIFIC_NUMBER pw = TransToScientific(speciesKmerInfo->HashTable[currentKmerTen]/double(speciesKmerInfo->totalKmer));
      //cout << "pw " << TransToReal(pw) << " pw0 " << TransToReal(pw0) << endl;
      if(pw0.value != 0)
      {
        avalue[speciesID] = SciMultiple(SciAddition(pw, SciNegative(pw0)), SciInverse(pw0));
      }
      //cout << TransToReal(pw) << " " << TransToReal(pw0) << endl;
      //cout << TransToReal(avalue[speciesID]) << endl;
      HAO_below[speciesID] = SciAddition(HAO_below[speciesID], SciPow(avalue[speciesID],2));
			
			// teeling
			SCIENTIFIC_NUMBER Nw = TransToScientific(speciesKmerInfo->HashTable[currentKmerTen]);
			//cout << "Nw " << TransToReal(Nw) << " Ew " << TransToReal(Ew) << " Vw " << TransToReal(SciPow(Vw, 0.5)) << endl;
			if(Vw.value != 0)
			{
				Zvalue[speciesID] = SciMultiple(SciAddition(Nw, SciNegative(Ew)), SciInverse(SciPow(Vw, 0.5)) );
			}
			//cout << TransToReal(pw) << " " << TransToReal(pw0) << endl;
			//cout << TransToReal(avalue[speciesID]) << endl;
			Teeling_below[speciesID] = SciAddition(Teeling_below[speciesID], SciPow(Zvalue[speciesID],2));
    }
    
    HAO_above = SciAddition(HAO_above, SciMultiple(avalue[0], avalue[1]));
		Teeling_above = SciAddition(Teeling_above, SciMultiple(Zvalue[0], Zvalue[1]));
    
  }
  //cout << TransToReal(HAO_above) << " " << TransToReal(HAO_below[0]) << " " << TransToReal(HAO_below[1]) << endl;
  
  SCIENTIFIC_NUMBER HAOvalueSCI = SciMultiple(HAO_above, SciInverse(SciPow(SciMultiple(HAO_below[0], HAO_below[1]) , 0.5)));                                                                    
  double HAOvalue = TransToReal(HAOvalueSCI);
	
	SCIENTIFIC_NUMBER TeelingValueSCI = SciMultiple(Teeling_above, SciInverse(SciPow(SciMultiple(Teeling_below[0], Teeling_below[1]) , 0.5)));
	double TeelingValue = TransToReal(TeelingValueSCI);
  //cout << HAOvalue << endl;
	
	vector<double> HaoTeelingvalues;
	HaoTeelingvalues.push_back(0.5*(1-HAOvalue));
	HaoTeelingvalues.push_back(0.5*(1-TeelingValue));
  return HaoTeelingvalues;
}




double JScomputeIID(vector<SCIENTIFIC_NUMBER>& iniProb)
{
  int size = iniProb.size();
  SCIENTIFIC_NUMBER entropy;
  entropy.value = 0; entropy.factor = 0;
  for(int currentIndex = 0; currentIndex < size; currentIndex++)
  {
    // compute JS distance
    if( TransToReal(iniProb[currentIndex]) != 0 ){
      
      SCIENTIFIC_NUMBER logTransProbSci = TransToScientific(log2(TransToReal(iniProb[currentIndex])));
      entropy = SciAddition(entropy, SciMultiple(iniProb[currentIndex], logTransProbSci));
    }else{
      //cout << "trans prob = 0, no addition on the entropy." << endl;
    }
  }
  return -TransToReal(entropy);
}






double JScomputeMC(vector<SCIENTIFIC_NUMBER>& iniProb, vector< vector<SCIENTIFIC_NUMBER> >& transMatrix)
{
  int transRowSize = iniProb.size();
  int transColSize = transMatrix[0].size();
  SCIENTIFIC_NUMBER entropyOverRow;
  entropyOverRow.value = 0; entropyOverRow.factor = 0;
  for(int currentRow = 0; currentRow < transRowSize; currentRow++)
  {
    // normalize trans matrix
    SCIENTIFIC_NUMBER entropyOverCol;
    entropyOverCol.value = 0; entropyOverCol.factor = 0;
    for(int currentCol = 0; currentCol < transColSize; currentCol++)
    {
      // compute JS distance
      if( TransToReal(transMatrix[currentRow][currentCol]) != 0 ){
        
        SCIENTIFIC_NUMBER logTransProbSci = TransToScientific(log2(TransToReal(transMatrix[currentRow][currentCol])));
        entropyOverCol = SciAddition(entropyOverCol, SciMultiple(transMatrix[currentRow][currentCol], logTransProbSci));
      }else{
        //cout << "trans prob = 0, no addition on the entropy." << endl;
      }
    }
    
    entropyOverRow = SciAddition(entropyOverRow, SciMultiple(iniProb[currentRow], entropyOverCol));
    
  }
  return -TransToReal(entropyOverRow);
  
}






// not necessary to make it to include complimentary
double S2compute(int ZI, int k, KMERINFO* speciesKmerInfoA, KMERINFO* speciesKmerInfoB)
{

  //cout << TransToReal(D2) << endl;
  SCIENTIFIC_NUMBER S2;
  S2.factor = 0; S2.value = 0;
  
  for(unsigned long currentKmerTen = 0; currentKmerTen < pow(ZI, k); currentKmerTen++)
  {
    vector<int> currentKmerFour = ten2four(currentKmerTen, k);
    //cout << currentKmerRevTen << endl;
    //printFour(currentKmerFour);
    
    SCIENTIFIC_NUMBER phi[2];
    
    for(int speciesID=0; speciesID<2; speciesID++)
    {
			KMERINFO* speciesKmerInfo;
			if (speciesID == 0)
				speciesKmerInfo = speciesKmerInfoA;
			else
				speciesKmerInfo = speciesKmerInfoB;
			
      phi[speciesID] = SciMultiple(TransToScientific( speciesKmerInfo->HashTable[currentKmerTen] / double(speciesKmerInfo->totalKmer) ), speciesKmerInfo->HashPw[currentKmerTen] ) ;
      //cout << HashTable[speciesID][currentKmerTen] << "," << totalKmer[speciesID] << ", pw: " << TransToReal(HashPw[speciesID][currentKmerTen]) << endl;
      
      //cout << pow(ZI, k) << ", " << currentKmerTen << ", phi: " << TransToReal(phi[speciesID]) << endl;

    }

    
    if( phi[0].value != 0 )
    {
      SCIENTIFIC_NUMBER tmp1 = SciLn(SciMultiple(phi[0], 2)) ;
      
      SCIENTIFIC_NUMBER tmp2 = SciLn(SciAddition(phi[0], phi[1])) ;
      
      SCIENTIFIC_NUMBER tmp3 = SciMultiple(phi[0], SciAddition( tmp1, SciNegative(tmp2))) ;
      tmp3.value = fabs(tmp3.value) ;
      
      //cout << "tmp1: " << TransToReal(tmp1) << ", tmp2: " << TransToReal(tmp2) << ", tmp1-tmp2: " << TransToReal(SciAddition( tmp1, SciNegative(tmp2))) << ", tmp3: " << TransToReal(tmp3) << endl;
      
      S2 = SciAddition( S2, tmp3 ) ;
      
    }
    
    if( phi[1].value != 0 )
    {
      
      SCIENTIFIC_NUMBER tmp4 = SciLn(SciMultiple(phi[1], 2)) ;
      SCIENTIFIC_NUMBER tmp2 = SciLn(SciAddition(phi[0], phi[1])) ;
      SCIENTIFIC_NUMBER tmp5 = SciMultiple(phi[1], SciAddition( tmp4, SciNegative(tmp2) )) ;
      tmp5.value = fabs(tmp5.value) ;
      
      S2 = SciAddition( S2, tmp5 ) ;
      
    }
    
    //cout << pow(ZI, k) << ", " << currentKmerTen << ", " << TransToReal(phi[0]) << ", " << TransToReal(phi[1]) << ", " << TransToReal(S2) << endl;
  }
  
  // so far S2 is very small, after dividing by pow(ZI, k), this term is almost 0!!!
  // so S2 is almost always 2*ln2!!! I dont get it!
  
  S2 = SciMultiple( S2, SciInverse(TransToScientific(pow(ZI, k))) ) ;
  S2 = SciAddition( S2, TransToScientific(2*log(2)) ) ;
  
  return TransToReal(S2);

}




void computeMultiStatsReturn(int ZI, int k, KMERINFO* speciesKmerInfoA, KMERINFO* speciesKmerInfoB, SPECIESINFO speciesInfoA, SPECIESINFO speciesInfoB, vector<string>& measureNames, vector<double>& measureValues)
{
	
	// 5. compute Eu, Ma, Ch
	//cout << "== compute the EuMaCh distances == " << endl;
	vector<double> EuMaChCombineDistvalues = EuMaChCombineDistNGS(ZI, k, speciesKmerInfoA, speciesKmerInfoB) ;
	//cout << "Eu, " << EuMaChCombineDistvalues[0] << endl ;
	//cout << "Ma, " << EuMaChCombineDistvalues[1] << endl ;
	//cout << "Ch, " << EuMaChCombineDistvalues[2] << endl ;
	//cout << "EuF, " << EuMaChCombineDistvalues[3] << endl ;
	
	measureNames.push_back("Eu");
	measureValues.push_back(EuMaChCombineDistvalues[0]);
	measureNames.push_back("Ma");
	measureValues.push_back(EuMaChCombineDistvalues[1]);
	measureNames.push_back("Ch");
	measureValues.push_back(EuMaChCombineDistvalues[2]);
	measureNames.push_back("EuF");
	measureValues.push_back(EuMaChCombineDistvalues[3]);
	
	//fout << "EuF, " << EuFDistvalues[0] << endl;
	
	// 4. compute JS distance
	//cout << "== compute the JSdistance == " << endl;
	double entropyRate1;
	double entropyRate2;
	double entropyRateAve;
	double JSdivergence;
	double JSdistance;
	int orderA = speciesInfoA.order;
	int orderB = speciesInfoB.order;
	if( orderA != orderB )
	{
		//cout << "ERROR: order 1 MUST equal order 2." << endl;
	}else{
		
		if(orderA == 0){
			
			vector<SCIENTIFIC_NUMBER> iidProbAve;
			for(int currentIndex = 0; currentIndex < ZI; currentIndex++)
			{
				//iidProb1.push_back(speciesKmerInfoA->HashPw[currentIndex]);
				//iidProb2.push_back(speciesKmerInfoB->HashPw[currentIndex]);
				iidProbAve.push_back(SciMultiple(SciAddition( (*(speciesKmerInfoA->probIIDPointer))[currentIndex], (*(speciesKmerInfoB->probIIDPointer))[currentIndex] ), 0.5));
			}
			
			entropyRate1 = JScomputeIID( (*(speciesKmerInfoA->probIIDPointer)) );
			entropyRate2 = JScomputeIID( (*(speciesKmerInfoB->probIIDPointer)) );
			entropyRateAve = JScomputeIID( iidProbAve );
			JSdivergence = entropyRateAve - 0.5 * entropyRate1 - 0.5 * entropyRate2;
			JSdistance = pow(JSdivergence, 0.5);
			
			
		}else{
			
			// order1 MUST equal order2
			//entropyRate1 = JScomputeMC(speciesKmerInfoA->iniProbPointer, speciesKmerInfoA->transMatrixPointer);
			//entropyRate2 = JScomputeMC(speciesKmerInfoB->iniProbPointer, speciesKmerInfoB->transMatrixPointer);
			entropyRate1 = JScomputeMC( (*(speciesKmerInfoA->iniProbPointer)), (*(speciesKmerInfoA->transMatrixPointer)) );
			entropyRate2 = JScomputeMC( (*(speciesKmerInfoB->iniProbPointer)), (*(speciesKmerInfoB->transMatrixPointer)) );
			//cout << TransToReal((*(speciesKmerInfoA->iniProbPointer))[0]) << endl;
			
			
			
			vector<SCIENTIFIC_NUMBER> iniProbAve;
			for(int currentIndex = 0; currentIndex < pow(ZI, orderA); currentIndex++)
			{
				iniProbAve.push_back(SciMultiple(SciAddition((*(speciesKmerInfoA->iniProbPointer))[currentIndex], (*(speciesKmerInfoB->iniProbPointer))[currentIndex]), 0.5));
				//cout << TransToReal(iniProbAve[currentIndex]) << endl;
			}
			
			vector< vector<SCIENTIFIC_NUMBER> > transMatrixAve(pow(ZI,orderA), vector<SCIENTIFIC_NUMBER>(ZI));
			for(int currentRow = 0; currentRow < pow(ZI, orderA); currentRow++)
			{
				for(int currentCol = 0; currentCol < ZI; currentCol++)
				{
					// compute JS distance
					transMatrixAve[currentRow][currentCol] = SciMultiple(SciAddition((*(speciesKmerInfoA->transMatrixPointer))[currentRow][currentCol], (*(speciesKmerInfoB->transMatrixPointer))[currentRow][currentCol]), 0.5);
					
				}
			}
			entropyRateAve = JScomputeMC(iniProbAve, transMatrixAve);
			
			JSdivergence = entropyRateAve - 0.5 * entropyRate1 - 0.5 * entropyRate2;
			JSdistance = pow(JSdivergence, 0.5);
			
		}
		
		//cout << entropyRate1 << ", " << entropyRate2 << ", ";
		//cout << entropyRateAve << ", " ;
		//cout << JSdivergence << ", " ;
		//cout << "JS, " << JSdistance << endl ;
		//fout << "JS, " << JSdistance << endl ;
		measureNames.push_back("JS");
		measureValues.push_back(JSdistance);
		
	}
	

	// 1. compute the D2 statistics
	string C2StatNames[3] = {"d2", "d2star", "d2shepp"} ;
	//cout << "== compute the D2 statistics, consider complementary == " << endl;
	vector<double> C2NGSvalues = D2C2computeNGS(ZI, k, speciesKmerInfoA, speciesKmerInfoB);
	//printFour(D2C2values);

	for(int i = 0; i < C2NGSvalues.size(); i++)
	{
		//cout << C2StatNames[i] << ", " << C2NGSvalues[i] << endl ;
		measureNames.push_back(C2StatNames[i]);
		measureValues.push_back(C2NGSvalues[i]);
		//cout << measureValues[1] << endl;
	}

	// 2. compute the S2 statistic
	//cout << "== compute the S2 statistics, considering single strand == " << endl;
	//double S2 = S2compute(ZI, k, speciesKmerInfoA, speciesKmerInfoB);
	//fout << "S2, " << S2 << endl;
	//cout << "S2, " << S2 << endl;

	// 3. compute HAO distance
	//cout << "== compute the HAO and Teeling statistic, considering single strand == " << endl;
	if(k < 3){
		//cerr << "ERROR: There is no Hao distance for k < 3!" << endl;
		
	}else{
		//cout << "==== compute the Hao statistics ==== " << endl;
		vector<double> HaoTeelingvalues = HAOTeelingcompute(ZI, k, speciesKmerInfoA, speciesKmerInfoB);
		//cout << "Hao, " <<  HaoTeelingvalues[0] << endl;
		//cout << "Teeling, " <<  HaoTeelingvalues[1] << endl;
		//fout << "hao, " << HAOvalue << endl;
		measureNames.push_back("Hao");
		measureValues.push_back(HaoTeelingvalues[0]);
		measureNames.push_back("Teeling");
		measureValues.push_back(HaoTeelingvalues[1]);
	}

	// 7. compute Willner
	//cout << "== compute the Willner distance == " << endl;
	if( k > 4)
	{
		//cout << "ERROR: no definition of Willner for k>4." << endl;
		
	}else{
		double willner = WillnerDiNGS(ZI, k, speciesKmerInfoA, speciesKmerInfoB) ;
		//cout << "Willner, " << willner << endl ;
		//fout << "willner, " << willner << endl ;
		measureNames.push_back("Willner");
		measureValues.push_back(willner);
	}
	

}









///////////////////////////////////////////////////////////////////////////////////////////
//KmerCount.out [k] [sample data file]
int main(int argc, char **argv)   //EDIT main(int argc, char *argv[])
{
	//bool doubleStrand = false;
	//char *speciesName[2];

	int k = 0;
	char kstr[5];
	
	string speciesInfoFilePathNameA;
	string speciesInfoFilePathNameB;
	string taxaFile;
	string outDIR;
	string htmlTmpDIR;
	
	vector<string> *measureNames = new vector<string>;
	
	//vector<int> order(2,0);
	
	//char *inputFileDir[2];
	//inputFileDir[0] = NULL;
	//inputFileDir[1] = NULL;
	
	
	/* getOptions from command line */
	int c;
	
	while (1)
	{
		static struct option long_options[] =
		{
			//{"doubleStrand",     no_argument,       0, 'd'},
			// necessary arguments: kvalue and inputFileName !!
			//{"species1",  required_argument, 0, 'a'},
			//{"order1",  required_argument, 0, 'b'},
			//{"species2",  required_argument, 0, 'c'},
			//{"order2",  required_argument, 0, 'd'},
			{"kvalue",  required_argument, 0, 'k'},
			{"inputFileDir1",  required_argument, 0, 'i'},
			{"inputFileDir2",  required_argument, 0, 'j'},
			{"taxaFile",  required_argument, 0, 't'},
			{"outDir",  required_argument, 0, 'o'},
			{0,         0,                 0,  0 }
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		c = getopt_long (argc, argv, "k:i:j:t:o:",
										 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1)
			break;
		
		switch (c)
		{

			//case 'd':
			//	doubleStrand = true;
			//	puts ("option -d: consider complimentary word. ");
			//	break;
				
			/* input kmer files should be singleStrand count! */
			 
			//case 'a':
				//speciesName[0] = optarg;
				//printf ("option -a, name of 1st species, with value `%s'\n", optarg);
				//break;
				
			//case 'b':
			//	order[0] = atoi(optarg);
				//printf ("option -b, order of 1st species, with value `%s'\n", optarg);
				//break;
				
			//case 'c':
				//speciesName[1] = optarg;
				//printf ("option -c, name of 2nd species, with value `%s'\n", optarg);
				//break;
				
			//case 'd':
				//order[1] = atoi(optarg);
				//printf ("option -d, order of 2nd species, with value `%s'\n", optarg);
				//break;
				
			case 'k':
				k = atoi(optarg) ;
				sprintf(kstr, "%d", k);
				//printf ("option -k, the value of k, with value `%s'\n", optarg);
				break;
				
			case 'i':
				speciesInfoFilePathNameA = string(optarg);
				//printf ("option -i, the input kmer count file directory, with value `%s'\n", optarg);
				break;
				
			case 'j':
				speciesInfoFilePathNameB = string(optarg);
				//printf ("option -j, the input kmer count file directory, with value `%s'\n", optarg);
				break;
				
			case 't':
				taxaFile = string(optarg);
				//printf ("option -t, the taxonomy file for host contigs, with value `%s'\n", optarg);
				break;
				
			case 'o':
				outDIR = string(optarg);
				//printf ("option -o, the output directory, with value `%s'\n", optarg);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				abort ();
		}
	}

	
	/////////// preparation for visualization: loading taxa files ////////////////
	// cp Yang's visualization files to outDIR
	string cmdString = string(argv[0]);
	size_t found = cmdString.find_last_of("/");
	string cmdDIR = cmdString.substr(0,found);

	string cpCssCMD = "cp -r " + cmdDIR + "/css" + " " + outDIR;
	system(cpCssCMD.c_str());
	string cpLibCMD = "cp -r " + cmdDIR + "/lib" + " " + outDIR;
	system(cpLibCMD.c_str());
	string cpLogoCMD = "cp -r " + cmdDIR + "/logo.jpg" + " " + outDIR;
	system(cpLogoCMD.c_str());
	
	htmlTmpDIR = outDIR + "/" + "tmp_html";
  system(("mkdir -p " + htmlTmpDIR).c_str());


	///////////////////////////////////////////////////////////////////////
	//////////////////// load the species information /////////////////////
	///////////////////////////////////////////////////////////////////////
		//cerr << "......computing measures......" << endl;
	
	vector<SPECIESINFO> speciesInfoListA;
	loadSpeciesInfo(speciesInfoFilePathNameA, speciesInfoListA);
	
	vector<SPECIESINFO> speciesInfoListB;
	loadSpeciesInfo(speciesInfoFilePathNameB, speciesInfoListB);
	
	///////////////////////////////////////////////////////////////////////
	/////////////////////// output matrix /////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	vector<vector<vector<double> > > resultMatrix (11, vector<vector<double> > (speciesInfoListA.size(), vector<double> (speciesInfoListB.size())));
	
//	vector<vector<vector<double> > > resultMatrix;
//	// Set up sizes. (HEIGHT x WIDTH)
//	resultMatrix.resize(11);
//	for (int i = 0; i < 11; ++i) {
//		resultMatrix[i].resize(speciesInfoListA.size());
//		
//		for (int j = 0; j < speciesInfoListA.size(); ++j)
//			resultMatrix[i][j].resize(speciesInfoListB.size());
//	}

	/////////////////////////////////////////////////////////////////////////////
	////////////////////////// load taxa file ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	vector<string> hostNCBIName;
	vector<string> hostName;
	vector<string> hostSuperkingdom;
	vector<string> hostPhylum;
	vector<string> hostClass;
	vector<string> hostOrder;
	vector<string> hostFamily;
	vector<string> hostGenus;
	vector<string> hostSpecies;
	int taxaLineNum = loadTaxaInfo(taxaFile, hostNCBIName, hostName, hostSuperkingdom, hostPhylum, hostClass, hostOrder, hostFamily, hostGenus, hostSpecies);
	
	if( taxaLineNum != speciesInfoListB.size() )
	{
		cerr << "ERROR: number of hosts in taxa file is not equal to number of host fasta files " << endl;
		return 0;
	}
	////////// sort the taxa info by NCBINames in the input hostList file /////////
	vector<string> hostName_sort;
	vector<string> hostSuperkingdom_sort;
	vector<string> hostPhylum_sort;
	vector<string> hostClass_sort;
	vector<string> hostOrder_sort;
	vector<string> hostFamily_sort;
	vector<string> hostGenus_sort;
	vector<string> hostSpecies_sort;
	for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
	{
		SPECIESINFO speciesInfoB = speciesInfoListB[IDB];
		string currentHostNCBIName = speciesInfoB.name;
		vector<string>::iterator it=find(hostNCBIName.begin(),hostNCBIName.end(),currentHostNCBIName);
		int pos = distance(hostNCBIName.begin(), it);
		if( pos == hostNCBIName.size() )
		{
			// can't find one!
			cerr << "ERROR: can't find taxa info for host " << currentHostNCBIName << "! \nThe host names in host taxa file should be the full file name of the host fasta file (including filename extensions). \nPlease double check. " << endl;
			return 0;
		}
		//cout << currentHostNCBIName << ":" << pos << endl;
		hostName_sort.push_back(hostName[pos]);
		hostSuperkingdom_sort.push_back(hostSuperkingdom[pos]);
		hostPhylum_sort.push_back(hostPhylum[pos]);
		hostClass_sort.push_back(hostClass[pos]);
		hostOrder_sort.push_back(hostOrder[pos]);
		hostFamily_sort.push_back(hostFamily[pos]);
		hostGenus_sort.push_back(hostGenus[pos]);
		hostSpecies_sort.push_back(hostSpecies[pos]);
	}
	
	
	///////////////////////////////////////////////////////////////////////
	//////////////////// load the kmerFiles and compute pwMC //////////////
	///////////////////////////////////////////////////////////////////////

  //vector<SCIENTIFIC_NUMBER> iniProbA;
  //vector<vector<vector<SCIENTIFIC_NUMBER> > > transMatrix (2, vector<vector<SCIENTIFIC_NUMBER> > (pow(ZI,order[1]), vector<SCIENTIFIC_NUMBER> (ZI)));
  
  // load the kmer count hashtables of the first file list (the fat list)
	// those files are only loaded once. They will be repeatedly used for many times
	// how many files you put on this list depend on the memory size
	
	ofstream outMatrixBinFile;
	string matrixBinFileName = outDIR + "/tmp/resultMatrix.bin";
	outMatrixBinFile.open(matrixBinFileName.c_str(), ios::out|ios::binary);

	vector<KMERINFO*> speciesKmerInfoListA( speciesInfoListA.size() );
	for(int IDA = 0; IDA < speciesInfoListA.size(); IDA++)
	{
		// initialization
		speciesKmerInfoListA[IDA] = new KMERINFO();
	}
	
	//cout << "=== start loading kmers & computing Pw in all species in ListA === " << endl;
	for(int IDA = 0; IDA < speciesInfoListA.size(); IDA++)
	{
		SPECIESINFO speciesInfoA = speciesInfoListA[IDA];

		// load k, k-1, k-2 mer info
		char k_1str[5];
		sprintf(k_1str, "%d", k-1);
		char k_2str[5];
		sprintf(k_2str, "%d", k-2);
		//cout << speciesKmerInfoA->totalK_1 << endl;
		
		//cout << "======= species " << speciesInfoA.name << " order " << speciesInfoA.order << endl;
		//KMERINFO* speciesKmerInfoA = new KMERINFO();
		std::string kmerFilePathNameA = speciesInfoA.dir + speciesInfoA.name + "_k" + kstr + "_ss_wc";
		loadKmerCountHash(kmerFilePathNameA, speciesKmerInfoListA[IDA], "kmerCount");
		
		// !!! 20150629, totalK_1/totalK_2 are not the same as totalOrder/totalOrder_1 !!!
		//cout << "==== HAO: load the kmer count for (k-1) = " << k_1str << " ==" << endl;
		std::string kmerFilek_1PathNameA = speciesInfoA.dir + speciesInfoA.name + "_k" + k_1str + "_ss_wc";
		//cout << "kmerFile: " << kmerFilek_1PathNameA << endl;
		loadKmerCountHash(kmerFilek_1PathNameA, speciesKmerInfoListA[IDA], "kmerCount_1");
		
		// load the kmer count for MC: k=order+1
		//cout << "==== HAO: load the kmer count for (k-2) = " << k_2str << " ==" << endl;
		std::string kmerFilek_2PathNameA = speciesInfoA.dir + speciesInfoA.name + "_k" + k_2str + "_ss_wc";
		//cout << "kmerFile: " << kmerFilek_2PathNameA << endl;
		loadKmerCountHash(kmerFilek_2PathNameA, speciesKmerInfoListA[IDA], "kmerCount_2");

		// load order, order+1 mer info, and compute Pw
		int orderA = speciesInfoA.order;
		char orderstrA[5];
		sprintf(orderstrA, "%d", orderA);
		char order1strA[5];
		sprintf(order1strA, "%d", (orderA+1));

		if(orderA == 0)
		{
			// load the 1-kmer count for MC: k=order+1
			//cout << "== load the kmer count for MC, order is " << speciesInfoListA[IDA]->order << " ==" << endl;
			std::string kmerFileOrder1PathNameA = speciesInfoA.dir + speciesInfoA.name + "_k" + order1strA + "_ss_wc";
			//cout << "kmerFile: " << kmerFileOrder1PathName << endl;
			loadKmerCountHash(kmerFileOrder1PathNameA, speciesKmerInfoListA[IDA], "kmerOrder+1");
			// compute the pw
			//cout << "== compute the pwIID for the words, for D2-type and S2 == " << endl;
			//20150727 JS IID correct
			vector<SCIENTIFIC_NUMBER> *probIIDA = new vector<SCIENTIFIC_NUMBER>;
			pwIID(ZI, k, speciesKmerInfoListA[IDA], (*probIIDA));
			
			speciesKmerInfoListA[IDA]->probIIDPointer = probIIDA;

		}else{
			
			// load the kmer count for MC: k=order
			//cout << "== load the kmer count for MC, order is " << order << " ==" << endl;
			std::string kmerFileOrderPathNameA = speciesInfoA.dir + speciesInfoA.name + "_k" + orderstrA + "_ss_wc";
			//cout << "kmerFile: " << kmerFileOrderPathName << endl;
			loadKmerCountHash(kmerFileOrderPathNameA, speciesKmerInfoListA[IDA], "kmerOrder");
			
			// load the kmer count for MC: k=order+1
			//cout << "== load the kmer count for MC, order+1 is " << order1str << " ==" << endl;
			std::string kmerFileOrder1PathNameA = speciesInfoA.dir + speciesInfoA.name + "_k" + order1strA + "_ss_wc";
			//cout << "kmerFile: " << kmerFileOrder1PathName << endl;
			loadKmerCountHash(kmerFileOrder1PathNameA, speciesKmerInfoListA[IDA], "kmerOrder+1");
			
			// compute the pw
			
			//cout << "== compute the pw for the words == " << endl;
			//vector<SCIENTIFIC_NUMBER> iniProbA;
			//iniProbA = new vector<SCIENTIFIC_NUMBER>;
			//vector<int> AA = new vector<int>(2,0);
			vector<SCIENTIFIC_NUMBER> *iniProbA = new vector<SCIENTIFIC_NUMBER>;
			vector< vector<SCIENTIFIC_NUMBER> > *transMatrixA = new vector< vector<SCIENTIFIC_NUMBER> >(pow(ZI, orderA), vector<SCIENTIFIC_NUMBER>(ZI));
			
			//vector< vector<SCIENTIFIC_NUMBER> > transMatrixA(pow(ZI, orderA), vector<SCIENTIFIC_NUMBER>(ZI));
			//cout << "here" << endl;
			pwMC(ZI, k, orderA, speciesKmerInfoListA[IDA], (*iniProbA), (*transMatrixA));
			
			speciesKmerInfoListA[IDA]->iniProbPointer = iniProbA;
			speciesKmerInfoListA[IDA]->transMatrixPointer = transMatrixA;
			//cout << "A:" << TransToReal((*(speciesKmerInfoListA[IDA]->iniProbPointer))[0]) << endl;
			
		}
	}
	
	
	//cout << "A1:" << TransToReal((*(speciesKmerInfoListA[0]->iniProbPointer))[0]) << " "<< TransToReal((*(speciesKmerInfoListA[0]->iniProbPointer))[1]) << " " << TransToReal((*(speciesKmerInfoListA[0]->iniProbPointer))[2]) << " "<< TransToReal((*(speciesKmerInfoListA[0]->iniProbPointer))[3]) << endl;
	//cout << "A2:" << TransToReal((*(speciesKmerInfoListA[1]->iniProbPointer))[0]) << " " << TransToReal((*(speciesKmerInfoListA[1]->iniProbPointer))[1]) << " " << TransToReal((*(speciesKmerInfoListA[1]->iniProbPointer))[2]) << " " << TransToReal((*(speciesKmerInfoListA[1]->iniProbPointer))[3]) << endl;
	
	// then for each file on the second list, (the light list)
	// compute the score between it and all the files in the first list
	for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++)
	{
		//cout << "=== start loading kmers & computePw in one of ListB === " << endl;
		
		// this varaible "speciesKmerInfoB" will be rewritten by the next file on the list B
		KMERINFO* speciesKmerInfoB = new KMERINFO();
		SPECIESINFO speciesInfoB = speciesInfoListB[IDB];
		
		// load k, k-1, k-2 mer info
		char k_1str[5];
		sprintf(k_1str, "%d", k-1);
		char k_2str[5];
		sprintf(k_2str, "%d", k-2);
		
		//cout << "======= species " << speciesInfoB.name << " order " << speciesInfoB.order << endl;
		std::string kmerFilePathNameB = speciesInfoB.dir + speciesInfoB.name + "_k" + kstr + "_ss_wc";
		loadKmerCountHash(kmerFilePathNameB, speciesKmerInfoB, "kmerCount");

		// !!! 20150629, totalK_1/totalK_2 are not the same as totalOrder/totalOrder_1 !!!
		//cout << "==== HAO: load the kmer count for (k-1) = " << k_1str << " ==" << endl;
		std::string kmerFilek_1PathNameB = speciesInfoB.dir + speciesInfoB.name + "_k" + k_1str + "_ss_wc";
		//cout << "kmerFile: " << kmerFilek_1PathNameB << endl;
		loadKmerCountHash(kmerFilek_1PathNameB, speciesKmerInfoB, "kmerCount_1");
		//cout << "HashTableK_1 " << " word " << currentKmerRemoveLastTen << " kmercount " << speciesKmerInfoB->HashTableK_1[currentKmerRemoveLastTen] << " totalK_1 " << speciesKmerInfoB->totalK_1 << endl;
		
		// load the kmer count for MC: k=order+1
		//cout << "==== HAO: load the kmer count for (k-2) = " << k_2str << " ==" << endl;
		std::string kmerFilek_2PathNameB = speciesInfoB.dir + speciesInfoB.name + "_k" + k_2str + "_ss_wc";
		//cout << "kmerFile: " << kmerFilek_2PathNameB << endl;
		loadKmerCountHash(kmerFilek_2PathNameB, speciesKmerInfoB, "kmerCount_2");
		
		
		// load order, order+1 mer info, and compute Pw
		int orderB = speciesInfoB.order;
		char orderstrB[5];
		sprintf(orderstrB, "%d", orderB);
		char order1strB[5];
		sprintf(order1strB, "%d", (orderB+1));
		
		vector<SCIENTIFIC_NUMBER> *probIIDB = new vector<SCIENTIFIC_NUMBER>;
		vector<SCIENTIFIC_NUMBER> *iniProbB = new vector<SCIENTIFIC_NUMBER>;
		vector< vector<SCIENTIFIC_NUMBER> > *transMatrixB = new vector< vector<SCIENTIFIC_NUMBER> >(pow(ZI, orderB), vector<SCIENTIFIC_NUMBER>(ZI));
		if(orderB == 0)
		{
			// load the 1-kmer count for MC: k=order+1
			//cout << "== load the kmer count for MC, order is " << speciesInfoListA[IDA]->order << " ==" << endl;
			std::string kmerFileOrder1PathNameB = speciesInfoB.dir + speciesInfoB.name + "_k" + order1strB + "_ss_wc";
			//cout << "kmerFile: " << kmerFileOrder1PathName << endl;
			loadKmerCountHash(kmerFileOrder1PathNameB, speciesKmerInfoB, "kmerOrder+1");
			// compute the pw
			//cout << "== compute the pwIID for the words, for D2-type and S2 == " << endl;
			
			//20150727 JS IID correct
			//vector<SCIENTIFIC_NUMBER> *probIIDB = new vector<SCIENTIFIC_NUMBER>;
			pwIID(ZI, k, speciesKmerInfoB, (*probIIDB));
			speciesKmerInfoB->probIIDPointer = probIIDB;

			
		}else{
			
			// load the kmer count for MC: k=order
			//cout << "== load the kmer count for MC, order is " << order << " ==" << endl;
			std::string kmerFileOrderPathNameB = speciesInfoB.dir + speciesInfoB.name + "_k" + orderstrB + "_ss_wc";
			//cout << "kmerFile: " << kmerFileOrderPathName << endl;
			loadKmerCountHash(kmerFileOrderPathNameB, speciesKmerInfoB, "kmerOrder");
			
			// load the kmer count for MC: k=order+1
			//cout << "== load the kmer count for MC, order+1 is " << order1str << " ==" << endl;
			std::string kmerFileOrder1PathNameB = speciesInfoB.dir + speciesInfoB.name + "_k" + order1strB + "_ss_wc";
			//cout << "kmerFile: " << kmerFileOrder1PathName << endl;
			loadKmerCountHash(kmerFileOrder1PathNameB, speciesKmerInfoB, "kmerOrder+1");
			
			//cout << "order1strB" << order1strB << endl;
			//unsigned long currentKmerRemoveLastTen = 0;
			//cout << "HashTableK_1 " << " word " << currentKmerRemoveLastTen << " kmercount " << speciesKmerInfoB->HashTableK_1[currentKmerRemoveLastTen] << " totalK_1 " << speciesKmerInfoB->totalK_1 << endl;
			
			// compute the pw
			//cout << "== compute the pw for the words == " << endl;
			//vector<SCIENTIFIC_NUMBER> iniProbB;
			//vector< vector<SCIENTIFIC_NUMBER> > transMatrixB(pow(ZI, orderB), vector<SCIENTIFIC_NUMBER>(ZI));
			//vector<SCIENTIFIC_NUMBER> *iniProbB = new vector<SCIENTIFIC_NUMBER>;
			//vector< vector<SCIENTIFIC_NUMBER> > *transMatrixB = new vector< vector<SCIENTIFIC_NUMBER> >(pow(ZI, orderB), vector<SCIENTIFIC_NUMBER>(ZI));
			pwMC(ZI, k, orderB, speciesKmerInfoB, (*iniProbB), (*transMatrixB));
			speciesKmerInfoB->iniProbPointer = iniProbB;
			speciesKmerInfoB->transMatrixPointer = transMatrixB;
			//cout << "B:" << TransToReal((*(speciesKmerInfoB->iniProbPointer))[0]) << endl;
			
		}
		
		
		//unsigned long currentKmerRemoveLastTen = 0;
		//cout << "HashTableK_1 " << " word " << currentKmerRemoveLastTen << " kmercount " << speciesKmerInfoB->HashTableK_1[currentKmerRemoveLastTen] << " totalK_1 " << speciesKmerInfoB->totalK_1 << endl;

		///////////////////////////////////////////////////////////////////////
		////////////////////// computing the multiple statistics //////////////
		///////////////////////////////////////////////////////////////////////
		cerr << " dissimilarity b/w host " << speciesInfoB.name << " - all viruses " << endl;
		
		for(int IDA = 0; IDA < speciesInfoListA.size(); IDA++)
		{
			KMERINFO* speciesKmerInfoA = speciesKmerInfoListA[IDA];
			SPECIESINFO speciesInfoA = speciesInfoListA[IDA];
			
			//cout << "=== start computing statistics === " << endl;
			//cout << "speciesA:" << IDA << endl;
			//cout << "speciesB:" << IDB << endl;
			//cerr << " pair " << speciesInfoA.name << " - " << speciesInfoB.name << endl;
			
			vector<double> *measureValues = new vector<double>;
			measureNames->clear();
			//cout << "check1" << endl;
			computeMultiStatsReturn(ZI, k, speciesKmerInfoA, speciesKmerInfoB, speciesInfoA, speciesInfoB, (*measureNames), (*measureValues));
			//cout << "check2" << endl;
			
			//vector<string> &measureNamesRef = *measureNames;
			//vector<double> &measureValuesRef = *measureValues;
			//cout << measureValuesRef[0] << endl;
			for( int statID = 0; statID < measureValues->size(); statID++ )
			{
				resultMatrix[statID][IDA][IDB] = measureValues->at(statID);
				outMatrixBinFile.write( (char*)&resultMatrix[statID][IDA][IDB], sizeof(double) );
			}
			//delete speciesKmerInfoA;
		
		}
		
		delete speciesKmerInfoB;
		delete probIIDB;
		delete iniProbB;
		delete transMatrixB;
		
//		ofstream outBinFile;
//		string binFileName = outDIR + "/resultMatrix.bin";
//		outBinFile.open(binFileName, ios::out | ios::binary);
//		outBinFile.write (reinterpret_cast<char*> (&resultMatrix), sizeof(int)) ;
//		outBinFile.close();
		
	}
	
	outMatrixBinFile.close();

//	/////////////////////////////////////////////////////////////////////////////
//	////////////////////////// load taxa file ////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////
//	vector<string> hostNCBIName;
//	vector<string> hostName;
//	vector<string> hostSuperkingdom;
//	vector<string> hostPhylum;
//	vector<string> hostClass;
//	vector<string> hostOrder;
//	vector<string> hostFamily;
//	vector<string> hostGenus;
//	vector<string> hostSpecies;
//	int taxaLineNum = loadTaxaInfo(taxaFile, hostNCBIName, hostName, hostSuperkingdom, hostPhylum, hostClass, hostOrder, hostFamily, hostGenus, hostSpecies);
//	
//	if( taxaLineNum != speciesInfoListB.size() )
//	{
//		cerr << "number of hosts in taxa file is not equal to number of host fasta files " << endl;
//		return 0;
//	}
//	////////// sort the taxa info by NCBINames in the input hostList file /////////
//	vector<string> hostName_sort;
//	vector<string> hostSuperkingdom_sort;
//	vector<string> hostPhylum_sort;
//	vector<string> hostClass_sort;
//	vector<string> hostOrder_sort;
//	vector<string> hostFamily_sort;
//	vector<string> hostGenus_sort;
//	vector<string> hostSpecies_sort;
//	for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
//	{
//		SPECIESINFO speciesInfoB = speciesInfoListB[IDB];
//		string currentHostNCBIName = speciesInfoB.name;
//		vector<string>::iterator it=find(hostNCBIName.begin(),hostNCBIName.end(),currentHostNCBIName);
//		int pos = distance(hostNCBIName.begin(), it);
//		//cout << currentHostNCBIName << ":" << pos << endl;
//		hostName_sort.push_back(hostName[pos]);
//		hostSuperkingdom_sort.push_back(hostSuperkingdom[pos]);
//		hostPhylum_sort.push_back(hostPhylum[pos]);
//		hostClass_sort.push_back(hostClass[pos]);
//		hostOrder_sort.push_back(hostOrder[pos]);
//		hostFamily_sort.push_back(hostFamily[pos]);
//		hostGenus_sort.push_back(hostGenus[pos]);
//		hostSpecies_sort.push_back(hostSpecies[pos]);
//	}
	
	
	/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// output //////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	for( int statID = 0; statID < measureNames->size(); statID++ )
	{
		string statName = measureNames->at(statID);
		//cerr << "...measure:" << statName << "..." << endl;
		
		ofstream csvOut;
		string csvFileName = outDIR + "/" + statName + "_k" + kstr + ".csv";
		csvOut.open(csvFileName.c_str());
		
		ofstream html2Out;
		string html2FileName = htmlTmpDIR + "/" + statName + "_k" + kstr + "_main.html_part2";
		html2Out.open(html2FileName.c_str());
		
		// create variables for Yang's visualization: virusNameArr
		html2Out << "var virusNameArr = [";
		for(int IDA = 0; IDA < speciesInfoListA.size(); IDA++)
		{
			SPECIESINFO speciesInfoA = speciesInfoListA[IDA];
			html2Out << "\"" << speciesInfoA.name << "\"";
			if(IDA != speciesInfoListA.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		

		// first row: statName and colnames
		csvOut << statName << ",";
		
		// create variables for Yang's visualization: hostNameArr
		html2Out << "var hostNameArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++)
		{
			SPECIESINFO speciesInfoB = speciesInfoListB[IDB];
			csvOut << speciesInfoB.name << ",";
			html2Out << "\"" << speciesInfoB.name << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		csvOut << endl;
		html2Out << "];" << endl;
		
		/////////////// output the matrix: csv files //////////////
		// create variables for Yang's visualization: virusHostMat
		html2Out << "var virusHostMat = [";
		for(int IDA = 0; IDA < speciesInfoListA.size(); IDA++)
		{
			SPECIESINFO speciesInfoA = speciesInfoListA[IDA];
			csvOut << speciesInfoA.name << ",";
			html2Out << "[";
			for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++)
			{
				SPECIESINFO speciesInfoB = speciesInfoListB[IDB];
				csvOut << resultMatrix[statID][IDA][IDB] << ",";
				html2Out << resultMatrix[statID][IDA][IDB];
				if( IDB != speciesInfoListB.size()-1 )
				{
					html2Out << ",";
				}
			}
			html2Out << "]";
			csvOut << endl;
			if(IDA != speciesInfoListA.size()-1 )
			{
				html2Out << "," << endl;
			}
		}
		csvOut.close();
		html2Out << "];" << endl;
		
		
		//////////////////// preparation for visualization ////////////////////////
		// create variable for Yang's visualization: host taxonomy

		
		
		html2Out << "var superkingdomArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostSuperkingdom_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var phylumArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostPhylum_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var classArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostClass_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var orderArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostOrder_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var familyArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostFamily_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var genusArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostGenus_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var speciesArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostSpecies_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		html2Out << "var organismArr = [";
		for(int IDB = 0; IDB < speciesInfoListB.size(); IDB++  )
		{
			html2Out << "\"" << hostName_sort[IDB] << "\"";
			if(IDB != speciesInfoListB.size()-1 ){
				html2Out << ", ";
			}
		}
		html2Out << "];" << endl;
		
		
		// cat html_part1, 2 and 3
		string html1FileName = cmdDIR + "/main.html_part1";
		string html3FileName = cmdDIR + "/main.html_part3";
		string htmlFileName = outDIR + "/" + statName + "_k" + kstr + "_main.html";
		string catCMD = "cat " + html1FileName + " " + html2FileName + " " + html3FileName + " > " + htmlFileName;
		system(catCMD.c_str());
		
	}
	
	system(("rm -rf " + htmlTmpDIR).c_str());
	
	return 0;
}






