/*
 *  Tools.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "Tools.h"
#include <math.h>
//#include <sys/time.h>
#include "Variables.h"
#include <stdint.h>
#include <float.h>
#include <climits>

/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size)
{
	//int flags[size];
	int * flags=new int[size];
	for (int i=0;i<size;i++) flags[i]=1;
	
	for (int i=0;i<size;i++)
	{
		int value=permutation[i];
		flags[value]=0;
	}
	
	int result,sum=0;
	for(int i=0;i<size;i++)
		sum+=flags[i];
	if (sum==0) result=true;
	else result=false;
	//delete [] flags;
	return result;
}

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int * permutation, int n)
{
	for (int i = 0; i < n; ++i)
	{
		int j = rand() % (i + 1);
		permutation[i] = permutation[j];
		permutation[j] = i;
	}
}

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr)
{
	size_t contains;
	
	contains = inputStr.find(searchStr);
	
	if(contains != string::npos)
		return true;
	else
		return false;
}

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double* array, int length, string text)
{
	cout<<text;
	for (int i=0;i<length;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(int * array, int size)
{
	for (int i=0;i<size;i++){
		cout<<array[i]<<" ";
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(float* array, int length, string text)
{
	cout<<text<<": ";
	for (int i=0;i<length;i++){
        printf(" %3.4f,",array[i]);
	}
	printf("\n ");

}

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(double* array, int length, string text)
{
    int i;
	cout<<text;
	for (i=0;i<length;i++)
		printf(" %3.4f ",array[i]);
	printf("\n ");
}

/*
 * Prints in the standard output given matrix.
 */
void PrintMatrix(int** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text;
	for (i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}

void PrintMatrix(float** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text;
	for (i=0;i<length;i++)
	{
		cout<<""<<endl;
		for (j=0;j<length2;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
	}
	cout<<" "<<endl;
}


/*
 * Prints in standard output 'lengthxlength' double elements of a given matrix.
 */
void PrintMatrixDouble(double** matrix, int length, int length2, string text)
{
    int i,j;
	cout<<text<<endl;
	for (i=0;i<length;i++)
	{
		for (j=0;j<length2;j++)
			printf("%3.9f, ",matrix[i][j]);
		printf("\n");
	}
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size)
{
    int i,dist;
    int * v,*composition,*invertedB,*m_aux;
    
    dist=0;
    v=new int[size-1];
    composition=new int[size];
    invertedB=new int[size];
    m_aux=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);
    
    for (i = 0; i < size-1; i++)
        dist += v[i];
    
    delete [] composition;
    delete [] invertedB;
    delete [] v;
    delete [] m_aux;
    
    return dist;
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux)
{
    int i,dist;
    int * v,*composition,*invertedB;
    
    dist=0;
    v=new int[size-1];
    composition=new int[size];
    invertedB=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);

    for (i = 0; i < size-1; i++)
        dist += v[i];
    
    delete [] composition;
    delete [] invertedB;
    delete [] v;
    
    return dist;
}

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple continuous executions.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux, int * invertedB, int *  composition, int * v)
{
    int i,dist;
    dist=0;
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    vVector_Fast(v,composition,size,m_aux);
    for (i = 0; i < size-1; i++)
        dist += v[i];

    return dist;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size)
{
    int * invertedB= new int[size];
    int * composition= new int[size];
    int * elemsToCycles= new int[size];
    int * maxPosInCycle= new int[size];
    int * freeCycle= new int[size];

    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    
    int index,cycle,distance;
    
    for(int i=0;i<size;i++)
    {
        elemsToCycles[i]=-1;
        maxPosInCycle[i]=-1;
        freeCycle[i]=1;
    }
    
    while((index=NextUnasignedElem(elemsToCycles,size))!=-1)
    {
        cycle=FindNewCycle(freeCycle,size);
        freeCycle[cycle]=0;
        do
        {
            elemsToCycles[index]=cycle;
            index = composition[index];//para permus de 1..n =>index = sigma[index]-1;
        }
        while(elemsToCycles[index] == -1);
    }
    distance=size-FindNewCycle(freeCycle,size);
    
    delete [] invertedB;
    delete [] composition;
    delete [] elemsToCycles;
    delete [] maxPosInCycle;
    delete [] freeCycle;
    
    return distance;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size,int * invertedB, int *  composition, int * elemsToCycles, int * maxPosInCycle, int * freeCycle)
{
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);
    
    int index,cycle,distance;
    
    for(int i=0;i<size;i++)
    {
        elemsToCycles[i]=-1;
        maxPosInCycle[i]=-1;
        freeCycle[i]=1;
    }
    
    while((index=NextUnasignedElem(elemsToCycles,size))!=-1)
    {
        cycle=FindNewCycle(freeCycle,size);
        freeCycle[cycle]=0;
        do
        {
            elemsToCycles[index]=cycle;
            index = composition[index];//para permus de 1..n =>index = sigma[index]-1;
        }
        while(elemsToCycles[index] == -1);
    }
    distance=size-FindNewCycle(freeCycle,size);
    
    return distance;
}


int FindNewCycle(int * freeCycle, int size){
    int i;
    
    
    for(i=0;i<size;i++)
        if(freeCycle[i])
            return i;
    return size;
}

int NextUnasignedElem(int * elemsToCycles, int size)
{
    int i;
    for(i=0;i<size;i++)
        if(elemsToCycles[i]==-1)
            return i;
    return -1;
}

/*
 * Calculates the Ulam distance between 2 permutations.
 */
int Ulam(int * permutationA, int * permutationB, int size){

    int * composition=new int[size];
    int * invertedB=new int[size];
    
    Invert(permutationB, size, invertedB);
    Compose(permutationA,invertedB,composition,size);

    int dist=(size- getLISLength(composition,size));
    
    delete [] composition;
    delete [] invertedB;
    
    return dist;
}

/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int*sigma,int size){
   
    // O(n log k)
    
    int i;
    vector<int> vc(1,sigma[0]);
    vector<int>::iterator vk;
    
    for (i=1;i<size;i++){
        for (vk=vc.begin(); vk != vc.end(); vk++)
            if (*vk>=sigma[i]) break;
        if (vk==vc.end())
            vc.push_back(sigma[i]);
        else *vk=sigma[i];
    }
    
    return (int)vc.size();
    
}
/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n)
{
    int i;
    for(i=0;i<n;i++)
        res[i]=s1[s2[i]];
}


/*
 * Calculates V_j-s vector.
 */
void vVector(int*v, int*permutation, int n)
{
    
    int i,j;
    for(i=0;i<n-1;i++)
        v[i]=0;
    
    for (i = n-2; i >= 0; i--)
        for (j = i+1; j < n; j++)
            if(permutation[i] > permutation[j])
                v[i]++;
}

/*
 *  Optimized version proposed by Leti for the calculation of the V_j-s vector.
 */
void vVector_Fast(int*v, int*permutation, int n, int * m_aux)
{
    int i,j, index;
    
    for(i=0;i<n-1;i++){
        v[i]=0;
        m_aux[i]=0;
    }
    m_aux[n-1]=0;
    for (j=0; j<n-1; j++){
        index=permutation[j];
        v[j]=index-m_aux[index];
        for (i=index; i<n; i++)
            m_aux[i]++;
    }
}

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted)
{
    int i;
    for(i=0; i<n; i++)
        inverted[permu[i]]=i;
}

/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys( int * a, float * criteriaValues, int size)
{
	bool * fixedValues= new bool[size];
	float criteria, min;
	int i, j;
	for (i=0;i<size;i++){
		fixedValues[i]=false;
		a[i]=0;
	}
	int minPos=0;
	for (i=0;i<size;i++)
	{
		min=INT_MAX;
		for (j=0;j<size;j++)
		{
			criteria=criteriaValues[j];
			if (!fixedValues[j] && min>criteria )
			{
				min=criteria;
				minPos=j;
			}
		}
        
		fixedValues[minPos]=true;
		//a[i]=minPos;// modification por el asunto ordering /ranking
    	a[minPos]=i;// original.
	}
	delete [] fixedValues;
}



/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val) {
    if(val <= 0) return 1;
    //long  N, b, c, p; // use int for fast calculation and small range of calculation..
    long   b, c;
    long double p, N;
    N=(long double)val;
    c = (long)N - 1;
    p = 1;
    while (c > 0) {
        p = 0;
        b = c;
        while (b > 0) {
            if (b & 1) {
                p += N; // p = p + N;
            }
            // if you would like to use double choose the alternative forms instead shifts
            // the code is fast even!
            // you can use the same tips on double or 64 bit int etc.... but you must... ;-)
            //b >>= 1; // b/=2; (b = b / 2;) ( b >> 1; a.s.r. is more efficent for int or long..!)
            b/=2;
            //N <<= 1; // N += N; N = N + N; N = N * 2; (N <<=1; a.s.l. is more efficent for int or long..!)
            N += N;
        } // end of: while(b>0)
        N = p;
        c--; // c = c - 1;
    } // end of: while(c > 0)
    //printf("[%d] is the factorial! \n", p);
    return p;
}

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void Swap(int * array, int i, int j)
{
	int aux=array[i];
	array[i]=array[j];
	array[j]=aux;
}
