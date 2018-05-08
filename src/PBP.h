/*
 *  PBP.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */
#ifndef __PBP_H__
#define __PBP_H__

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>

using namespace std;
using std::stringstream;
using std::string;

class PBP
{
public:
		
    /*
     * Auxiliary vector for inversion.
     */
    int * m_aux;
    
	int* identityPerm;
	//int** identityPermutations;
	vector<int*> identityPermutations;
	int nChanges, pbsize;
	double* identityPermutationChanges; // percentage of the maxFes at which a change happens
	int nextChangeIndex;
	string dynProfilePath;

	/*
	 * The constructor.
	 */
	PBP();
	
	/*
	 * The destructor. It frees the memory allocated.
	 */
	virtual ~PBP();
	
	/*
	 * Virtual evaluation function.
	 */
	virtual double Evaluate(int*genes)=0;

    /*
	 * Virtual evaluation function of the inverted solution.
	 */
	virtual double EvaluateInv(int*genes)=0;
    
    /*
	 * Virtual instance reading function.
	 */
	virtual int Read(string filename)=0;
    
    /*
     * Returns the size of the problem.
     */
    virtual int GetProblemSize()=0;
    

	virtual void setIdentityPermutationChanges(string dynamicfilename);

	virtual bool changeIdentityPermutation(int fes, int maxfes);

	virtual void setPbsize(int size);

private:
    

};
#endif
