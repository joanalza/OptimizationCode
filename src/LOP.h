/*
 *  LOP.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _LOP_H__
#define _LOP_H__

#include "Tools.h"
#include "PBP.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>

using std::ifstream;
using std::ofstream;
using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

class LOP : public PBP
{
	
public:
	
    /*
     * Entries matrix of the LOP.
     */
	int ** m_matrix;
    
	/*
	 * The size of the problem.
	 */
	int m_problemsize;
		

    /*
     * The constructor.
     */
	LOP();
	
    /*
     * The destructor.
     */
    virtual ~LOP();
	
	/*
	 * Read LOP instance file.
	 */
	int Read(string filename);
	
	/*
	 * This function evaluates the solution for the LOP problem.
	 */
	double Evaluate(int * genes);
    
    /*
     * This function evaluates the inverted solution of the given solution for the LOP problem.
     */
    double EvaluateInv(int * genes);
	
    /*
     * Returns the size of the problem.
     */
    int GetProblemSize();
    
private:
    
};
#endif





