//
//  InfiniteRankingEDA.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 14/03/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__InfiniteRankingEDA__
#define __RankingEDAsCEC__InfiniteRankingEDA__

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "PBP.h"
#include "Tools.h"
#include "InfinitePopulation.h"
#include "RankingModel.h"

using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

class InfiniteRankingEDA
{
	
public:
	
    /*
     * Problem
     */
    PBP * m_problem;
    
    /*
     * Problem size
     */
    int m_problem_size;
    
    /*
     * Maximum number of evaluations permitted in the execution.
     */
    long int m_max_evaluations;
    
    /*
     * The number of evaluations performed at the moment.
     */
    int m_evaluations;

    /*
     * The size of the population.
     */
    int m_pop_size;
    
    //The type of the metric to use in the model of the EDA.
    char m_metric_type[10];
    
    //The type of the model used in the EDA.
    char m_model_type[10];
    
   
    /*
     * The population
     */
    CInfinitePopulation * m_population;
    
    /*
     * The probalistic model.
     */
    CRankingModel * m_model;
    
    /*
     * The constructor.
     */
	InfiniteRankingEDA(PBP * problem, int problem_size, long int max_evaluations, char * model_type, char * metric_type, int inverse);
	
    /*
     * The destructor.
     */
    virtual ~InfiniteRankingEDA();
	
    /*
     * Running function
     */
	int Run();
    
    
private:
    
	
};
#endif
