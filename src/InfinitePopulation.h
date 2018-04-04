//
//  InfinitePopulation.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 14/03/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__InfinitePopulation__
#define __RankingEDAsCEC__InfinitePopulation__

#include <vector>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "Individual.h"
#include "PBP.h"

using std::istream;
using std::ostream;
using std::stringstream;
using std::string;
using namespace std;
class CInfinitePopulation;


class CInfinitePopulation
{
    
public:
    
    struct Better{
        bool operator()(CIndividual * a, CIndividual * b) const{
            return a->Value()>b->Value();
        }
    } Better;
    
    /*
     * Vector of individuals that constitute the population.
     */
    vector<CIndividual *> m_individuals;
    
    /*
     * Array of probabilities associated to the solutions in the infinite population,
     */
    float * m_probabilities;
    
    /*
     * Size of the population
     */
    int m_pop_size;
    
    /*
     * Constructor function.
     */
    CInfinitePopulation(int pop_size, int individual_size, PBP * problem);
    
	/*
     * The destructor.
     */
	virtual ~CInfinitePopulation();
    
    /*
     * Function to set an individual in a specific position of the population
     */
    void SetToPopulation(int * genes, int index, long int fitness);
    
    /*
     * Sorts the individuals in the population in decreasing order of the fitness value.
     */
    void SortPopulation();
    
	/*
	 * Prints the current population.
	 */
	void Print();
    
private:
    
};

#endif /* defined(__RankingEDAsCEC__InfinitePopulation__) */
