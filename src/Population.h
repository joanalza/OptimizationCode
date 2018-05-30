//
//  Population.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef _POPULATION_
#define _POPULATION_

#include <vector>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "Individual.h"

using std::istream;
using std::ostream;
using std::stringstream;
using std::string;
using namespace std;
class CPopulation;


class CPopulation
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
     * Size of the vector of individuals
     */
    int m_size;
    
    /*
     * Size of the population
     */
    int m_pop_size;
    
    /*
     * Size of the population of offsprings
     */
    int m_offspring_size;

	/*
	 * Sum of all the fitness of the population
	 */
    
	long int m_all_fitness;

	/*
     * The constructor. It creates an empty list.
     */
    CPopulation(int pop_size, int offspring_size, int individual_size);
    
	/*
     * The destructor.
     */
	virtual ~CPopulation();
    
	/*
	 * Function to add an individual to the population
	 */
    void AddToPopulation(int * genes, int index, long int fitness);
    
    /*
     * Function to set an individual in a specific position of the population
     */
    void SetToPopulation(int * genes, int index, long int fitness);
    
    /*
     * Sorts the individuals in the population in decreasing order of the fitness value.
     * mode = 0 means sort only the population.
     * mode = 1 means sort population and offsprings population together.
     */
    void SortPopulation(int mode);
    
	/*
	 * Prints the current population.
	 */
	void Print();
	
	/*
	 * Determines if the individuals in the population are equal.
	 */
	bool Same(int size);
    
    /*
     * Calculates the average fitness of the first 'size' solutions in the population
     */
    float AverageFitnessPopulation();

private:
    
};

#endif


