//
//  InfinitePopulation.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 14/03/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#include "InfinitePopulation.h"
#include "Tools.h"

/*
 * Constructor function.
 */
CInfinitePopulation::CInfinitePopulation(int pop_size, int individual_size, PBP * problem)
{
    //1. Initialize infinite population with empty solutions
    m_pop_size=pop_size;
    m_individuals.resize(m_pop_size);
    for (int i=0;i<m_pop_size;i++)
        m_individuals[i]= new CIndividual(individual_size);
    m_probabilities= new float[m_pop_size];
    
    //2. Build initial population. All the solutions of the search space.
    
	//int solution[individual_size];
	int* solution = new int[individual_size];
	
    for (int i=0;i<individual_size;i++) solution[i]=i;
    sort(solution,solution+individual_size);
    int index=0;
    do{
        SetToPopulation(solution, index, problem->EvaluateInv(solution));
        index++;
    }
  	while ( next_permutation (solution,solution+individual_size) );
    SortPopulation();

    //3. Generate initial random vector of probabilities.
    double acumul=0;
    for (int i=0;i<m_pop_size;i++){
        m_probabilities[i]=(double)(rand()%100);
        acumul+=m_probabilities[i];
    }
    acumul=acumul-m_probabilities[m_pop_size-1];
    m_probabilities[m_pop_size-1]=200000;
    acumul+= m_probabilities[m_pop_size-1];
    
    for (int i=0;i<m_pop_size;i++)
        m_probabilities[i]=m_probabilities[i]/acumul;
    
    cout<<"initialized "<<m_pop_size<<" solutions."<<endl;
}

/*
 * Destructor function.
 */
CInfinitePopulation::~CInfinitePopulation()
{
    for (int i=0;i<m_pop_size;i++)
    {
        delete m_individuals[i];
    }
    m_individuals.clear();
    
    delete [] m_probabilities;
}

/*
 * Function to set an individual in a specific position of the population
 */
void CInfinitePopulation::SetToPopulation(int * genes, int index, long int fitness)
{
    m_individuals[index]->SetGenes(genes);
    m_individuals[index]->SetValue(fitness);
}

/*
 * Prints the current population.
 */
void CInfinitePopulation::Print()
{
	for(int i=0;i<m_pop_size; i++)
		cout<<m_probabilities[i]<<"     [ "<<m_individuals[i]<<" ]"<<endl;
}

/*
 * Sorts the individuals in the population in decreasing order of the fitness value.
 */
void CInfinitePopulation::SortPopulation()
{
    sort(m_individuals.begin(), m_individuals.end(), Better);
}
