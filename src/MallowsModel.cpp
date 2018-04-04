//
//  MallowsModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "MallowsModel.h"
#include "Kendall.h"
#include "Cayley.h"
#include "Ulam.h"
#include "Ulam2.h"
#include "Tools.h"
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
//#include <unistd.h>
#include <list>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;


/*
 * Class constructor.
 */
CMallowsModel::CMallowsModel(int problem_size, int sel_size, char * metric_type)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    if (((string)metric_type)=="K")
        m_distance_model= new Kendall_Model(problem_size,sel_size);
    else if (((string)metric_type)=="C"){
        m_distance_model= new Cayley_Model(problem_size,sel_size);
    }
    else{
       // m_distance_model = new Ulam_Model(problem_size,sel_size); //implementacion con el bratelli.
        m_distance_model = new Ulam_Model2(problem_size,sel_size); //implementacion con ficheros de datos.
    }
 
}

/*
 * Class destructor.
 */
CMallowsModel::~CMallowsModel()
{
	delete m_distance_model;
    delete [] m_samples;
}

/*
 * Virtual learning function.
 */
bool CMallowsModel::Learn(CPopulation * population, int size)
{
    for(int i=0;i<size;i++){
		m_samples[i] = population->m_individuals[i]->Genes();
    }
    return m_distance_model->Learn(m_samples, size);
}

/*
 * Given an infinite population of samples, it learns a Mallows model.
 */
bool CMallowsModel::Learn(CInfinitePopulation * population, float * selection_probs){
    
    for(int i=0;i<population->m_pop_size;i++){
		m_samples[i] = population->m_individuals[i]->Genes();
    }
    return m_distance_model->Learn(m_samples, selection_probs, population->m_pop_size);
}

/*
 * Builds the Mallows model with the given CR and theta parameters.
 */
bool CMallowsModel::Learn(int * consensus_ranking, double theta){
    return m_distance_model->Learn(consensus_ranking,theta);
}

/*
 * Virtual sampling function.
 */
void CMallowsModel::Sample(int * genes)
{
    m_distance_model->Sample(genes);
}

/*
 * Calculates the probability of the solution given the probabilistic model.
 */
float CMallowsModel::Probability(int * solution)
{
    return m_distance_model->Probability(solution);
}

