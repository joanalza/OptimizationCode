//
//  GeneralizedMallowsModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "GeneralizedMallowsModel.h"

#include "GeneralizedKendall.h"
#include "GeneralizedCayley.h"
//#include "Cayley.h"
//#include "Ulam.h"
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
CGeneralizedMallowsModel::CGeneralizedMallowsModel(int problem_size, int sel_size, char * metric_type)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
    m_samples= new int*[m_sample_size];

    if (((string)metric_type)=="K"){
        m_distance_model= new GeneralizedKendall_Model(problem_size);
    }
    else if (((string)metric_type)=="C"){

         m_distance_model= new GeneralizedCayley_Model(problem_size);
    }
    else{
        cout<<"There is no Generalized Mallows implementation for the metric "<<(string)metric_type<<"."<<endl;
    }
    
}

/*
 * Class destructor.
 */
CGeneralizedMallowsModel::~CGeneralizedMallowsModel()
{
	delete m_distance_model;
    delete [] m_samples;
}

/*
 * Virtual learning function.
 */
bool CGeneralizedMallowsModel::Learn(CPopulation * population, int size)
{

    for(int i=0;i<size;i++){
		m_samples[i] = population->m_individuals[i]->Genes();
    }

    return m_distance_model->Learn(m_samples, size);
}

/*
 * Virtual sampling function.
 */
void CGeneralizedMallowsModel::Sample(int * genes)
{
    m_distance_model->Sample(genes);
}

/*
 * Calculates the probability of the solution given the probabilistic model.
 */
float CGeneralizedMallowsModel::Probability(int * solution)
{
    return m_distance_model->Probability(solution);
}
