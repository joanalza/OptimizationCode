//
//  GeneralizedMallowsModel.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__GeneralizedMallowsModel__
#define __RankingEDAsCEC__GeneralizedMallowsModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include "GeneralizedDistance.h"
#include <list>
#include <vector>


class CGeneralizedMallowsModel : public CRankingModel
{
	
public:
	
    /*
     * Distance model.
     */
    GeneralizedDistance_Model * m_distance_model;
    
    
    /*
     * The constructor.
     */
	CGeneralizedMallowsModel(int problem_size, int sel_size, char * metric_type);
	
    /*
     * The destructor.
     */
    virtual ~CGeneralizedMallowsModel();
	
    /*
     * Given a population of samples, it learns a Generalized Mallows model from the data.
     */
    bool Learn(CPopulation * population, int size);
    
    /*
     * Creates a new individual sampling the generalized mallows probabilistic model.
     */
    CIndividual * Simulate();
    
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	float Probability(int * individual);
    
    /*
     * Given the consensus ranking, it samples a new individual given the s
     */
    void Sample(int * permutation);
	
private:
    
    /*
     * Problem size.
     */
    int m_problem_size;
    
    /*
     * Sample size.
     */
    int m_sample_size;
    
    /*
     * Sample of solutions.
     */
    int ** m_samples;
    
    
	
};
#endif /* defined(__RankingEDAsCEC__GeneralizedMallowsModel__) */
