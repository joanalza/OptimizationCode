//
//  MallowsModel.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDA__MallowsModel__
#define __RankingEDA__MallowsModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include "Distance.h"
#include <list>
#include <vector>


class CMallowsModel : public CRankingModel
{
	
public:
	
    /*
     * Distance model.
     */
    Distance_Model * m_distance_model;
    
    /*
     * The constructor.
     */
	CMallowsModel(int problem_size, int sel_size, char * metric_type);
	
    /*
     * The destructor.
     */
    virtual ~CMallowsModel();
	
    /*
     * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
    
    /*
     * Given an infinite population of samples, it learns a Mallows model.
     */
    bool Learn(CInfinitePopulation * population, float * selection_probs);
    
    /*
     * Given a population of samples, it learns a Mallows model from the data.
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

#endif /* defined(__RankingEDA__MallowsModel__) */
