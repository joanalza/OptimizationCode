//
//  Distances.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __DISTANCE_MODEL_H__
#define __DISTANCE_MODEL_H__
#include "InfinitePopulation.h"

class Distance_Model
{
public:
	
	/*
	 * Problem size.
	 */
	int m_problem_size;
	
    /*
     * Sample size.
     */
    int m_sample_size;
    
	/*
	 * Spread parameters vector.
	 */
	float m_theta_parameter;
		
	/*
	 * The consensus ranking.
	 */
	int * m_consensus_ranking;
	
	/*
	 * Psi normalization constant.
	 */
	float * m_psis;
	
    /*
     * Defines upper theta bound.
     */
    float m_upper_theta_bound;
    
    /*
     * Defines lower theta bound.
     */
    float m_lower_theta_bound;
    
	/*
	 * The constructor.
	 */
	Distance_Model();
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~Distance_Model();

    /*
	 * Learns a Mallows model based on the Kendall distance and the individuals in the model.
	 */
	virtual bool Learn(int ** samples, int size)=0;
    
    /*
	 * Learns a Mallows model based on the Kendall distance and the individuals in the model.
	 */
	virtual bool Learn(int ** samples, float * probabilities, int size)=0;
    
    /*
	 * Learns a Mallows model based on the Kendall distance and consensus ranking and theta parameters given.
	 */
	virtual bool Learn(int * consensus_ranking, double theta)=0;
    
    /*
	 * Given the consensus ranking, it samples a new individual.
	 */
	virtual void Sample(int * permutation)=0;
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	virtual double Probability(int * individual)=0;

    /*
     * Determines if the given indiviual is a local optima for the problem and the current distance.
     */
    virtual bool isLocalOptima(int * indiviual, PBP * problem)=0;
private:
		

};
#endif