//
//  GeneralizedDistances.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__GeneralizedDistances__
#define __RankingEDAsCEC__GeneralizedDistances__

class GeneralizedDistance_Model
{
public:
	
	/*
	 * Problem size.
	 */
	int m_problem_size;
	    
	/*
	 * Spread parameters vector.
	 */
	float * m_theta_parameters;
    
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
	GeneralizedDistance_Model();
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~GeneralizedDistance_Model();
    
	/*
	 * Learns a Mallows model based on the Kendall distance and the individuals in the model.
	 */
	virtual bool Learn(int ** samples, int size)=0;
        
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
	virtual float Probability(int * individual)=0;
    
private:
    
    
};


#endif /* defined(__RankingEDAsCEC__GeneralizedDistances__) */
