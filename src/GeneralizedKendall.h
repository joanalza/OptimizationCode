//
//  GeneralizedKendall.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__GeneralizedKendall__
#define __RankingEDAsCEC__GeneralizedKendall__

#include "GeneralizedDistance.h"


typedef struct kendall_data
{
    int n;
    float * Vjsmean;
    int j;
    
}kendall_data;


class GeneralizedKendall_Model : public GeneralizedDistance_Model
{
public:
    
	/*
	 * The constructor.
	 */
	GeneralizedKendall_Model(int problem_size);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~GeneralizedKendall_Model();
	
	/*
	 * Learns the Mallows kendall model give a population of samples.
	 */
    bool Learn(int ** samples, int size);
	
    /*
     * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
    
	/*
	 * Given the consensus ranking, it samples a new individual.
	 */
	void Sample(int * permutation);
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	float Probability(int * individual);
	
private:
	
    /*
     * Matrix of V probabilities.
     */
	float ** m_vprobs;
    
    /*
     * Frequency matrix employed when calculating the consensus ranking.
     */
    int ** m_frequency_matrix;
    
    /*
     * Auxiliary data structures for sampling stage.
     */
    int * aux;
    int * aux_v;
    int * aux_n;
    
    /*
     *
     */
    float * consensusVector;
    
    /*
     * Auxiliary data structures for learning stage.
     */
    float * VjsMean;
	int * VjsNonMean;
	int * Vjs;
	int * invertedB;
	int * composition;
    
    /*
     * Data structure for the generic input of parameters into NewtonRaphson.
     */
    kendall_data * m_newton_d;    
    
	/*
     * Calculates de consensus permutation from the given population cases.
     */
	void CalculateConsensusRanking(int** samples, int sample_size, int* consensus_ranking);
	
    /*
     * Calculates the distances of the given solution to the sample.
     */
    float CalculateDistancetoSample(int * solution, int ** samples, int cases_num);

    /*
     * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
     */
    void CalculateThetaParameters(int*consensus_ranking, int** samples, int samples_num,float * theta_parameters);
    
    /*
     * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
     */
    void CalculatePsiConstants(float* thetas, float* psi);
    
	/*
	 * Generates a permutation from a v vector.
	 */
	void GeneratePermuFromV(int * v, int * permu);
    	
};

#endif /* defined(__RankingEDAsCEC__GeneralizedKendall__) */
