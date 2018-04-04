//
//  Kendall.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef _KENDALL_
#define _KENDALL_
#include "Distance.h"
#include "NewtonRaphson.h"


typedef struct kendall_data
{
    int n;
    float * Vjsmean;
    int j;
    
}kendall_data;


class Kendall_Model : public Distance_Model
{
public:

	/*
	 * The constructor.
	 */
	Kendall_Model(int problem_size, int sel_size);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~Kendall_Model();
	
	/*
	 * Learns the Mallows kendall model give a population of samples.
	 */
    bool Learn(int ** samples, int size);
    
    /*
     * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
	
    /*
	 * Learns a Mallows model based under the Kendall distance for the given infinite individual.
	 */
	bool Learn(int ** samples, float * probabilities, int size);
    
	/*
	 * Given the consensus ranking, it samples a new individual.
	 */
	void Sample(int * permutation);
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	double Probability(int * individual);
	
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
     * Auxiliary vector used for calculating the consensus ranking.
     */
    float * consensusVector;
    
    /*
     * Auxiliary data structures for learning stage.
     */
    float * VjsMean;
	float * VjsNonMean;
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
     * Calculates de consensus permutation from the given population cases.
     */
    void CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, float * weights, int* consensus_ranking);
    
    /*
     * Calculates the sum of Kendall distances of the given solution to the sample of solutions.
     */
    float CalculateDistancetoSample(int * solution, int ** samples, int cases_num);
    
    /*
     * Calculates the sum of Kendall distances of the given solution to the sample of solutions.
     */
    float CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, float * weights);
    
    /*
     * Calculates the Kendall distance between two permutations.
     */
    int CalculateKendallDistance(int * perm_a, int * perm_b);

    /*
     * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
     */
    float CalculateThetaParameter(int*consensus_ranking, int** samples, int samples_num);
    
    /*
     * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
     */
    float CalculateThetaParameter_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, float * weights);
    
    /*
     * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
     */
    void CalculatePsiConstants(float theta, float* psis);
    
	/*
	 * Generates a permutation from a v vector.
	 */
	void GeneratePermuFromV(int * v, int * permu);

    /*
     * Determines if the given indiviual is a local optima for the problem and the current distance.
     */
    bool isLocalOptima(int * indiviual, PBP * problem);
	
};

#endif
