//
//  Ulam2.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 03/02/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__Ulam2__
#define __RankingEDAsCEC__Ulam2__

#include "Ferrers_diagram2.h"
#include "Distance.h"
#include <iostream>
//return values of generation functions
#define GEN_NEXT  0 //ok, print and continue
#define GEN_TERM  1 //ok, terminate
#define GEN_EMPTY 2 //ok, print EMPTY SET and continue
#define GEN_ERROR 3 //an error occured, print an error message and terminate

typedef struct ulam_data2
{
    int n;
    float dist_avg;
    long double * count;
    
}ulam_data2;

class Ulam_Model2 : public Distance_Model
{
public:
    
	/*
	 * The constructor.
	 */
	Ulam_Model2(int problem_size, int sel_size);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the ulam model.
	 */
	virtual ~Ulam_Model2();
	
    /*
     * Learns a Mallows model based on the Ulam distance and the individuals in the model.
     */
	bool Learn(int ** samples, int size);
    
    /*
     * Builds the Mallows model for the Ulam distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
	
    /*
     * Learns a Mallows model based under the Ulam distance for the given infinite individual.
     */
    bool Learn(int ** samples, float * probabilities, int size);

    /*
	 * Calculates de consensus permutation from the given population cases.
	 */
    int CalculateConsensusRanking(int ** samples, int cases_num, int*consensusPermutation);
	
    /*
     * Calculates de consensus permutation from the given population cases.
     */
    float CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, float * weights, int* consensus_ranking);
        
	/*
	 * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
	 */
    float CalculateThetaParameter(int*consensus, int ** samples, int cases_num, int dist_avg);
    
    /*
     * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
     */
    float CalculateThetaParameter_WeightedSamples(int*consensus, int ** samples, int cases_num, float average_dist);
    
	/*
	 * Given the consensus ranking, it samples a new individual.
	 */
	void Sample(int * permutation);
    
	/*
	 * Given the consensus ranking, it samples a new individual.
	 */
	void Sample(int * permutation, bool print);
	
	/*
	 * Calculates the probability of the given individual in the learnt Mallows probabilistic model with the Cayley distance.
	 */
	double Probability(int * individual);
	
private:

    char m_str_base_path[500];
    
    /*
     * Indicates the number of permutations at each Ulam distance.
     */
    long double * m_num_perms_distance;
    
    /*
     * The average distance of the permutations in the sample to the consensus ranking.
     */
    float m_averageDistance;
    
    /*
     * Array of probabilities used to sample the new solutions.
     */
    double * m_probabilities;

    /*
     * Auxiliary vector for composition
     */
    int * m_composed;
    
    /*
     * Auxiliary vector for inversion.
     */
    int * m_inverted;
    
    /*
     * Data structure for the generic input of parameters into NewtonRaphson.
     */
    ulam_data2 * m_newton_d;
    
    /*
     * Reads from file the number of permutations to each Ulam distance for n size of permutations.
     */
    bool ReadPermusPerDistance(int n);
    
    /*
     * Reads...
     */
    void ReadPermusPerShape(int d, int bound, int * shape, int*shape_len);

    /*
     * Calculates the number of permutations at each distance Ulam, and saves into files.
     */
    void SaveCounts_toFile();
	//////////////////////////////////////////ESTIMATION OF CONSENSUS RANKING//////////////////////////////////
    
    /*
     * Calculates the sum of Ulam distances of the given solution to the sample of solutions.
     */
    int CalculateDistancetoSample(int * solution, int ** samples, int cases_num);
    
    /*
     * Calculates the sum of Ulam distances of the given solution to the sample of solutions.
     */
    float CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, float * weights);
    
    //////////////////////////////////////////ESTIMATION OF THETA PARAMETER//////////////////////////////////
    
    int gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k);
    int gen_part_next(unsigned char *vector, unsigned char *k, int bound);
    

    ////////////////////////////////////////SAMPLING///////////////////////////////
    
    int * m_shape;
    int * m_col_index;
    int * m_row_index;
    
    /*
     * Given a LIS size, generates a permu at that distance.
     */
    void GeneratePermuAtLIS(int l, int *sigma);
    
    /*
     * Determines if the given indiviual is a local optima for the problem and the current distance.
     */
    bool isLocalOptima(int * individual, PBP * problem);
    
};



#endif /* defined(__RankingEDAsCEC__Ulam2__) */
