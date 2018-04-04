/*
 *  Cayley.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/18/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#ifndef __CAYLEY_MODEL_H__
#define __CAYLEY_MODEL_H__

#include "Distance.h"

#define ALF 1.0e-4
#define TOLX 1.0e-7
#define MAXITS 200
#define EPS 1.0e-4
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define STPMX 100.0
#include <cmath>
#include <iostream>
#include <math.h>


#define FREERETURN {free_matrix(fjac,1,n,1,n);free_vector(fvec,1,n);\
free_vector(p,1,n);free_ivector(indx,1,n);return;}

typedef struct cayley_data
{
    int n;
    float dist_avg;
    int j;
    
}cayley_data;

class Cayley_Model : public Distance_Model
{
public:

	/*
	 * The constructor.
	 */
	Cayley_Model(int problem_size, int sel_size);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the Cayley model.
	 */
	virtual ~Cayley_Model();
	
	/*
     * Learns a Mallows model based on the Cayley distance and the individuals in the model.
     */
	bool Learn(int ** samples, int size);
	
    /*
	 * Learns a Mallows model based on the Cayley distance and the individuals in the model.
	 */
    bool Learn(int ** samples, float * probabilities, int size);
    
    /*
     * Builds the Mallows model for the Cayley distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
    
	/*
	 * Given the consensus ranking, it samples a new individual.
	 */
	void Sample(int * permutation);
	
	/*
	 * Calculates the probability of the given individual in the learnt Mallows probabilistic model with the Cayley distance.
	 */
	double Probability(int * individual);
    
    /*
     * Determines if the given indiviual is a local optima for the problem and the Cayley distance.
     */
    bool isLocalOptima(int * indiviual, PBP * problem);
	
private:
    
	/*
	 * Xj probabilities matrix.
	 */
	float* m_xprobs;
	
    /*
     * Stirling matrix.
     */
    long double **stirling_matrix_;

    /*
     * Auxiliary vector for composition
     */
    int * m_composed;
    
    /*
     * Auxiliary vector for inversion.
     */
    int * m_inverted;

    /*
     * Auxiliary vector for the x.
     */
    int * m_x;

    /*
     * Auxiliary vector for the x.
     */
    int * m_aux;
    
    /*
     * Auxiliary variable for the distance;
     */
    bool * m_visited;
    
    /*
     * Data structure for the generic input of parameters into NewtonRaphson.
     */
    cayley_data * m_newton_d;
    
	/*
	 * Calculates de consensus permutation from the given population cases.
	 */
    void CalculateConsensusRanking(int ** samples, int cases_num, int*consensusPermutation);
	
    /*
     * Calculates de consensus permutation from the given population cases.
     */
    void CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, float * weights, int* consensus_ranking);
    
	/*
	 * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
	 */
    float CalculateThetaParameter(int*consensus, int ** samples, int cases_num);
    
    /*
     * Calculates the spread theta parameter from the m_consensus_ranking and the individuals in the population.
     */
    float CalculateThetaParameter_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, float * weights);
    
    /*
	 * Calculates the total Psi normalization constant from the ThetaParameters and psi-s vector.
	 */
	float CalculatePsiConstants(float theta, float* psi);
	
    /*
     * Calculates the sum of Cayley distances of the given solution to the sample of solutions.
     */
    float CalculateDistancetoSample(int * solution, int ** samples, int cases_num);
    
    /*
     * Calculates the sum of Cayley distances of the given solution to the sample of solutions.
     */
    float CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, float * weights);
    
    //////////////////////////////////////////OBTAIN PERMUTATION FROM A X VECTOR//////////////////////////////////
    /*
     * Generates a permutation from a x vector.
     */
    void GeneratePermuFromX(int *x, int *permu);

    inline bool item_closes_cycle(int pos, int item, int *sigma, int *sigma_inv);

    /*
     * Obtains the most probable cycle when building the permutation from a given vector of X.
     */
    int GetMostProbCycle(int ind, int **cycles, int len, int *leng_cycles);
    
    /*
     * Generates a random permutaiton with K cycles.
     */
    int * GeneratePermutation_withKcycles(int n, int k);
    
    /*
     * Calculates the distance of the given solution with respect to the identity permutation
     */
    int CalculateDistanceAndX(int * sigma, int *x);
    
    //////////////////////////////////////////ESTIMATION OF THETA PARAMETER//////////////////////////////////
    
    /*
     * Theta parameter estimation function.
     */
    double f(double theta, double dist_avg_);
    
    /*
     * Theta parameter estimation function derivation.
     */
    double fdev(double theta);
    
    //////////////////////////////////////////ESTIMATION OF CONSENSUS RANKING//////////////////////////////////
    
    void Approximated_MLE_Consensus( int m, int **samples, int *sigma_0, double *best_likeli);
    void approx_consensus_mm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_distance);
    void variable_neighborhood_search(int m, int **samples, int *sigma, double *f_eval);
    void local_search_swap_mm(int m, int **samples, int *sigma_0, double *f_eval);
    void local_search_insert (int m, int **samples, int *sigma_0, double *f_eval);
    int get_cycles(int *sigma, int *cycle_items, int *cycle_indices);
    int distance_to_sample(int m, int **samples, int *sigma);
};

#endif

