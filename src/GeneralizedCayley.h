// aa
//  GeneralizedCayley.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/22/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__GeneralizedCayley__
#define __RankingEDAsCEC__GeneralizedCayley__


#include "GeneralizedDistance.h"

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

class GeneralizedCayley_Model : public GeneralizedDistance_Model
{
public:
    
	/*
	 * The constructor.
	 */
	GeneralizedCayley_Model(int problem_size);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~GeneralizedCayley_Model();
	
	/*
	 * Learns a model for the mixture with the given CR, theta parameter and population.
	 */
	bool Learn(int ** samples, int size);
	
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
	float Probability(int * individual);
	
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
    
    int *m_x_acumul;
    /*
     * Auxiliary variable for the distance;
     */
    bool * m_visited;
    
	/*
	 * Calculates de consensus permutation from the given population cases.
	 */
    void CalculateConsensusRanking(int ** samples, int cases_num, int*consensusPermutation);
	
	/*
	 * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
	 */
    void CalculateThetaParameter(int*consensus, int ** samples, int cases_num, float * thetas);
    
    /*
	 * Calculates the total Psi normalization constant from the ThetaParameters and psi-s vector.
	 */
	void CalculatePsiConstants(float * thetas, float* psi);
	
    /*
     * Calculates the sum of Cayley distances of the given solution to the sample of solutions.
     */
    float CalculateDistancetoSample(int * solution, int ** samples, int cases_num);
    
    //////////////////////////////////////////OBTAIN PERMUTATION FROM A X VECTOR//////////////////////////////////
    /*
     * Generates a permutation from a x vector.
     */
    void GeneratePermuFromX(int *x, int *permu, int n);
    
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
    float get_theta_log_likelihood(int m, int *x_acumul, int *x_acumul_variation, float *theta_estim);
    
    //////////////////////////////////////////ESTIMATION OF CONSENSUS RANKING//////////////////////////////////
    
    void Approximated_MLE_Consensus( int m, int **samples, int *sigma_0, double *best_likeli);
    void approx_consensus_gmm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_likeli);
};

#endif /* defined(__RankingEDAsCEC__GeneralizedCayley__) */
