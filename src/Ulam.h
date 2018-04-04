//
//  Ulam.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/25/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__Ulam__
#define __RankingEDAsCEC__Ulam__
#include "Ferrers_diagram.h"
#include "Distance.h"
#include <iostream>
//return values of generation functions
#define GEN_NEXT  0 //ok, print and continue
#define GEN_TERM  1 //ok, terminate
#define GEN_EMPTY 2 //ok, print EMPTY SET and continue
#define GEN_ERROR 3 //an error occured, print an error message and terminate

typedef struct ulam_data
{
    int n;
    float dist_avg;
    long double * count;
    
}ulam_data;

class Ulam_Model : public Distance_Model
{
public:
    
	/*
	 * The constructor.
	 */
	Ulam_Model(int problem_size, int sel_size);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the ulam model.
	 */
	virtual ~Ulam_Model();
	
    /*
     * Learns a Mallows model based on the Ulam distance and the individuals in the model.
     */
	bool Learn(int ** samples, int size);
	
    /*
	 * Calculates de consensus permutation from the given population cases.
	 */
    int CalculateConsensusRanking(int ** samples, int cases_num, int*consensusPermutation);
	
	/*
	 * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
	 */
    float CalculateThetaParameter(int*consensus, int ** samples, int cases_num, float upper_bound, int dist_avg);
    
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
	
    float m_averageDistance;
private:
    /*
     * Vector of ferrer diagrams
     */
    std::vector<Ferrers_diagram*> * m_bratteli;
    
    
    long double * m_count_permus;
    double * m_probabilities;
    int m_bound_lis;
    long double*m_acumul_ld;
    double * m_acumul;
    int m_num_partitions;
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
    ulam_data * m_newton_d;

	//////////////////////////////////////////ESTIMATION OF CONSENSUS RANKING//////////////////////////////////

    /*
     * Calculates the sum of Ulam distances of the given solution to the sample of solutions.
     */
    int CalculateDistancetoSample(int * solution, int ** samples, int cases_num);

    
    //////////////////////////////////////////ESTIMATION OF THETA PARAMETER//////////////////////////////////
    
    void fill_bratteli_diagram(int bound_lis);
    void insert_shapes_bound(int current_n, bool has_parents, int bound_from, int bound_to);
    int gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k);
    int gen_part_next(unsigned char *vector, unsigned char *k, int bound);
    
    ////////////////////////////////////////SAMPLING///////////////////////////////
    
    int*m_col_index;
    int*m_row_index;
    
    /*
     * Given a LIS size, generates a permu at that distance.
     */
    void GeneratePermuAtLIS(int l, int *sigma);
    
    /*
     * Generates a random Young Standard Tableaux given a ferrer diagram.
     */
    void GenerateRandomYST(Ferrers_diagram *dad, int **out, int*row_index, int*col_index);
    
    /*
     * Prints the Bratteli diagram.
     */
    void print_bratteli_diagram();
};


#endif /* defined(__RankingEDAsCEC__Ulam__) */
