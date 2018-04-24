/*
 *  PlackettLuce.h
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 7/16/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */


#ifndef PLACKETTLUCE_H__
#define PLACKETTLUCE_H__

#include <list>
#include <vector>
#include "RankingModel.h"
#include "Population.h"
#include "PBP.h"
class CPlackettLuceModel : public CRankingModel
{
public:
    
	/*
	 * The constructor.
	 */
	CPlackettLuceModel(int problem_size, int sel_size);
	
    /*
	 * The destructor.
	 */
    virtual ~CPlackettLuceModel();
	
    /*
     * Given a population of samples, it learns a Plackett-Luce model from the data.
     */
    bool Learn(CPopulation * population, int size);
    
    /*
     * Learns a Mallows model under the Ulam distance from a boltzmann distribution, and returns the full Mallows distribution
     */
    double * Learn_fromBoltzmann(int ** search_space, double * boltzmann_probabilities, int size);
    
	/*
	 * Simulates a new individual sampling the probabilistic model.
	 */
    int Sample(CPopulation * population, int num_samples, int inverse, PBP * problem);
    
	/*
     * Calculates ThetaParameters average value.
     */
	double GetWeightsAverage();
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	double Probability(int * individual);
    
    /*
     * This method updates the best solution of the optimization processes if the most probable solutions
     * are better than those solutions obtained by the population.
     */
    void PostProcesses();
    
private:
    
    /*
     * Samples individuals from the distribution
     */
    void Sample_BoltzmannDistribution(int ** search_space,double * distribution, int size, int sample_size);

    /*
     * It applies random shake_power times a perturbation over the given individual.
     */
    void ShakeByInsert(int * permutation, int shake_power);
    
    /*
	 * It estimates the Plackett-Luce model from the given cases with the MM algorithm.
	 */
	bool LearnMM(CPopulation * population, int sel_total);
    
    /*
     * It estimates the Plackett-Luce model from the given Boltzmann distribution over the search space with the MM algorithm.
     */
    bool LearnMM_fromBoltzmann(int ** search_space, double * boltzmann_probabilities, int size);

    /*
     * Samples the plackett-luce model learnt by shaking the most probable solution.
     */
    void SampleByShake(int * permutation, int * most_probable);

    
    CIndividual * SimulateInverse();
    
    /*
     * Calculates the most probable solution given the vector of weights
     */
    void CalculateMostProbableSolution();

    /*
     * Returns the most probable solution given the vector of weights
     */
    CIndividual * GetMostProbableSolution();
    
    int ** m_huge_sample;
	
	/*
     * Problem size.
     */
    int m_problem_size;
    int m_sample_size_huge;

    /*
     * Number of samples.
     */
    int m_samples_size;
    
    /*
     * Number of maximum iterations to perform by the MM algorithm.
     */
    int m_MM_iterations;

	/*
     * Parameter vector v
     */
	double * m_weights;
    	
    /* 
     * Euclidean norm of the weights.
     */
    double m_norm;
    
    /*
     * Most probable individual according to the straight model.
     */
    CIndividual * m_mostProbable;
    
    
    /*
     * The auxiliary vector for sampling solutions.
     */
    int * m_sampling_permutation;

    /*
     * Auxiliary data.
     */
    int * m_aux1;
    int * m_aux2;
    
    // Number of items fixed in the template.
    int m_k;
    
    //Parameters for learning
    int m_M;
	int m_N; //number of rankings in the population.
	int m_P;
    int ** f;
    int ** r;
    double ** r2;
    double ** g;
    double * dgamma;
	double * gamma;
    double * w;
	int * pp;
    double **aux;
    double * newgamma;
    double * sr2;
    int * m_last;
    int m_estimation_type;
    
    //parameters for sampling
    double * distribution;
	int * objects;
    
    //forward-backward sampling arrays
    int * m_located;
    int * m_p1;
    int * m_p2;

};
#endif
