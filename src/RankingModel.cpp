//
//  RankingModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "RankingModel.h"
#include "Population.h"


/*
 * The constructor.
 */
CRankingModel::CRankingModel()
{
	
}

/*
 * The destructor.
 */
CRankingModel::~CRankingModel()
{

}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(CPopulation * population, int size)
{
    return true;
}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(CInfinitePopulation * population, float * selection_probs)
{
    return true;
}


/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(int * consensus_ranking, double theta_parameter)
{
    return true;
}


/*
 * Virtual sampling function.
 */
void CRankingModel::Sample(int * genes)
{
    
}

/*
 * Calculates the probability of the solution given the probabilistic model.
 */
float CRankingModel::Probability(int * solution)
{
    return 0;
}