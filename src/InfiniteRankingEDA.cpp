//
//  InfiniteRankingEDA.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 14/03/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#include "InfiniteRankingEDA.h"
#include "MallowsModel.h"
#include "Cayley.h"
#include "Kendall.h"
#include "GeneralizedMallowsModel.h"
#include <iomanip>

/*
 * The constructor.
 */
InfiniteRankingEDA::InfiniteRankingEDA(PBP * problem, int problem_size, long int max_evaluations, char* model_type, char* metric_type, int inverse)
{
    //1. standard initializations
    m_problem=problem;
    m_problem_size=problem_size;
    m_max_evaluations=max_evaluations;
    m_evaluations=0;
    m_pop_size=factorial(m_problem_size);
    memcpy(m_metric_type,metric_type,sizeof(char)*10);
    memcpy(m_model_type,model_type,sizeof(char)*10);
    
    //2. initialize the infinite population
    m_population= new CInfinitePopulation(m_pop_size, m_problem_size, m_problem);
    
    //3. Build model structures
    m_model=new CMallowsModel(m_problem_size, m_pop_size, metric_type);
}

/*
 * The destructor. It frees the memory allocated..
 */
InfiniteRankingEDA::~InfiniteRankingEDA()
{
    delete m_population;
    delete m_model;
}

/*
 * Running function
 */
int InfiniteRankingEDA::Run(){
    
    //cout<<"Running..."<<endl;
    //variable initializations.
    int iterations=1;
    int i=0;
    float * selection_prob= new float[m_pop_size];
    float prob_max=0;
    int index_best=0;
    for (i=0;i<m_pop_size;i++){
        if (m_population->m_probabilities[i]>prob_max){
            index_best=i;
            prob_max=m_population->m_probabilities[i];
        }
    }
    m_population->Print();
    cout<<"Initial: "<<m_population->m_individuals[index_best]<<" Prob: "<<m_population->m_probabilities[index_best]<<endl;

    cout<<endl<<"-----------------------------"<<endl;
        cout<<endl<<"-----------------------------"<<endl;
        cout<<endl<<"-----------------------------"<<endl;
    
    //EDA iteration. Stopping criterion is the maximum number of evaluations performed
    while (iterations<20 ){
        
        //apply selection probabilities
        double sumatory=0;
        for (i=m_pop_size-1;i>=0;i--){
            selection_prob[i]=pow(m_population->m_probabilities[i],2) + 2*m_population->m_probabilities[i]*sumatory;
            sumatory+=m_population->m_probabilities[i];
        }
        
        //learn model
        m_model->Learn(m_population,selection_prob);
        
        //recalculate probabilities
        for (i=0;i<m_pop_size;i++){
            m_population->m_probabilities[i]=m_model->Probability(m_population->m_individuals[i]->Genes());
        }
       //     m_population->Print();
        cout<<endl<<"-----------------------------"<<endl;
        iterations++;
    }
    prob_max=0;
    index_best=0;
    for (i=0;i<m_pop_size;i++){
        if (m_population->m_probabilities[i]>prob_max){
            index_best=i;
            prob_max=m_population->m_probabilities[i];
        }
    }
    cout<<"Best: "<<m_population->m_individuals[index_best]<<" Prob: "<<m_population->m_probabilities[index_best]<<endl;
    cout<<"fitness: "<<m_problem->EvaluateInv(m_population->m_individuals[index_best]->Genes())<<endl;
    if (((CMallowsModel*)m_model)->m_distance_model->isLocalOptima(m_population->m_individuals[index_best]->Genes(),m_problem)){
        cout<<"IS local optima!"<<endl;
    }
    else{
        cout<<"is NOT local optima!"<<endl;
    }
        
    return 0;
}



