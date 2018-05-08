//
//  GeneralizedKendall.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "GeneralizedKendall.h"
#include "Variables.h"
#include "Tools.h"
#include "NewtonRaphson.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <climits>
using std::cerr;
using std::cout;
using std::endl;

/*
 * The constructor of the model.
 */
GeneralizedKendall_Model::GeneralizedKendall_Model(int problem_size)
{

	m_problem_size=problem_size;
	m_theta_parameters=new float[m_problem_size-1];
	m_psis= new float[m_problem_size-1];
	m_consensus_ranking=new int[m_problem_size];
    
	//reservar memoria para las probabilidades de Vs.
	m_vprobs=new float*[m_problem_size-1];
	for (int i=0;i<m_problem_size-1;i++)
	{
		m_vprobs[i]=new float[m_problem_size];
	}
    
    
    m_frequency_matrix= new int*[m_problem_size];
    for (int i=0;i<m_problem_size;i++)
        m_frequency_matrix[i]=new int[m_problem_size];
    
    m_lower_theta_bound=0.001;
    m_upper_theta_bound=12;
    m_newton_d = (kendall_data *) malloc(sizeof(kendall_data));
    m_newton_d->n=m_problem_size;
    
    //sampling stage auxiliary data structures
    aux= new int[m_problem_size];
    aux_v= new int[m_problem_size]; //last position always 0.
    aux_n= new int [m_problem_size];
    //learning stage auxiliary data structures
    VjsMean= new float[m_problem_size-1];
	VjsNonMean= new int[m_problem_size-1];
	Vjs = new int[m_problem_size-1];
	invertedB = new int[m_problem_size];
	composition = new int[m_problem_size];
    consensusVector=new float[m_problem_size];
    
}

/*
 * The destructor of the model.
 */
GeneralizedKendall_Model::~GeneralizedKendall_Model()
{
    delete [] m_consensus_ranking;
    
    for (int i=0;i<m_problem_size-1;i++)
	{
		delete [] m_vprobs[i];
	}
	delete [] m_vprobs;
	delete [] m_psis;
    
    for (int i=0;i<m_problem_size;i++)
        delete [] m_frequency_matrix[i];
    delete [] m_frequency_matrix;
    delete [] m_theta_parameters;
    
    //sampling stage auxiliary data structures
    delete [] aux;
    delete [] aux_v;
    delete [] aux_n;
    
    //learning stage auxiliary data structures
	delete [] VjsMean;
	delete [] VjsNonMean;
	delete [] Vjs;
	delete [] invertedB;
	delete [] composition;
    
    delete[] consensusVector;
    delete m_newton_d;
}


/*
 * Calculates de consensus permutation from the given population cases.
 */
bool GeneralizedKendall_Model::Learn(int ** samples, int size)
{
    //1.- Calculate the consensus ranking/permutation.
    CalculateConsensusRanking(samples, size, m_consensus_ranking);
    
    //2.- Calculate theta parameters
	CalculateThetaParameters(m_consensus_ranking, samples, size, m_theta_parameters);

    
	//3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameters, m_psis);
	
	//4.- Calculate Vj probabilities matrix.
	float upper, lower;
	int j,r;
	for(j=0;j<m_problem_size-1;j++)
	{
		for(r=0;r<m_problem_size-j;r++)
		{ //v(j,i) proba dql v_j tome valor i
			upper=exp((-1) * r * m_theta_parameters[j]);
			lower=m_psis[j];
			m_vprobs[j][r] =  upper/lower;
		}
	}
    return true;
}

/*
 * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
 */
bool GeneralizedKendall_Model::Learn(int * consensus_ranking, double theta)
{
    //1.- Copy the consensus ranking
    memcpy(m_consensus_ranking, consensus_ranking, sizeof(int)*m_problem_size);
    
    //2.- Copy the theta parameter.
    for (int i=0;i<m_problem_size;i++){
        m_theta_parameters[i]= theta;
    }
    
	//3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameters, m_psis);
	
	//4.- Calculate Vj probabilities matrix.
	float upper, lower;
	int j,r;
	for(j=0;j<m_problem_size-1;j++)
	{
		for(r=0;r<m_problem_size-j;r++)
		{ //v(j,i) proba dql v_j tome valor i
			upper=exp((-1) * r * m_theta_parameters[j]);
			lower=m_psis[j];
			m_vprobs[j][r] =  upper/lower;
		}
	}
    return true;
    
    return true;
}

/*
 * Calculates the probability of a given individuals in the current model.
 */
float GeneralizedKendall_Model::Probability(int * solution)
{
    Invert(m_consensus_ranking, m_problem_size, invertedB);
	Compose(solution,invertedB,composition,m_problem_size);
	vVector_Fast(Vjs,composition,m_problem_size,aux_v);
	
    float exponent=0;
	float Psi=1;
	for (int i=0;i<m_problem_size-1;i++){
		exponent+=m_theta_parameters[i]*Vjs[i];
		Psi=Psi*m_psis[i];
	}
    
	float probability=exp(-exponent)/Psi;
	return probability;
}



/*
 * It samples a new individual.
 */
void GeneralizedKendall_Model::Sample(int * permutation)
{
    int i, index, limit;
	float randVal, acumul;
	
	for (i=0;i < m_problem_size;i++) aux[i]=-1;
	
	//cout << "2.1. Sample a Vj."<<endl;
	for(i=0;i<m_problem_size-1;i++)
	{
		//muestreo de las n-1 posiciones del vector V que definen la permutacion
		randVal=(float)rand()/((float)RAND_MAX+1);
		acumul=m_vprobs[i][0];
        index=0;
        limit=m_problem_size-1-i;
		for(;(index<limit && acumul<randVal);index++)
			acumul += m_vprobs[i][index+1];
		aux_v[i]=index;
	}
    aux_v[m_problem_size-1]=0;
	GeneratePermuFromV(aux_v,aux);
    Compose(aux,m_consensus_ranking,permutation,m_problem_size); //original
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void GeneralizedKendall_Model::CalculateConsensusRanking(int** samples, int sample_size, int* consensus_ranking)
{

    int option=1;
    if (option==0){
        int i, j;
        
        //Initialize matrix to zeros.
        for (i=0;i<m_problem_size;i++)
        {
            for (j=0;j<m_problem_size;j++)
                m_frequency_matrix[i][j]=0;
        }
        
        int geneValue, genePosition;
        
        //Fill frecuency matrix reviewing individuals in the population.
        //para cada permutación muestreada:
        for (i = 0; i < sample_size; i++)
        {
            for (j=0;j<m_problem_size;j++)
            {
                geneValue=j; //original
                genePosition=samples[i][j]; //original
                
                m_frequency_matrix[genePosition][geneValue]++;
            }
        }
        
        //cout<<"Calculate consensus vector: "<<endl;
        //Calculate consensus vector.
        
        int job,position;
        float positionMean;
        for (job=0;job<m_problem_size;job++)
        {
            positionMean=0;
            for (position=0;position<m_problem_size;position++)
            {
                positionMean=positionMean+position*m_frequency_matrix[position][job];
            }

            consensusVector[job]=positionMean;
        }
        
        RandomKeys(m_consensus_ranking, consensusVector, m_problem_size);
    }
    else{
        //Option 3. Choose the solution that minimizes the distance to the sample.
        int best_index=0;
        int best_distance= INT_MAX;
        int dist;
        for (int i=0;i<sample_size;i++){
            dist=CalculateDistancetoSample(samples[i], samples, sample_size);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
        memcpy(consensus_ranking, samples[best_index],sizeof(int)*m_problem_size);
    }
}

/*
 * Calculates sum of the Kendall distances of the solution to the sample of solutions.
 */
float GeneralizedKendall_Model::CalculateDistancetoSample(int * solution, int ** samples, int cases_num){
    
    float dist=0;
    int i,j;
    for ( i=0;i<cases_num;i++){
        Invert(samples[i],m_problem_size,invertedB);
        Compose(solution, invertedB, composition, m_problem_size);
        vVector_Fast(Vjs,composition,m_problem_size,aux_v);
        
        for (j = 0; j < m_problem_size-1; j++)
            dist += Vjs[j];
    }
    return dist;
}


/*
 * Theta parameter estimation function.
 */
float f_GeneralizedKendall(float theta, void *d )
{
    kendall_data * data= (kendall_data*)d;
	float oper1 = ( (data->n - data->j + 1) /(float)(exp((data->n - data->j + 1)*theta) - 1 ) );
	float oper2 = data->Vjsmean[data->j-1];
	float results= ( 1 /(float)( exp(theta) -1 ) ) -oper1 -oper2;
	
	return results;
}

/*
 * Theta parameter estimation function derivation.
 */
float fdev_GeneralizedKendall(float theta, void *d )
{
    kendall_data * data= (kendall_data*)d;
	float oper1=((-1)*exp(theta) / pow(exp(theta)-1, 2));
	float oper2 =( pow(data->n-data->j+1,2) * exp(theta*(data->n-data->j+1))) / (float)pow(exp(theta*(data->n-data->j+1))-1,2);
	if (isnan(oper2)==true)
	{
		oper2=0;
	}
    
	return oper1+oper2;
}

/*
 * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
 */
void GeneralizedKendall_Model::CalculateThetaParameters(int*consensus_ranking, int** samples, int samples_num,float * theta_parameters)
{
	int i,j;
    for (i=0;i<m_problem_size-1;i++) {VjsNonMean[i]=0;}
    
	//Calculate Vjmean vector from the population
	//cout << "Calculate Vjsmean vector"<<endl;
	Invert(consensus_ranking,m_problem_size,invertedB);
	int * individua;
	for (i=0;i<samples_num;i++)
	{
		individua=samples[i];
		
        Compose(individua,invertedB,composition,m_problem_size);
        
		vVector_Fast(Vjs, composition, m_problem_size,aux_v);
        
		for (j=0;j<m_problem_size-1;j++)
		{
			VjsNonMean[j]=VjsNonMean[j]+Vjs[j];
		}
	}
	
	for (i=0;i<m_problem_size-1;i++)
	{
		VjsMean[i]=(float)VjsNonMean[i]/samples_num;
	}
    
    m_newton_d->Vjsmean=VjsMean;
    for (j=1;j<m_problem_size;j++)
	{
        m_newton_d->j=j;
		theta_parameters[j-1]= Newton(0.0001, m_upper_theta_bound, &f_GeneralizedKendall, &fdev_GeneralizedKendall, ((void *)m_newton_d));
	}
    
    for (i=0;i<m_problem_size-1;i++){
        m_theta_parameters[i]=MAX(m_theta_parameters[i],m_lower_theta_bound);
        m_theta_parameters[i]=MIN(m_theta_parameters[i],m_upper_theta_bound);
    }
}

/*
 * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
 */
void GeneralizedKendall_Model::CalculatePsiConstants(float* thetas, float* psi)
{
	//en el param psi se dejan los valores de las ctes d normalización psi_j
	//returns psiTotal: la multiplicacion de los n-1 psi_j
	int n=m_problem_size;
	//calculate psi
	int i, j;
	for(i=0;i< n - 1;i++)
	{
		j=i+1;//el indice en el paper es desde =1 .. n-1
		psi[i] = (1 - exp((-1)*(n-j+1)*(thetas[i])))/(1 - exp((-1)*(thetas[i])));
	}
}

/*
 * Generates a permutation from a v vector.
 */
void GeneralizedKendall_Model::GeneratePermuFromV(int*v,int*permu)
{
 /*
     int i,j;
     permu[v[0]]=0;
     for (i=1;i<n_;i++){
     j=0;
     while (v[i]>=0)
     {
     v[i] = v[i] - (permu[j]==-1);
     j++;
     }
     permu[j-1]=i;
     }
    */

    int val, i, index;
    for(i = 0 ; i < m_problem_size ; i ++) aux_n[i]=i;
    for(i = 0 ; i < m_problem_size - 1 ; i++){
        val = v[i];
        index = 0;
        while( !(aux_n[ index ] != -1 && val == 0))
            if(aux_n[ index ++ ] != -1)
                val --;
        permu[ i ] = index  ;
        aux_n[ index ] = -1 ;
    }
    index=0;
    while(aux_n[ index ] == -1 )index++;
    permu[ m_problem_size - 1 ] = index ;

}
