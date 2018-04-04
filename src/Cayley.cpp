/*
 *  Cayley.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/18/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "Cayley.h"
#include "Tools.h"
#include "Variables.h"
#include "NewtonRaphson.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using std::cerr;
using std::cout;
using std::endl;


/*
 * The constructor of the model.
 */
Cayley_Model::Cayley_Model(int problem_size, int sel_size)
{
	m_problem_size=problem_size;
    m_sample_size=sel_size;
	m_psis= new float[m_problem_size-1];
	m_consensus_ranking=new int[m_problem_size];
    m_theta_parameter=1;
    
	//reservar memoria para las probabilidades de Xs.
	m_xprobs=new float[m_problem_size-1];
    
    m_upper_theta_bound=10;
    m_lower_theta_bound=0.001;
    m_newton_d = (cayley_data *) malloc(sizeof(cayley_data));
    m_newton_d->n=m_problem_size;
    
    m_composed = new int[ m_problem_size ];
    m_inverted = new int[ m_problem_size ];
    m_x= new int [m_problem_size];
    m_visited = new bool[ m_problem_size ];
    m_aux= new int [m_problem_size];
    

}

/*
 * The destructor of the model.
 */
Cayley_Model::~Cayley_Model()
{
	delete [] m_xprobs;
	delete [] m_psis;
	delete [] m_consensus_ranking;
    delete [] m_composed;
    delete [] m_inverted;
    delete [] m_x;
    delete [] m_aux;
    delete [] m_visited;
    delete m_newton_d;
}

/*
 * Learns a Mallows model based on the Cayley distance and the individuals in the model.
 */
bool Cayley_Model::Learn(int ** samples, int size)
{

	//1. calculate global consensus ranking
    CalculateConsensusRanking(samples, size, m_consensus_ranking);

    //2. calculate theta parameters
	m_theta_parameter= CalculateThetaParameter(m_consensus_ranking, samples, size);
    
    //3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameter, m_psis);
		
	//4.- Calculate Xj probabilities vector.
	for(int j=0;j<m_problem_size-1;j++)
        m_xprobs[j]=(float)1/m_psis[ j ];/*(double) exp(-m_theta_parameter)/m_psis[ j ];*/

	return true;
}

/*
 * Learns a Mallows model based under the Kendall distance for the given infinite individual.
 */
bool Cayley_Model::Learn(int ** samples, float * probabilities, int size){
    
    //1.- Calculate the consensus ranking/permutation.
    CalculateConsensusRanking_WeightedSamples(samples, size, probabilities, m_consensus_ranking);
    
    //2.- Calculate theta parameter
	m_theta_parameter= CalculateThetaParameter_WeightedSamples(m_consensus_ranking, samples, size,probabilities);
    
    cout<<"Theta: "<<m_theta_parameter<<endl;
	//3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameter, m_psis);
	
    return true;
}


/*
 * Builds the Mallows model for the Cayley distance with the given CR and theta parameters.
 */
bool Cayley_Model::Learn(int * consensus_ranking, double theta)
{
    //1.- Copy the consensus ranking
    memcpy(m_consensus_ranking, consensus_ranking, sizeof(int)*m_problem_size);
    
    //2.- Copy the theta parameter.
	m_theta_parameter= theta;
    
    //3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameter, m_psis);
    
	//4.- Calculate Xj probabilities vector.
	for(int j=0;j<m_problem_size-1;j++)
        m_xprobs[j]=(float)1/m_psis[ j ];/*(double) exp(-m_theta_parameter)/m_psis[ j ];*/
    
	return true;
}

/*
 * Given the consensus ranking, it samples a new individual.
 */
void Cayley_Model::Sample(int * permutation)
{
    for(int i= 0;i < m_problem_size - 1; i ++ ){
//        probaX_j = (float) 1/m_psis[ i ];
        if(((float)rand() / RAND_MAX) < m_xprobs[i])
            m_x[ i ] = 0;
        else
            m_x[ i ] =1;
    }
    m_x[m_problem_size-1] = 0;
    
    GeneratePermuFromX(m_x, m_aux);
    Compose(m_aux, m_consensus_ranking, permutation, m_problem_size);
}

/*
 * Calculates the sum of Cayley distances of the given solution to the sample of solutions.
 */
float Cayley_Model::CalculateDistancetoSample(int * solution, int ** samples, int cases_num){

    float dist=0;
    for (int i=0;i<cases_num;i++){
        Invert(samples[i],m_problem_size,m_inverted);
        Compose(solution, m_inverted, m_composed, m_problem_size);
        dist+=CalculateDistanceAndX(m_composed, NULL);
    }
    return dist;
}

/*
 * Calculates the sum of Cayley distances of the given solution to the sample of solutions.
 */
float Cayley_Model::CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, float * weights){
    
    float dist=0;
    for (int i=0;i<cases_num;i++){
        Invert(samples[i],m_problem_size,m_inverted);
        Compose(solution, m_inverted, m_composed, m_problem_size);
        dist+=CalculateDistanceAndX(m_composed, NULL)*weights[i];
    }
    return dist;
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void Cayley_Model::CalculateConsensusRanking(int ** samples, int cases_num, int*consensusPermutation)
{
    int option=2;
    if (option==0){
        //Option 1. Choose the best solution
        for (int i=0;i<m_problem_size;i++)
            consensusPermutation[i]=samples[0][i];
    }
    else if (option==1){
        //Option 2. Calculate the MLE consensus. Ekhine's method.
        double fitness=0;
       Approximated_MLE_Consensus(cases_num, samples, consensusPermutation, &fitness);
    }
    else{
        //Option 3. Choose the solution that minimizes the distance to the sample.
        int best_index=0;
        int best_distance= INT_MAX;
        int dist;
        for (int i=0;i<cases_num;i++){
            dist=CalculateDistancetoSample(samples[i], samples, cases_num);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
        memcpy(consensusPermutation, samples[best_index],sizeof(int)*m_problem_size);
    }

}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void Cayley_Model::CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, float * weights, int* consensus_ranking)
{
    //Option 1. Choose the solution that minimizes the distance to the sample.
    int best_index=0;
    float best_distance= INT_MAX;
    float dist;
    for (int i=0;i<sample_size;i++){
        dist=CalculateDistancetoSample_WeightedSamples(samples[i], samples, sample_size, weights);
        
        if (dist<best_distance)
        {
            best_distance=dist;
            best_index=i;
        }
    }
    memcpy(consensus_ranking, samples[best_index],sizeof(int)*m_problem_size);
    cout<<"CR index: "<<best_index<<",  prob: "<<weights[best_index]<< "       ";    PrintArray(consensus_ranking,m_problem_size,"");
    
}


/*
 * Theta parameter estimation function.
 */
float f_Cayley(float theta, void * d) {
    cayley_data * data = (cayley_data*)d;
    double sum=0, ex,denom;
    for(int j = 1 ; j < data->n ; j++){
        ex = exp(theta);
        denom = j+ex;
        sum += (double)j/denom;
        //sum += (double) j /(double)(exp(theta)+(double)j);
    }
    return (double)(sum - data->dist_avg);
}

/*
 * Theta parameter estimation function derivation.
 */
float fdev_Cayley(float theta, void * d) {
        cayley_data * data = (cayley_data*)d;
    double sum=0;
    for(int j = 1 ; j < data->n ; j++)
        sum += (double)( - j * exp( theta ))/pow(exp(theta) + j, 2);
    return sum;
}

/*
 * Calculates the spread theta parameter from the m_consensus_ranking and the individuals in the population.
 */
float Cayley_Model::CalculateThetaParameter(int*consensus, int ** samples, int cases_num)
{

        //Option 2. Calculate the theta parameter according to the spread of the solutions around the consensus ranking.
        int s, distance= 0;
        Invert(consensus, m_problem_size, m_inverted);
        for(s = 0 ; s < cases_num ; s ++){
            Compose(samples[s], m_inverted, m_composed, m_problem_size);
            distance += CalculateDistanceAndX(m_composed, NULL);
        }

        m_newton_d->dist_avg=(double)distance/(double)cases_num;
    
    m_theta_parameter=Newton(0.0001,m_upper_theta_bound, &f_Cayley, &fdev_Cayley, ((void *)m_newton_d));
    m_theta_parameter=MAX(m_theta_parameter,m_lower_theta_bound);
    m_theta_parameter=MIN(m_theta_parameter,m_upper_theta_bound);
    
    return m_theta_parameter;
    
}

/*
 * Calculates the spread theta parameter from the m_consensus_ranking and the individuals in the population.
 */
float Cayley_Model::CalculateThetaParameter_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, float * weights)
{
    
    //Option 2. Calculate the theta parameter according to the spread of the solutions around the consensus ranking.
    int s, distance= 0;
    Invert(consensus_ranking, m_problem_size, m_inverted);
    for(s = 0 ; s < samples_num ; s ++){
        Compose(samples[s], m_inverted, m_composed, m_problem_size);
        distance += CalculateDistanceAndX(m_composed, NULL)*weights[s];
    }
    
    m_newton_d->dist_avg=distance;
    
    m_theta_parameter=Newton(0.0001,100, &f_Cayley, &fdev_Cayley, ((void *)m_newton_d));

    return m_theta_parameter;
    
}

/*
 * Calculates the total Psi normalization constant from the m_theta_parameter and psi-s vector.
 */
float Cayley_Model::CalculatePsiConstants(float theta, float* psi)
{
    int i,j;
    for (i = 0 ; i <  m_problem_size - 1; i ++ ){
        j=i+1;
        psi[ i ] = 1 + ( m_problem_size - j) *(double) exp(-theta);
    }
	return 0;
}

/*
 * Calculates the probability of the given individual in the learnt Mallows probabilistic model with the Cayley distance.
 */
double Cayley_Model::Probability(int * individual){
    
    Invert(m_consensus_ranking, m_problem_size, m_inverted);
    Compose(individual, m_inverted, m_composed, m_problem_size);
    CalculateDistanceAndX(m_composed, m_x);
    
    double proba=1;
    
    for(int i = 0 ; i < m_problem_size - 1 ; i ++ )
        proba *= exp(-m_theta_parameter*m_x[ i ])/m_psis[ i ];

    return proba;
 
}

/*
 * Generates a permutation from a x vector.
 */
void Cayley_Model::GeneratePermuFromX(int *x, int *sigma){

    int item ;//from 0..n-1 !!!!
    for (int i = 0 ; i < m_problem_size ; i++) m_inverted [ i ] = -1;
    for (int pos = 0 ; pos < m_problem_size ; pos ++) {
        if (x[ pos ] == 0){// an item <= pos can close the cycle
            item = pos ;
            while ( m_inverted [ item ] != -1 || ! item_closes_cycle(pos, item, sigma, m_inverted)) item --;
        }else{
            item = - 1 ;
            int random = rand() % (m_problem_size - pos - 1);
            while ( random  >= 0 ){
                item++;
                if (m_inverted [ item  ] == -1 &&  ! item_closes_cycle(pos , item, sigma, m_inverted)  )     random --;
            }
        }
        m_inverted[ item ] = pos;
        sigma[ pos ] = item;
    }
}

inline bool Cayley_Model::item_closes_cycle(int pos, int item, int *sigma, int *sigma_inv){
    //    item -= FIRST_ITEM;
    while (item < pos ) item = sigma[ item ];
    if ( item == pos ) return true;
    return false;
}

int Cayley_Model::CalculateDistanceAndX(int * sigma, int *x){
    //also updates the x vector if it isnot null
    
    if(x!=NULL)for (int i = 0 ; i < m_problem_size; i ++ )x[ i ] =1;
    
    int num_cycles=0, num_visited=0, item= 0;

    for (int i = 0 ; i < m_problem_size; i ++ )
        m_visited[ i ] =false;
    while(num_visited<m_problem_size){
        item=num_cycles;
        while(m_visited[item])item++;
        num_cycles++;
        int maxItemInCycle= 0;
        do{
            if(item>maxItemInCycle)maxItemInCycle=item;
            m_visited[item] =true;
            num_visited++;
            item=sigma[item];
        }while(!m_visited[item]);
        if(x!=NULL)x[maxItemInCycle] = 0;
    }
    return (m_problem_size-num_cycles);
}


/*
 * Obtains the most probable cycle when building the permutation from a given vector of X.
 */
int Cayley_Model::GetMostProbCycle(int ind, int **cycles, int len, int *leng_cycles){
    //returns the cycle index
    int i,j;
    int cont = 0;
    while(ind>cycles[cont][0]) cont++;
    //each cycles from cycle[cont] to cycle[len-1] has a particuler
    //probability of being selected as the terget
    //this prob is
	float* proba_cycles = new float[len];
	float* proba_cyclesNoAcum = new float[len];

    /*float proba_cycles[len];
    float proba_cyclesNoAcum[len];*/

    for (i = 0 ; i < len; i ++ ) proba_cycles[ i ] = 0;
    float acumul= 0;
    float multi;
    for(i=cont ; i < len; i ++ ){
        for(j=cont;j<len;j++){
            //if(j==i)proba_cycles[ i ] *= (float)gen.factorial(lengCycles[ j ]);
            //else proba_cycles[ i ] *= (float) gen.factorial(lengCycles[ j ]-1);

            if(j==i) multi= (float)factorial(leng_cycles[ j ]);
            else  multi= (float) factorial(leng_cycles[ j ]-1);
            
            if(proba_cycles[ i ] ==0) proba_cycles[ i ] =multi;
            else proba_cycles[ i ] *= multi;
        }
        acumul += proba_cycles[ i ];
    }
    
    for (i = 0 ; i < len; i ++ ) proba_cyclesNoAcum[ i ] =(float)proba_cycles[ i ]/acumul;
    for(i=cont ; i < len; i ++ )
        if(i!=0)proba_cycles[ i ] = (float)proba_cycles[ i ]/acumul + (float)proba_cycles[i-1];
        else proba_cycles[ i ] = (float)proba_cycles[ i ]/acumul ;
    float ran = (float)rand()/(RAND_MAX);

    for(i=cont ; i < len; i ++ )
        if(proba_cycles[ i ]>ran){
            //*partialProb=proba_cyclesNoAcum[ i ];
            return i;
        }
    return len-1;
}

/*
 * Generates a random permutaiton with K cycles.
 */
int * Cayley_Model::GeneratePermutation_withKcycles(int n, int k){//generateUARPermuWithKCyles
    int *sigma = new int[n];
    if(n==k) for (int i = 0 ; i < n; i ++ )sigma[ i ] =i+1;
    else
        if(k==1){
            int *cycle = new int[n];
            GenerateRandomPermutation(cycle, n);

            for (int i = 0 ; i < n-1; i ++ ) sigma[cycle[ i ]] = cycle[i+1]+1;
            sigma[cycle[n-1]] =cycle[0]+1;
            delete []cycle;
        }else{
            long double  ran1 =((long double)rand() / ((double )RAND_MAX +1 ) );
            int ran2=-1;
            if(ran1 < (long double)(stirling_matrix_[ n - 1 ][ k - 1 ] / stirling_matrix_[ n ][ k ])){//stirling1(n-1, k-1) / stirling1(n,k))){////n is in a cycle by itself
                /* long double ran1 = rand() % (long double)stirling1(n,k),  ran2=-1; if(ran1 < stirling1(n-1, k-1)){ */
                int *small=GeneratePermutation_withKcycles(n-1, k-1);
                for (int i = 0 ; i < n-1; i ++ )sigma[ i ] =small[ i ];
                sigma[n-1] =n;
                delete []small;
            }else{
                int *small=GeneratePermutation_withKcycles(n-1, k);
                for (int i = 0 ; i < n-1; i ++ )sigma[ i ] =small[ i ];
                ran2 = rand() % (n-1) ;//0..n-2
                sigma[n-1] =sigma[ran2];
                sigma[ran2] =n;
                delete []small;
            }
        }
    return sigma;
}


void Cayley_Model::Approximated_MLE_Consensus( int m, int **samples, int *sigma_0, double *best_likeli){
    int **  samples_inv = new int*[ m ];
    int ** samples_copy = new int*[ m ];
    for (int i = 0 ; i < m; i++) {
        samples_copy[ i ] = new int[ m_problem_size ];
        samples_inv [ i ] = new int[ m_problem_size ];
        for(int j = 0 ; j < m_problem_size ; j ++) {
            samples_inv[ i ][ samples[ i ][ j ]] = j;
            samples_copy[ i ][ j ] = samples[ i ][ j ];
        }
    }
    approx_consensus_mm(m, samples_copy, samples_inv, sigma_0, best_likeli);
    
    //variable_neighborhood_search(m, samples, sigma_0 , best_likeli);
    for (int i = 0 ; i < m; i++) {
        delete [] samples_inv[ i ];
        delete [] samples_copy[ i ];
    }
    delete [] samples_copy;
    delete [] samples_inv;
}

void Cayley_Model::approx_consensus_mm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_distance){
    // pre: sigma, sigmaInv =-1
    int distance_increase = 0, remaining = m_problem_size;
    int **freq = new int*[ m_problem_size ];for (int i = 0; i < m_problem_size; i++) freq[ i ] = new int[ m_problem_size ];
    
    for (int i = 0; i < m_problem_size; i++) sigma_0[ i ] = -1;
    do{ //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  1"<<endl;}
        for (int i = 0; i < m_problem_size; i++)
            for (int j = 0; j < m_problem_size; j++)freq[ i ][ j ] = 0;
        
        int max_freq = 0;
        int pi = -1, pj = -1; // 0..n-1
        bool dirty_items = true ; //si hay alguno q aparezca >m/2 veces?
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  2"<<endl;}
        for (int s = 0; dirty_items && s < m; s++)
            for (int i = 0; ( dirty_items && i < m_problem_size ) ; i++){
                if (sigma_0[ i ] == -1) (freq[ i ][samples_copy[ s ][ i ] ])++;
                if ((freq[ i ][samples_copy[ s ][ i ] ]) > max_freq){
                    max_freq = (freq[ i ][samples_copy[ s ][ i ] ]);
                    pi = i;
                    pj = samples_copy[ s ][ i ] ;
                    if (max_freq > m/2 ) dirty_items = false ;
                }
            }
        //for (int i = 0 ; i < n; i ++ ){for (int j = 0 ; j < n;j++)cout<<freq[ i ][ j ]<<" ";cout<<endl;}
        sigma_0[pi] = pj ;
        for (int s = 0; s < m; s++) {
            if (samples_copy[ s ][pi] != pj ) {//////////////swap
                //swapSamplesDo(s, pi, pj, samples, samples_inv, posSwap);
                int x = samples_copy[ s ][pi];
                int y = samples_inv[ s ][pj] ;
                samples_copy[ s ][pi] = pj ;
                samples_copy[ s ][y] = x;
                samples_inv[ s ][pj] = pi ;
                samples_inv[ s ][x ] = y ;
                distance_increase++;
            }
        }
        remaining--;
    }while(remaining > 0) ;
    (*best_distance) = distance_increase;
    for (int i = 0 ; i < m_problem_size ; i ++) delete [] freq[ i ];
    delete [] freq;
}

void Cayley_Model::variable_neighborhood_search(int m, int **samples, int *sigma, double *f_eval){
  //  bool improve;
  //  do{
        double f_eval_ini = (*f_eval);
     //   cout<<" distance ini"<<(*f_eval)<<endl;
        local_search_swap_mm  (m, samples, sigma, f_eval);
        
     //   cout<<" distance new swap "<<(*f_eval)<<endl;
        //local_search_insert(m, samples, sigma, f_eval);
        //cout<<" distance new ins "<<(*f_eval)<<endl;
    //    improve=false ;
  //      if( (f_eval_ini) > *f_eval ) improve = true;
    //}while(improve);
}

int Cayley_Model::get_cycles(int *sigma, int *cycle_items, int *cycle_indices){
    bool*visited = new bool[ m_problem_size ];
    for (int i = 0 ; i < m_problem_size; i ++ )visited[ i ] =false;
    int item_index, cont=0, cycle_index= 0;
    while(cont<m_problem_size) {
        item_index= 0;
        while(visited[item_index])item_index++;
        while (!visited[item_index]) {
            visited[item_index] = true;
            cycle_items[cont] =item_index;
            cycle_indices[cont] =cycle_index;
            item_index=sigma[item_index];
            cont++;
        }
        cycle_index++;
    }
    delete [] visited;
    return cycle_index;
}

void Cayley_Model::local_search_swap_mm(int m, int **samples, int *sigma_0, double *f_eval){
    int ** samples_comp = new int*[ m ];
    for (int s = 0 ; s < m; s++)  samples_comp[ s ] = new int[ m_problem_size ];
    int *sigma_0_inv = new int[ m_problem_size ];
    int *cycle_items = new int [m_problem_size ], *cycle_index = new int[ m_problem_size];
    int **same_cycle = new int*[ m_problem_size ];
    for (int i = 0 ; i < m_problem_size ; i ++){ same_cycle[ i ] = new int[ m_problem_size ]; for (int j = 0 ; j < m_problem_size ; j++) same_cycle[ i][ j ] = 0;}

    int index_i,index_j, distance_variation = 0;
    for(int i = 0 ; i < m_problem_size ; i ++) sigma_0_inv[ sigma_0 [ i ] ] = i ;
    bool improve;
    do {//iterating use just sigma_0_inv
        int max_freq = 0;
        for (int s = 0 ; s < m; s++) {
            for(int j = 0 ; j < m_problem_size ; j ++ ) samples_comp[ s ][ j ] = samples[ s ][ sigma_0_inv [ j ]  ];

            //cout<<"Cycles ";gen.print_int_vector(cycle_items, n_); gen.print_int_vector(cycle_index, n_);
            for ( int i = 0 ; i < m_problem_size ; i ++){
                for (int j = i + 1; j < m_problem_size && (cycle_index[ i ] == cycle_index[ j ]) ; j++){
                    int max, min; // for a triangular matrix
                    if (cycle_index [ i ] > cycle_index[ j ]){ max = cycle_items[ i ] ; min  = cycle_items[ j ];}
                    else {min  = cycle_items[ i ]; max  = cycle_items[ j ];}
                    same_cycle[ min ][ max ] ++;
                    if (max_freq < same_cycle[ min ][ max ] ) {
                        max_freq = same_cycle[ min ][ max ];
                        index_i = min;
                        index_j = max;
                    }
                }
            }
        }
        //there are m oermus, max_freq of them are going to decrease dist in 1. (m-max_freq) are going to incresse
        distance_variation = m - 2 * max_freq;
        
        improve = false;
        if( distance_variation < 0 ) {
            improve = true;
            int aux = sigma_0_inv[index_i];
            sigma_0_inv[index_i] = sigma_0_inv[index_j];
            sigma_0_inv[index_j] = aux;
            
            (*f_eval) += distance_variation;
            for(int i = 0 ; i < m_problem_size ; i ++) sigma_0[ sigma_0_inv[ i ]  ] = i ;
            
        }
        
        //cout<<"trace in local search swap, distance variation "<<distance_variation<<" mxa freq "<<max_freq<<endl;
    } while (improve );
    
    delete [] sigma_0_inv;
    delete [] cycle_index;
    delete [] cycle_items;
    for (int i = 0 ; i < m_problem_size ; i++) delete [] same_cycle[ i ];
    delete [] same_cycle;
    for (int i = 0 ; i < m ; i++)  delete [] samples_comp[ i ];
    delete [] samples_comp;
}

int Cayley_Model::distance_to_sample(int m, int **samples, int *sigma){
    int distance= 0;
    int *comp = new int[ m_problem_size ], *sigma_inv = new int[ m_problem_size ];
    for(int j = 0 ; j < m_problem_size ; j ++) sigma_inv[sigma[ j ] ] = j ;
    for(int s = 0 ; s < m ; s ++){
        for(int i = 0 ; i < m_problem_size ; i ++) comp[ i ] = samples[ s ][ sigma_inv [ i ]  ];
        distance += CalculateDistanceAndX(comp, NULL);
    }
    delete []sigma_inv;
    delete []comp;
    return distance ;
}

/*
 * Determines if the given indiviual is a local optima for the problem and the current distance.
 */
bool Cayley_Model::isLocalOptima(int * individual, PBP * problem){
    
    //asumiendo que la distancia de Cayley es el numero de swaps, vamos a comprar la solucion del individuo con todas aquellas
    //que estan a un swap de distancia.
    
    double cost_neigh;
    int j,k,aux;
    double cost_individual=problem->Evaluate(individual);
    int * solution = new int [m_problem_size];
    memcpy(solution, individual, sizeof(int)*m_problem_size);
    for ( k=0; k<m_problem_size-1; k++) {
        for (j=k+1; j<m_problem_size; j++) {
            aux=solution[k];
            solution[k]=solution[j];
            solution[j]=aux;
            cost_neigh=problem->Evaluate(solution);
            if (cost_neigh>cost_individual) {
                delete[] solution;
                return false;
            }
            aux= solution[j];
            solution[j]=solution[k];
            solution[k]=aux;
        }
    }
    delete [] solution;
    return true;
}

