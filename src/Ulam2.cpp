//
//  Ulam2.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 03/02/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#include "Ulam2.h"
#include "Tools.h"
#include "NewtonRaphson.h"
#include <vector>
#include <set>
#include <list>
#include "Variables.h"

/*
 * The constructor.
 */
Ulam_Model2::Ulam_Model2(int problem_size, int sel_size)
{
    
    m_problem_size=problem_size;
    m_sample_size=sel_size;
	//m_psis= new float[m_problem_size-1];
	m_consensus_ranking=new int[m_problem_size];
    m_theta_parameter=1;
    m_lower_theta_bound=0.001;
    m_upper_theta_bound=10;
    m_probabilities= new double[m_problem_size];
    m_composed = new int[ m_problem_size ];
    m_inverted = new int[ m_problem_size ];
    m_averageDistance=0;
 
    memcpy(m_str_base_path,"./",sizeof(char)*3);

 
    if (ReadPermusPerDistance(m_problem_size)==false){
        //for (int i=10;i<60;i=i+10){
        //    m_problem_size=i;
            SaveCounts_toFile();
        //}
        exit(1);
        ReadPermusPerDistance(m_problem_size);
    }
    
    Generic gen;
    
    m_col_index=new int[m_problem_size];
    m_row_index=new int[m_problem_size];
    m_shape = new int[ m_problem_size ];
    m_newton_d = (ulam_data2 *) malloc(sizeof(ulam_data2));
    m_newton_d->n=m_problem_size;

}
/*
 * The destructor.
 */
Ulam_Model2::~Ulam_Model2()
{
    delete [] m_probabilities;
	delete [] m_consensus_ranking;
    //delete [] m_psis;
    delete [] m_composed;
    delete [] m_inverted;
    
    delete [] m_num_perms_distance;
    
    
    delete [] m_col_index;
    delete [] m_row_index;
    delete m_newton_d;
}

/*
 * Reads from file the number of permutations to each Ulam distance for n size of permutations.
 */
bool Ulam_Model2::ReadPermusPerDistance(int n){
    if ( m_num_perms_distance == NULL ){
        m_num_perms_distance = new long double[ n ];
        char integer_string[5];
        sprintf(integer_string, "%d", n);
        char str_permus_per_dist[500];
        strcpy(str_permus_per_dist, m_str_base_path);
        strcat(str_permus_per_dist, "permus_per_dist_");
        strcat(str_permus_per_dist, integer_string);  // other_string now contains "Integer: 1234"
        ifstream  file;
        file.open(str_permus_per_dist);
        if(!((file ))){
            cout << "Cannot read input file - permus_per_dist "<<endl;
            //exit(1);
            return false;
        }
        for (int i = 0 ; i < n ; i ++)
            file >> m_num_perms_distance[ i ];
        file.close();
    }
    return true;
}

/*
 * Reads...
 */
void Ulam_Model2::ReadPermusPerShape(int d, int bound, int * shape, int*shape_len){
    string line;
    long double acum;
    bool end = false;
    char integer_string[5];
    
    sprintf(integer_string, "%d", m_problem_size);
    char str_permus_per_shape[500] ;
    strcpy (str_permus_per_shape, m_str_base_path);
    strcat(str_permus_per_shape, "permus_per_shape_");
    strcat(str_permus_per_shape, integer_string);
    strcat(str_permus_per_shape, "_");
    sprintf(integer_string, "%d", d);
    strcat(str_permus_per_shape, integer_string);
    
    ifstream  file;
    file.open(str_permus_per_shape);
    if(!((file ))){
        cout << "Cannot read input file - permus_per_dist "<<endl;
        exit(1);
    }
    do{
        file >> acum ;
        //cout<<acum<<endl;
        if ( acum < bound )
            file.ignore(m_problem_size * 2 , '\n');
        else {
            getline(file, line);
            end = true;
        }
    }while ( ! end );
    file.close();
    
    stringstream tokenizer( line );
    int lambda_i;
    while( tokenizer >> lambda_i ){
        //cout << "Token: " << lambda_i << endl;
        shape[ (*shape_len)++ ] =lambda_i;
    }
}

/*
 * Calculates the number of permutations at each distance Ulam, and saves into files.
 */
void Ulam_Model2::SaveCounts_toFile(){
    /* adapted from ZS1: Zoghbi, Antoine and Stojmenovic, Ivan: Fast Algorithms for Generating Integer Partitions. International Journal of Computer Mathematics, 70, 1998, 319-332.
     generaes partitions in  anti-lexicographic order */
    ofstream file_permus_per_dist;
    ofstream file_shapes_per_dist;
    ofstream file_permus_per_shape;
    
    char integer_string[5];
    sprintf(integer_string, "%d", m_problem_size);
    char str_permus_per_shape[500] ;//= "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_shape_"; // un file e estos por cada n y distance
    strcpy(str_permus_per_shape, m_str_base_path);
    strcat(str_permus_per_shape, "permus_per_shape_");
    
    char str_permus_per_shape_n_d[600];
    char str_permus_per_dist[600] ;//= "/Users/eki/Dropbox/permus/prj/perms_mallows/permus_per_dist_";
    strcpy(str_permus_per_dist, m_str_base_path);
    strcat(str_permus_per_dist, "permus_per_dist_");
    
    strcat(str_permus_per_shape, integer_string);
    strcat(str_permus_per_dist, integer_string);
    
    unsigned char k;                  //length of figures
    unsigned char *vector     = NULL; //where the current figure is stored
    int           gen_result;         //return value of generation functions
    int           prev_distance = -1 , dist ;
    int           part_len;
    long double     permus_per_shape = 0;
    vector = (unsigned char *)malloc(sizeof(unsigned char) * m_problem_size );
    if(vector == NULL)    {
        fprintf(stderr, "error: insufficient memory\n");
        exit(EXIT_FAILURE);
    }
    m_num_perms_distance = new long double[ m_problem_size ];
    for (int i = 0 ; i < m_problem_size ;i++) m_num_perms_distance[ i ] = 0 ;
    
    Ferrers_diagram2 *f;
    gen_result = gen_part_init(vector, m_problem_size, &k);
    while(gen_result == GEN_NEXT ) {
        part_len = (int)k;
        int*part = new int[part_len];//DO NOT delete, member of Ferrers_diagram
        for(int i = 0 ; i < part_len; i++) part[i]=(int)vector[i];
        dist = part[ 0 ];
        
        f = new Ferrers_diagram2(m_problem_size, part , part_len);
        f->calculate_hook_length_greene_niejenhuis_wilf();
        dist = f->get_resulting_distance();
        m_num_perms_distance[ dist ] += f->get_num_permus();
        if ( dist != prev_distance ){
            permus_per_shape = 0;
            cout<<"Generating shape at distance "<<dist<<endl;
            if (file_permus_per_shape.is_open() )
                file_permus_per_shape.close();
            strcpy(str_permus_per_shape_n_d, str_permus_per_shape);
            strcat(str_permus_per_shape_n_d, "_");
            sprintf(integer_string, "%d", dist );
            strcat(str_permus_per_shape_n_d, integer_string);
            file_permus_per_shape.open(str_permus_per_shape_n_d);
        }
        permus_per_shape += f->get_num_permus();
        file_permus_per_shape << permus_per_shape <<" ";
        for ( int i =  0; i < part_len ; i++ )
            file_permus_per_shape<<part[ i ]<<" ";
        file_permus_per_shape<<endl;
        prev_distance = dist;
        gen_result = gen_part_next(vector, &k, 0);
        delete f ;
    }
    free(vector);
    
    file_permus_per_shape.close();
    file_permus_per_dist.open (str_permus_per_dist);
    for (int i = 0 ; i < m_problem_size; i++ ){
        file_permus_per_dist<< m_num_perms_distance[ i ]<<endl ;
    }
    delete [] m_num_perms_distance;
    file_permus_per_dist.close();
}


/*
 * Learns a Mallows model based on the Ulam distance and the individuals in the model.
 */
bool Ulam_Model2::Learn(int ** samples, int size)
{
	//1. calculate global consensus ranking
    int min_dist=CalculateConsensusRanking(samples, size, m_consensus_ranking);
    
    //2. calculate theta parameters
	m_theta_parameter= CalculateThetaParameter(m_consensus_ranking, samples, size, min_dist);
    
    //3. Calculate probabilities for sampling.
    m_probabilities[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < m_problem_size ; i++)//acumulate the number of permus at each distance
        m_probabilities[i] = m_num_perms_distance[i] * exp ( -m_theta_parameter * i ) + m_probabilities[ i - 1 ];
    
  //      double prob=1/m_probabilities[m_problem_size-1];
  //  cout<<"prob: "<<prob<<endl;exit(1);
	return true;
}

/*
 * Learns a Mallows model based under the Ulam distance for the given infinite individual.
 */
bool Ulam_Model2::Learn(int ** samples, float * probabilities, int size){
    
    //1.- Calculate the consensus ranking/permutation.
    float min_dist=CalculateConsensusRanking_WeightedSamples(samples, size, probabilities, m_consensus_ranking);
    
    //2.- Calculate theta parameter
	m_theta_parameter= CalculateThetaParameter_WeightedSamples(m_consensus_ranking, samples, size,min_dist);
    
    cout<<"Theta: "<<m_theta_parameter<<endl;
    
    //3. Calculate probabilities for sampling.
    m_probabilities[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < m_problem_size ; i++)//acumulate the number of permus at each distance
        m_probabilities[i] = m_num_perms_distance[i] * exp ( -m_theta_parameter * i ) + m_probabilities[ i - 1 ];

	
    return true;
}

/*
 * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
 */
bool Ulam_Model2::Learn(int * consensus_ranking, double theta)
{
    //1.- Copy the consensus ranking
    memcpy(m_consensus_ranking, consensus_ranking, sizeof(int)*m_problem_size);
    
    //2.- Copy the theta parameter.
	m_theta_parameter= theta;
    
    //3. Calculate probabilities for sampling.
    m_probabilities[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < m_problem_size ; i++)//acumulate the number of permus at each distance
        m_probabilities[i] = m_num_perms_distance[i] * exp ( -m_theta_parameter * i ) + m_probabilities[ i - 1 ];
    
    //      double prob=1/m_probabilities[m_problem_size-1];
    //  cout<<"prob: "<<prob<<endl;exit(1);
	return true;
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
int Ulam_Model2::CalculateConsensusRanking(int **samples, int cases_num, int * consensusPermutation){
    int min_dist = m_problem_size * cases_num;
    int best_index = -1;
    int dist=0;
    for(int i = 0; i < cases_num ; i++){
        dist = CalculateDistancetoSample(samples[i], samples, cases_num);
        if (dist < min_dist) {
            min_dist = dist;
            best_index = i;
        }
    }
    memcpy(consensusPermutation, samples[best_index], sizeof(int)*m_problem_size);
    return min_dist;
}



/*
 * Calculates de consensus permutation from the given population cases.
 */
float Ulam_Model2::CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, float * weights, int* consensus_ranking){
    float min_dist = m_problem_size * sample_size;
    int best_index = -1;
    int dist=0;
    for(int i = 0; i < sample_size ; i++){
        dist= CalculateDistancetoSample_WeightedSamples(samples[i], samples, sample_size, weights);
        if (dist < min_dist) {
            min_dist = dist;
            best_index = i;
        }
    }
    memcpy(consensus_ranking, samples[best_index], sizeof(int)*m_problem_size);
        cout<<"CR index: "<<best_index<<",  prob: "<<weights[best_index]<< "       ";    PrintArray(consensus_ranking,m_problem_size,"");
    return min_dist;
}


/*
 * Calculates the probability of the given individual in the learnt Mallows probabilistic model with the Ulam distance.
 */
double Ulam_Model2::Probability(int * individual){
    int dist=Ulam(individual, m_consensus_ranking, m_problem_size);
     return  exp(-dist * m_theta_parameter )/m_probabilities[ m_problem_size - 1 ];
}

/*
 * Theta parameter estimation function.
 */
float f_Ulam2(float theta, void * d){
    ulam_data2 * data = (ulam_data2*)d;
    double numer = 0, denom = 0;
    double aux;
    for (int d = 0 ; d < data->n - 1; d++){
        aux = data->count[d ] * exp(-theta *d ) ;
        numer += aux * d;
        denom += aux;
    }
    return numer / denom - data->dist_avg;
}

/*
 * Theta parameter estimation function derivation.
 */
float fdev_Ulam2(float theta, void * d){
    ulam_data2 * data = (ulam_data2*)d;
    double numer1 = 0, numer2 = 0, numer3 = 0, denom = 0;
    double aux;
    for (int d = 0 ; d < data->n - 1; d ++){
        aux = data->count[d] * exp(-theta *d ) ;
        numer1 += aux * d * d ;
        numer2 += aux;
        numer3 += aux * d;
        denom += aux;
    }
    return (-numer1 * numer2 - numer3*numer3)/(denom*denom);
}

/*
 * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
 */
float Ulam_Model2::CalculateThetaParameter(int*consensus, int ** samples, int cases_num, int dist_avg){
    
    m_newton_d->count=m_num_perms_distance;
    m_newton_d->dist_avg=(double)dist_avg/cases_num;
    m_theta_parameter=Newton(0.0001,m_upper_theta_bound, &f_Ulam2, &fdev_Ulam2, (void *)m_newton_d);

    m_theta_parameter=MAX(m_theta_parameter,m_lower_theta_bound);
    m_theta_parameter=MIN(m_theta_parameter,m_upper_theta_bound);
    return m_theta_parameter;
}

/*
 * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
 */
float Ulam_Model2::CalculateThetaParameter_WeightedSamples(int*consensus, int ** samples, int cases_num, float average_dist){
    
    m_newton_d->count=m_num_perms_distance;
    m_newton_d->dist_avg=average_dist;
    m_theta_parameter=Newton(0.0001,20, &f_Ulam2, &fdev_Ulam2, (void *)m_newton_d);
    
    return m_theta_parameter;
}

/*
 * Given the consensus ranking, it samples a new individual.
 */
void Ulam_Model2::Sample(int * permutation)
{
    double rand_distance = ((double)rand() / (double)(RAND_MAX)) * m_probabilities[m_problem_size-1];
    int target_distance = 0;
    while(m_probabilities[ target_distance ] < rand_distance)//ekhine pone <=
        target_distance++;

    GeneratePermuAtLIS( m_problem_size - target_distance, m_composed);
    
   // for (int i=0;i<m_problem_size;i++) m_composed[i]--;
    Compose(m_composed,m_consensus_ranking,permutation,m_problem_size);

}

/*
 * Given a LIS size, generates a permu at that distance.
 */
void Ulam_Model2::GeneratePermuAtLIS(int l, int *sigma){
    
    //initializations
    int     d = m_problem_size - l;
    int     to_insert;
    int     col, row, aux, new_col, new_row;
    

    double target_permus = (double)rand() / (double)(RAND_MAX) * m_num_perms_distance[ d ];
    int shape_len = 0;
    ReadPermusPerShape(d , target_permus, m_shape, &shape_len);
    
    int* shape1 = new int [ shape_len ];//member of f1
    int* shape2 = new int [ shape_len ];
    memcpy(shape1,m_shape,sizeof(int)*shape_len);
    memcpy(shape2,m_shape,sizeof(int)*shape_len);
    
    Ferrers_diagram2 * f1 = new Ferrers_diagram2(m_problem_size, shape1 , shape_len);
    Ferrers_diagram2 * f2 = new Ferrers_diagram2(m_problem_size, shape2 , shape_len);
    
    f1->random_SYT_greene_niejenhuis_wilf();
    f2->random_SYT_greene_niejenhuis_wilf();
    int ** tableau1 = f1->get_syt();
    int ** tableau2 = f2->get_syt();
    
    for (int i = 0 ; i < f2->get_ferrers_shape_length() ; i++){
        for (int j =  0 ; j < f2->get_ferrers_shape()[i] ; j++) {
            m_row_index[ tableau2[ i ][ j ] - 1 ] = i;
            m_col_index[ tableau2[ i ][ j ] - 1 ] = j;
        }
    }
    
    for (int index = m_problem_size - 1 ; index >= 0 ; index --){
        col = m_col_index[ index ];
        row = m_row_index[ index ];
        to_insert = tableau1[ row ][ col ];
        while (row != 0) {
            new_col=0, new_row = row - 1;
            while (f1->get_ferrers_shape()[new_row] > new_col+1
                   && tableau1[new_row][new_col + 1 ] < to_insert)
                new_col++;
            aux = tableau1[new_row][new_col];
            tableau1[new_row][new_col] = to_insert;
            to_insert = aux;
            row = new_row;
            col = new_col;
        }
        sigma[index ] = to_insert-1;
        tableau1[m_row_index[ index ]][ m_col_index[ index ]] = m_problem_size + 1;
        //gen.print_int_vector(sigma, n_);
    }

    delete f1;
    delete f2;
}

/*
 * Calculates the sum of Ulam distances of the given solution to the sample of solutions.
 */
int Ulam_Model2::CalculateDistancetoSample(int * solution, int ** samples, int cases_num){
    int dist = 0;
    for (int s= 0; s < cases_num; s ++) {
        Invert(samples[s],m_problem_size,m_inverted);
        Compose(solution, m_inverted, m_composed, m_problem_size);
        //for(int i = 0 ; i < m_problem_size ; i++) m_composed[ i ] = samples[ s ] [ sigma_inv [ i ]];
        dist += (m_problem_size- getLISLength(m_composed,m_problem_size));
    }
    return dist;
}

/*
 * Calculates the sum of Ulam distances of the given solution to the sample of solutions.
 */
float Ulam_Model2::CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, float * weights){
    float dist = 0;
    for (int s= 0; s < cases_num; s ++) {
        Invert(samples[s],m_problem_size,m_inverted);
        Compose(solution, m_inverted, m_composed, m_problem_size);
        //for(int i = 0 ; i < m_problem_size ; i++) m_composed[ i ] = samples[ s ] [ sigma_inv [ i ]];
        dist += (m_problem_size- getLISLength(m_composed,m_problem_size))*weights[s];
    }
    return dist;
}


int Ulam_Model2::gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k){
    int j; //index
    //test for special cases
    if(n == 0)    {
        (*k) = 0;
        return(GEN_EMPTY);
    }
    //initialize: vector[0] = n, vector[1, ..., n - 1] are 1
    vector[0] = n;
    for(j = 1; j < n; j++)        vector[j] = 1;
    (*k) = 1;
    return(GEN_NEXT);
}

int Ulam_Model2::gen_part_next(unsigned char *vector, unsigned char *k, int bound){
    //bound == 0 => no bound
    static int j = 0; //remember index of the rightmost part which is greater than 1
    int        r;     //temporary remainder
    int        temp;  //auxiliary element
    
    //easy case
    if(vector[j] == 2){
        vector[j] = 1;
        j--;
        (*k)++;
        //terminate if the num of columns is smaller than 'bound'
        if ((int)vector[0] < bound){
            j=0;
            return (GEN_TERM);
        }
        return(GEN_NEXT);
    }
    //terminate if all parts are 1
    if(vector[0] == 1){
        j = 0;
        return(GEN_TERM);
    }
    
    //decrease
    vector[j]--;
    temp = vector[j];
    r = *k - j;
    //set right-hand elements
    while(r > temp){
        j++;
        vector[j] = temp;
        r -= temp;
    }
    *k = j + 2;
    //set rightmost element
    if(r > 1){
        j++;
        vector[j] = r;
    }
    //terminate if the num of columns is smaller than 'bound'
    if ((int)vector[0] < bound){
        j=0;
        return (GEN_TERM);
    }
    return(GEN_NEXT);
}

/*
 * Determines if the given indiviual is a local optima for the problem and the current distance.
 */
bool Ulam_Model2::isLocalOptima(int * individual, PBP * problem){
    
    //asumiendo que la distancia de Ulam es el numero de inserts, vamos a comprar la solucion del individuo con todas aquellas
    //que estan a un insert de distancia.
    
       double cost_neigh;    int i,j;
    int * solution = new int [m_problem_size];
    int * aux = new int [m_problem_size];
    memcpy(solution, individual, sizeof(int)*m_problem_size);
    double cost_individual=problem->Evaluate(individual);

    
        //forward step
        for (i=0;i<(m_problem_size-1);i++)
        {
            memcpy(aux,solution, sizeof(int)*m_problem_size);
            Swap(aux,i,i+1);
            cost_neigh=problem->Evaluate(aux);
            if (cost_neigh>cost_individual)
            {
                delete [] aux;
                delete [] solution;
                return false;
            }
            
            for (j=i+1; j<(m_problem_size-1);j++)
            {
                Swap(aux,j,j+1);
                cost_neigh=problem->Evaluate(aux);
                if (cost_neigh>cost_individual)
                {
                    delete [] aux;
                    delete [] solution;
                    return false;
                }
            }
        }
        
        //backward step
        Swap(solution,0,1);
        for (i=2;i<m_problem_size;i++)
        {
            Swap(solution,0,i);
            memcpy(aux, solution, sizeof(int)*m_problem_size);
            cost_neigh=problem->Evaluate(aux);
            if (cost_neigh>cost_individual)
            {
                delete [] aux;
                delete [] solution;
                return false;
            }
            for (j=0;j<(i-2);j++)
            {
                Swap(aux,j,j+1);
                cost_neigh=problem->Evaluate(aux);
                if (cost_neigh>cost_individual)
                {
                    delete [] aux;
                    delete [] solution;
                    return false;
                }
            }
        }
    return true;
}

