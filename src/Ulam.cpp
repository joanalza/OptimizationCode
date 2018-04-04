//
//  Ulam.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/25/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "Ulam.h"
#include "Tools.h"
#include "NewtonRaphson.h"
#include <vector>
#include <set>
#include <list>
#include "Variables.h"

/*
 * The constructor.
 */
Ulam_Model::Ulam_Model(int problem_size, int sel_size)
{
 
    m_problem_size=problem_size;
    m_sample_size=sel_size;
	m_psis= new float[m_problem_size-1];
	m_consensus_ranking=new int[m_problem_size];
    m_theta_parameter=1;
    
	//reservar memoria para las probabilidades de Xs.
    
    m_upper_theta_bound=10;
    m_lower_theta_bound=1;
    m_probabilities= new double[m_problem_size];
    m_composed = new int[ m_problem_size ];
    m_inverted = new int[ m_problem_size ];
    m_averageDistance=0;
    
    
    //initialize bratteli diagrams.
    m_bound_lis=m_problem_size-m_problem_size/3;
    m_count_permus = new long double[m_problem_size];
    
    m_bratteli= new std::vector<Ferrers_diagram*>[ m_problem_size + 1];
    Generic gen;
    int i,l,j;
    fill_bratteli_diagram(m_bound_lis);
    
    for(i = 0 ; i < m_problem_size ; i++)
        m_count_permus[i] = 0;
    
    for (i = 0; i < m_bratteli[m_problem_size].size(); i++) {
        l = m_bratteli[m_problem_size].at(i)->get_num_cols();
        m_count_permus[m_problem_size-l] +=  m_bratteli[m_problem_size].at(i)->get_num_permus();
    }
    
    //PrintArray(m_count_permus, m_problem_size, "count permutations at dist ");
    
    long double acumul_permus = 0, remaining;
    for (i = 0 ;(m_count_permus[i] >0 &&  i < m_problem_size); i++)
        acumul_permus +=m_count_permus[i];
    
    if (acumul_permus != m_problem_size ) {
        //the remaining permus are split. The last pos = 1 (because S_u(n,n-1)=1)
        remaining = (gen.factorial(m_problem_size) - acumul_permus - 1 )/(m_problem_size-i-1);
        for(j = i ; j< m_problem_size - 1 ; j ++)
            m_count_permus[j]=remaining;
        m_count_permus[m_problem_size - 1 ] = 1;
        
        // PrintArray(m_count_permus, m_problem_size, "count permutations at dist ");
    }
    m_num_partitions=m_bratteli[m_problem_size].size();
    m_acumul_ld = new long double[m_num_partitions];
    m_acumul= new double[m_num_partitions];
    m_col_index=new int[m_problem_size];
    m_row_index=new int[m_problem_size];
    
    m_newton_d = (ulam_data *) malloc(sizeof(ulam_data));
    m_newton_d->n=m_problem_size;
}
/*
 * The destructor.
 */
Ulam_Model::~Ulam_Model()
{
    delete [] m_probabilities;
	delete [] m_consensus_ranking;
    delete [] m_composed;
    delete [] m_inverted;
    for(int i = 0; i < m_problem_size + 1 ; i ++)
        for(int j=0; j < m_bratteli[i].size() ; j++)
            delete m_bratteli[i].at(j);
    delete [] m_bratteli;
    
    delete [] m_count_permus;
    
    
    delete [] m_col_index;
    delete [] m_row_index;
    delete [] m_acumul_ld;
    delete [] m_acumul;
    delete m_newton_d;
}


/*
 * Learns a Mallows model based on the Ulam distance and the individuals in the model.
 */
bool Ulam_Model::Learn(int ** samples, int size)
{
	//1. calculate global consensus ranking
    int min_dist=CalculateConsensusRanking(samples, size, m_consensus_ranking);

    //2. calculate theta parameters
	m_theta_parameter= CalculateThetaParameter(m_consensus_ranking, samples, size, m_upper_theta_bound, min_dist);

    
    m_probabilities[ 0 ] = 1; // exp(-theta*d) = exp (-theta *0)
    for (int i = 1 ; i < (m_problem_size-m_bound_lis) ; i++)//acumulate the number of permus at each distance
        m_probabilities[i] = m_count_permus[i] * exp ( -m_theta_parameter * i ) + m_probabilities[ i - 1 ];
    for (int i = (m_problem_size-m_bound_lis) ; i < (m_problem_size) ; i++)//acumulate the number of permus at each distance
        m_probabilities[i] =m_probabilities[ i - 1 ];

    
	return true;
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
int Ulam_Model::CalculateConsensusRanking(int **samples, int cases_num, int * consensusPermutation){
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
 * Calculates the probability of the given individual in the learnt Mallows probabilistic model with the Ulam distance.
 */
double Ulam_Model::Probability(int * individual){
    
    /* int *x = new int[ m_problem_size ];
     
     Invert(m_consensus_ranking, m_problem_size, m_inverted);
     Compose(individual, m_inverted, m_composed, m_problem_size);
     CalculateDistanceAndX(m_composed, x);
     
     double proba=1;
     
     for(int i = 0 ; i < m_problem_size - 1 ; i ++ ) proba *= exp(-m_theta_parameter*x[ i ])/m_psis[ i ];
     
     delete []x;
     return proba;*/
    return 0;
}

/*
 * Theta parameter estimation function.
 */
float f_Ulam(float theta, void * d){
    ulam_data * data = (ulam_data*)d;
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
float fdev_Ulam(float theta, void * d){
    ulam_data * data = (ulam_data*)d;
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
float Ulam_Model::CalculateThetaParameter(int*consensus, int ** samples, int cases_num, float upper_bound, int dist_avg){
    
    m_newton_d->count=m_count_permus;
    m_newton_d->dist_avg=(double)dist_avg/cases_num;
    return Newton(5.001,m_upper_theta_bound, &f_Ulam, &fdev_Ulam, (void *)m_newton_d);
}


/*
 * Given the consensus ranking, it samples a new individual.
 */
void Ulam_Model::Sample(int * permutation)
{
    double rand_distance = ((double)rand() / (double)(RAND_MAX)) * m_probabilities[m_problem_size-1];
    int target_distance = 0;
    //cout<<"rand_ "<<rand_distance<<endl;
    while(m_probabilities[ target_distance ] < rand_distance)//ekhine pone <=
        target_distance++;
    
    GeneratePermuAtLIS( m_problem_size - target_distance, m_composed);
    Compose(m_composed,m_consensus_ranking,permutation,m_problem_size);
}

/*
 * Prints the Bratteli diagram.
 */
void Ulam_Model::print_bratteli_diagram(){
    for(int i= m_problem_size  ;i > m_problem_size-1; i--){
        cout<<"----- ni = "<<i<<" -------"<<endl;
        for(int j=0;j<m_bratteli[i].size();j++){
            m_bratteli[i].at(j)->toStr();
        }
    }
}

/*
 * Given a LIS size, generates a permu at that distance.
 */
void Ulam_Model::GeneratePermuAtLIS(int l, int *sigma){
    
    //initializations
    Generic gen;
    Ferrers_diagram*shape;
    
    long double nperms;
    
    int i;
    
    //randomly choose a shape for the Young Standard Tableaux (YST)
    // cout<<"num_partitions: "<<num_partitions<<" lis:  "<<l<<endl;
    // print_bratteli_diagram();
    for(i=0;i<m_num_partitions;i++){
        // cout<<"hola: cols="<<m_bratteli[m_problem_size].at( i )->get_num_cols() <<endl;
        nperms=0;
        if(m_bratteli[m_problem_size].at( i )->get_num_cols() == l){
            nperms=m_bratteli[m_problem_size].at( i  ) -> get_num_permus();
        }
        m_acumul_ld [i] = nperms;
        if(i != 0)
            m_acumul_ld [i] += m_acumul_ld [i-1];
    }
    //PrintArray(acumul_ld, num_partitions, "acumul_ld: ");
    if(m_acumul_ld[m_num_partitions - 1 ] == 0){
        cout<<"Permu generation not possible. Refill Bratteli diagram to current target distance. "<<endl;
        exit(1);
    }
    
    for (i = 0 ; i < m_num_partitions ; i++) m_acumul[i]=m_acumul_ld[i] / m_acumul_ld[m_num_partitions-1];
    double bound =( (double)rand()/RAND_MAX);
    int pos=0;
    while(m_acumul[pos] < bound) pos++; //ekhine pone <=
    shape = m_bratteli[m_problem_size].at(pos);
    
    //create 2 YST with the ferrer diagram guiven by "shape, len"
    int**tableau1=new int*[shape->get_ferrers_shape_length()];
    int**tableau2=new int*[shape->get_ferrers_shape_length()];
    for(i = 0 ; i < shape->get_ferrers_shape_length(); i ++){
        tableau1[i]=new int[shape->get_num_cols()];
        tableau2[i]=new int[shape->get_num_cols()];
    }//todo ini 0
    
    GenerateRandomYST(shape, tableau1, NULL, NULL);
    GenerateRandomYST(shape, tableau2 , m_row_index , m_col_index);
    
    //generate permutation
    int index, col,row,to_insert,aux, new_col, new_row;
    for (index = m_problem_size - 1 ; index >= 0 ; index --){
        col = m_col_index[index];
        row = m_row_index[index];
        to_insert = tableau1[row][col];
        while (row != 0) {
            new_col=0;
            new_row = row - 1;
            while (shape->get_ferrers_shape()[new_row] > new_col+1 && tableau1[new_row][new_col + 1 ] < to_insert)
                new_col++;
            aux = tableau1[new_row][new_col];
            tableau1[new_row][new_col] = to_insert;
            to_insert = aux;
            row = new_row;
            col = new_col;
        }
        sigma[index ] = to_insert-1;//nosotros de 0 a n-1
        tableau1[m_row_index[index ]][m_col_index[index ]] = m_problem_size + 1;
    }
    
    for(int i = 0 ; i < shape->get_ferrers_shape_length(); i ++){
        delete [] tableau1[i];
        delete [] tableau2[i];
    }
    delete [] tableau1;
    delete [] tableau2;
}

/*
 * Generates a random Young Standard Tableaux given a ferrer diagram.
 */
void Ulam_Model::GenerateRandomYST(Ferrers_diagram *dad, int **out, int*row_index, int*col_index){
    Ferrers_diagram * son;
    int num_children;
    int i,j, target, pos;
    long double bound;
    long double * num_yst_child_acum;
    for( pos = m_problem_size ; pos > 1 ; pos-- ){
        
        //son <- randomly select one descendant
        num_children = dad->get_descendants()->size();
        num_yst_child_acum =new long double[num_children];
        num_yst_child_acum[0] = dad->get_descendants()->at(0)->get_hook_len();
        for(i =  1 ; i < num_children ; i ++)
            num_yst_child_acum[i] = num_yst_child_acum[i-1] + dad->get_descendants()->at(i)->get_hook_len();
        bound = (long double)rand() / (long double)(RAND_MAX) * num_yst_child_acum[num_children - 1];
        target=0;
        while(num_yst_child_acum[target] < bound) target++;//ekhine pone <=
        
        //if (target==num_children) target--;
        
        son = dad->get_descendants()->at(target);
        
        //look for the positions that differ in the son and father
        i=0, j=0;
        if ( dad->get_ferrers_shape_length() != son ->get_ferrers_shape_length()) {
            j = 0;
            i = son->get_ferrers_shape_length();
        }else{
            while(dad->get_ferrers_shape()[i] == son->get_ferrers_shape()[i])
                i ++;
            j = son->get_ferrers_shape()[i];
        }
        out[i][j] = pos;
        if(col_index != NULL){
            row_index[ pos - 1 ] = i;
            col_index[ pos - 1 ] = j;
        }
        delete [] num_yst_child_acum;
        dad=son;
    }
    out[0][0]=1;
    if(col_index != NULL){
        row_index[0] = 0;
        col_index[0] = 0;
    }
}

/*
 * Calculates the sum of Ulam distances of the given solution to the sample of solutions.
 */
int Ulam_Model::CalculateDistancetoSample(int * solution, int ** samples, int cases_num){
    int dist = 0;
    for (int s= 0; s < cases_num; s ++) {
        Invert(samples[s],m_problem_size,m_inverted);
        Compose(solution, m_inverted, m_composed, m_problem_size);
        //for(int i = 0 ; i < m_problem_size ; i++) m_composed[ i ] = samples[ s ] [ sigma_inv [ i ]];
        dist += (m_problem_size- getLISLength(m_composed,m_problem_size));
    }
    return dist;
}



void Ulam_Model::fill_bratteli_diagram(int bound_lis){
    //terminate if the num of columns is smaller than 'bound'
    //   if(m_bratteli == NULL){
    int i;
    for(i= m_problem_size ; i > 0 ; i--) //insert tableax for which the col number is <= n and >= bound
        insert_shapes_bound(i, (i != m_problem_size ), m_problem_size,  bound_lis-m_problem_size+i );
    for( i = 0 ; i < m_bratteli[m_problem_size].size() ; i ++)
        m_bratteli[m_problem_size].at(i)->calculate_hook_length();
    //print_m_bratellidiagram();
    /*    }else {
     //TODO insert til new bound
     if ((m_bratteli[m_problem_size].at(m_bratteli[m_problem_size].size()-1)->get_num_cols() > bound_lis && bound_lis != 0 )||
     (m_bratteli[m_problem_size].at(m_bratteli[m_problem_size].size()-1)->get_num_cols() != 1 && bound_lis == 0)) {
     //cout<<"Bratelli is full up to "<<m_bratelli[m_problem_size].at(m_bratelli[m_problem_size].size()-1)->get_num_cols()<<" columns. Run again to match new bound."<<endl;
     for(int i= m_problem_size ; i > 0 ; i--){
     int last_inserted_cols = m_bratteli[ i ].at(m_bratteli[ i ].size()-1)->get_num_cols();
     insert_shapes_bound(i, (i != m_problem_size ), last_inserted_cols - 1 ,  bound_lis-m_problem_size+i );
     }
     for( int i = 0  ; i < m_bratteli[m_problem_size].size() ; i ++) m_bratteli[m_problem_size].at(i)->calculate_hook_length();
     //print_m_bratellidiagram();
     
     }
     }
     */
}

void Ulam_Model::insert_shapes_bound(int current_n, bool has_parents, int bound_from, int bound_to){
    /*ZS1: Zoghbi, Antoine and Stojmenovic, Ivan: Fast Algorithms for Generating Integer Partitions. International Journal of Computer Mathematics, 70, 1998, 319-332.
     generaes partitions in  anti-lexicographic order */
    //terminate if the num of columns of the FS is smaller than 'bound'
    unsigned char k;                  //length of figures
    unsigned char *vector     = NULL; //where the current figure is stored
    int           gen_result;         //return value of generation functions
    // int           x;                  //iterator
    int           i, part_len;
    
    //alloc memory for vector
    vector = (unsigned char *)malloc(sizeof(unsigned char) * current_n);
    if(vector == NULL)    {
        fprintf(stderr, "error: insufficient memory\n");
        exit(EXIT_FAILURE);
    }
    gen_result = gen_part_init(vector, current_n, &k);
    while(gen_result == GEN_NEXT ) {//&& (int)vector[0] >= bound
        //for(x = 0; x < k; x++)            printf("___%u ", vector[x]); cout<<endl;
        if (vector[ 0 ] <= bound_from) {
            part_len = (int)k;
            int*part = new int[part_len];//DO NOT delete, member of Ferrers_diagram
            for(i = 0 ; i < part_len; i++)
                part[i]=(int)vector[i];
            Ferrers_diagram*f=new Ferrers_diagram(current_n, part , part_len);
            if( has_parents )
                for(i=0;i<m_bratteli[current_n + 1].size();i++){//TODO: since they are ordered thereis no need of searching over the whole array
                    if(m_bratteli[current_n + 1].at(i)->is_parent_of(f))
                        m_bratteli[current_n + 1 ].at(i)->relate_with(f);
                }
            if( bound_from != m_problem_size){//if it was included afterwards it must be related with the sons on the tree
                for(i = 0 ; i < m_bratteli[current_n - 1].size() ; i++){//TODO: since they are ordered thereis no need of searching over the whole array
                    if( f->is_parent_of(m_bratteli[current_n - 1].at(i)))
                        f->relate_with(m_bratteli[current_n - 1].at(i));
                }
            }
            m_bratteli[current_n].push_back(f);
            //cout<<"generated "<<endl;f->toStr();
        }
        gen_result = gen_part_next(vector, &k, bound_to);
    }
    free(vector);
    
}

int Ulam_Model::gen_part_init(unsigned char *vector, const unsigned char n, unsigned char *k){
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

int Ulam_Model::gen_part_next(unsigned char *vector, unsigned char *k, int bound){
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


