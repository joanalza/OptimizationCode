/*
 *  PlackettLuce.cpp
 *  MixtureMallowsEDA
 *
 *  Created by Josu Ceberio Uribe on 7/16/12.
 *  Copyright 2012 University of the Basque Country. All rights reserved.
 *
 */

#include "PlackettLuce.h"
#include "Tools.h"
#include <math.h>
#include <complex>
using namespace std;
ofstream fweights;
/*
 * Class constructor.
 */
CPlackettLuceModel::CPlackettLuceModel(int problem_size, int sel_size)
{
	m_problem_size=problem_size;
    m_samples_size=sel_size;
    
    m_MM_iterations=10000;

	m_weights= new float[m_problem_size];
    m_mostProbable= new CIndividual(m_problem_size);

    
    int *aux1 = new int[m_problem_size];
    for (int i=0;i<m_problem_size;i++)
        aux1[i]=i;
    m_mostProbable->SetGenes(aux1);
    
    m_aux1= new int[m_problem_size];
    m_aux2= new int[m_problem_size];
    m_sampling_permutation= new int[m_problem_size];

    //creating sampling structures.
    distribution= new double[m_problem_size];
    objects= new int[m_problem_size];
    m_last= new int[m_problem_size];
    m_sample_size_huge= m_samples_size * 10; //1000000;
    /*m_huge_sample= new int*[m_sample_size_huge];
    for (int i=0;i<m_sample_size_huge;i++){
        m_huge_sample[i]=new int[m_problem_size];
    }*/

        //create learning structures.
        m_M = m_problem_size;

//        m_N = m_sample_size_huge; //number of rankings in the population.
        m_N = m_samples_size; //number of rankings in the population.
        m_P = m_problem_size;
    
        r2 = new double*[m_M];
        f= new int*[m_P];
        r= new int*[m_M];
        aux= new double*[m_P];
        g= new double*[m_P];
        for (int i=0;i<m_P;i++)
        {
            f[i]= new int[m_N];
            r[i]= new int[m_N];
            r2[i]= new double[m_N];
            aux[i]= new double[m_N];
            g[i]= new double[m_N];
        }
        
        dgamma= new double[m_M];
        gamma= new double[m_M];
        newgamma= new double[m_M];
    
        w= new double[m_M];
        pp= new int[m_N];
        sr2= new double[m_M];
   
    m_p1=new int[m_problem_size];
    m_p2= new int[m_problem_size];
    m_located= new int[m_problem_size];
}

/*
 * Class destructor.
 */
CPlackettLuceModel::~CPlackettLuceModel()
{
    delete [] m_sampling_permutation;
	delete [] m_weights;
    delete m_mostProbable;
    
    /*for (int i=0;i<m_sample_size_huge;i++){
        delete [] m_huge_sample[i];
    }
    delete [] m_huge_sample;*/
    
    delete [] m_last;
    delete [] distribution;
    delete [] objects;
    delete [] m_aux1;
    delete [] m_aux2;
    
        for (int i=0;i<m_P;i++)
        {
            delete [] f[i];
            delete [] r[i];
            delete [] r2[i];
            delete [] aux[i];
            delete [] g[i];
        }
        delete [] f;
        delete [] r;
        delete [] r2;
        delete [] aux;
        delete [] g;
    
        delete [] dgamma;
        delete [] gamma;
        delete [] newgamma;
        delete [] w;
        delete [] pp;
        delete [] sr2;
    
    delete [] m_p1;
    delete [] m_p2;
    delete [] m_located;
}

void cumsum(double * values, int size, double * resul)
{
	resul[0]=values[0];
	for (int i=1;i<size;i++)
	{
		resul[i]=resul[i-1]+values[i];
	}
}
void cumsum(double ** values,int size1, int size2, double ** resul)
{
	for (int j=0;j<size2;j++)
		resul[0][j]=values[0][j];
	
	for(int i=1;i<size1;i++)
		for (int j=0;j<size2;j++)
			resul[i][j]=values[i][j]+resul[i-1][j];
}

double sum (double * values, int size)
{
	double res=0;
	for (int i=0;i<size;i++) { res+=values[i];}
	return res;
}
double norm(double * values, int size)
{
	double result=0;
	for (int i=0;i<size;i++)
		result+=pow(values[i],2);
	return sqrt(result);
}

bool CPlackettLuceModel::Learn(CPopulation * population, int size)
{

    bool result=LearnMM(population,size);

    return result;
}

/*
 * Learns a Mallows model under the Ulam distance from a boltzmann distribution, and returns the full Mallows distribution
 */
double * CPlackettLuceModel::Learn_fromBoltzmann(int ** search_space, double * boltzmann_probabilities, int size){
    
    //1.- Calculate the consensus ranking/permutation.
    LearnMM_fromBoltzmann(search_space,boltzmann_probabilities,size);

    double * full_distribution= new double[size];
    for (int i=0;i<size;i++){
        full_distribution[i]=Probability(search_space[i]);
//        cout<<full_distribution[i]<<endl;
    }
    //CalculateMostProbableSolution();
    //PrintArray(m_mostProbable->Genes(),m_problem_size,"most probable ");
    return full_distribution;
}

/*
 * Samples individuals from the distribution
 */
void CPlackettLuceModel::Sample_BoltzmannDistribution(int ** search_space, double * distribution, int size, int sample_size){

    int* samples_indices = new int[size];

    for (int i=0;i<size;i++) samples_indices[i]=0;
    
    double random_value;
    int index;
    double acumul;
    for (int i=0;i<sample_size;i++){
        
        random_value=((double)rand() / (double)(RAND_MAX));
        index=0;
        acumul=0;
        while (acumul<random_value) {
            acumul+=distribution[index];
            index++;
        }
        samples_indices[index]++;
    }

    //copy permus to the huge sample.
    
    int sampled_permus=0;
    
    /*for (int i=0;i<size;i++){
        for (int j=0;j<samples_indices[i];j++){
            memcpy(m_huge_sample[sampled_permus],search_space[i],sizeof(int)*m_problem_size);
            sampled_permus++;
        }
    }*/
}

/*
 * It estimates the Plackett-Luce model from the given Boltzmann distribution over the search space with the MM algorithm.
 */
bool CPlackettLuceModel::LearnMM_fromBoltzmann(int ** search_space, double * boltzmann_probabilities, int size)
{

    m_N=m_sample_size_huge;
    Sample_BoltzmannDistribution(search_space, boltzmann_probabilities, size, m_N);
    //This code is a parse of the code of Hunter et al. 2004 of the MM algorithm implemented in matlab.
    

    //Initialization of auxiliary data.
    int i,j, aux_i,aux_j;
    //	int M = ProblemSize;
    //	int N = sel_total; //number of rankings in the population.
    //	int P = ProblemSize;
    
    int * individua;
    
    //extract individuals from population
    for (i = 0; i < m_N-1; i++)
    {
        /*individua = m_huge_sample[i];
        for (j = 0; j < m_M; j++)
        {
            //we consider that individuals in the population have an ordering description
            //and these matrix consider ordering description of the solutions.
            f[j][i]=individua[j];
            r[individua[j]][i]= j + m_P*i;
        }*/
    }
    
    
    GenerateRandomPermutation( m_last,m_problem_size);
    for (j = 0; j < m_M; j++)
    {
        //we consider that individuals in the population have an ordering description
        //and these matrix consider ordering description of the solutions.
        f[j][i]=m_last[j];
        r[m_last[j]][i]= j + m_P*i;
    }
    
    //create a copy of r in r2.
    for (i=0;i<m_M;i++) for(j=0;j<m_N;j++)r2[i][j]=r[i][j];
    
    
    for (i=0;i<m_M;i++) w[i]=0;
    for (i=0;i<m_N;i++) pp[i]=m_M;
    
    for (i=0;i<m_N;i++)
        for (j=0;j<m_P-1;j++)
            w[f[j][i]]++;
    
    
    for (i=0;i<m_N;i++) pp[i]+= i*m_P;
    
    //Starts the MM algorithm
    
    for (i=0;i<m_M;i++) {gamma[i]=1;dgamma[i]=1;}
    
    int iterations=0;
    //double ** aux= new double*[P];
    m_norm=1;
    
    //cout<<"Starting MM algorithm...."<<endl;
    while(m_norm>0.01 && iterations<m_MM_iterations)
    {
       // cout<<"norm: "<<m_norm<<"iterations: "<<iterations<<endl;
        //PRIMERA PARTE
        iterations++;
        
        // for (i=0;i<N;i++) for(j=0;j<N;j++)g[i][j]=f[i][j]; //esta linea en matlab hace falta, pero aqui no.
        
        for(i=0;i<m_P;i++)
            for(j=0;j<m_N;j++)
                g[i][j]=gamma[f[i][j]];
        
        //for (j=0;j<N;j++)
        //		aux[P-1][j]=g[P-1][j];
        
        memcpy(aux[m_P-1], g[m_P-1], sizeof(double)*m_N);
        
        for(i=m_P-2;i>=0;i--)
            for (j=0;j<m_N;j++)
                aux[i][j]=g[i][j]+aux[i+1][j];
        
        for(i=m_P-1;i>=0;i--)
            memcpy(g[i],aux[i],sizeof(double)*m_N);
        
        
        for(j=0;j<m_N;j++)
            g[m_P-1][j]=0;
        
        
        for(i=0;i<m_M-1;i++)
            for(j=0;j<m_N;j++)
                g[i][j]=(double)1/g[i][j]; //en este for anidado va a dar error porque la ultima fila de g es de zeros
        //como debe, y no estamos seguros d que se controle. Hay que ver si el N-1 esta bien puesto en el primer
        //for
        
        cumsum(g,m_P,m_N,g);
        
        //SEGUNDA PARTE
        for (j=0;j<m_N;j++)
        {
            for(i=0;i<m_M;i++)
            {
                aux_i= r[i][j] % m_P;
                aux_j= r[i][j] / m_P;
                r2[i][j]=g[aux_i][aux_j];
            }
        }
        
        for (i=0;i<m_M;i++)
        {
            sr2[i] = sum(r2[i],m_N);
            newgamma[i] = w[i] / sr2[i];
            dgamma[i] = newgamma[i] - gamma[i];
        }
        
        memcpy(gamma, newgamma, sizeof(double)*m_M);
        
        m_norm=norm(dgamma, m_M);
    }
    
    
    double sumGamma=sum(gamma,m_M);
    
    for (int i=0;i<m_M;i++) m_weights[i]=gamma[i]/sumGamma;
    
    return true;
}

/*
 * It estimates the Plackett-Luce model from the given cases. The learning process consists
 * of finding the MLE w parameters. For that task, Minorise-Maximize (MM) algorithm is utilized.
 */
bool CPlackettLuceModel::LearnMM(CPopulation * population, int sel_total)
{
	//This code is a parse of the code of Hunter et al. 2004 of the MM algorithm implemented in matlab.
	
	//Initialization of auxiliary data.
	int i,j, aux_i,aux_j;
//	int M = ProblemSize;
//	int N = sel_total; //number of rankings in the population.
//	int P = ProblemSize;

	int * individua;
	//extract individuals from population
	for (i = 0; i < m_N-1; i++)
	{
		individua = population->m_individuals[i]->Genes();
		for (j = 0; j < m_M; j++)
		{
            //we consider that individuals in the population have an ordering description
			//and these matrix consider ordering description of the solutions.
            f[j][i]=individua[j];
            r[individua[j]][i]= j + m_P*i;
		}
	}
    
    
    GenerateRandomPermutation( m_last,m_problem_size);
    for (j = 0; j < m_M; j++)
    {
        //we consider that individuals in the population have an ordering description
        //and these matrix consider ordering description of the solutions.
        f[j][i]=m_last[j];
        r[m_last[j]][i]= j + m_P*i;
    }
 
    //create a copy of r in r2.
	for (i=0;i<m_M;i++) for(j=0;j<m_N;j++)r2[i][j]=r[i][j];

    
	for (i=0;i<m_M;i++) w[i]=0;
	for (i=0;i<m_N;i++) pp[i]=m_M;

	for (i=0;i<m_N;i++)
		for (j=0;j<m_P-1;j++)
			w[f[j][i]]++;	
		
	
	for (i=0;i<m_N;i++) pp[i]+= i*m_P;
	
	//Starts the MM algorithm

    for (i=0;i<m_M;i++) {gamma[i]=1;dgamma[i]=1;}

	int iterations=0;
	//double ** aux= new double*[P];
	m_norm=1;

	//cout<<"Starting MM algorithm...."<<endl;
	while(m_norm>0.001 && iterations<m_MM_iterations)
	{
		//cout<<"norm: "<<m_norm<<"iterations: "<<iterations<<endl;
		//PRIMERA PARTE
		iterations++;

		// for (i=0;i<N;i++) for(j=0;j<N;j++)g[i][j]=f[i][j]; //esta linea en matlab hace falta, pero aqui no.
 
		for(i=0;i<m_P;i++)
			for(j=0;j<m_N;j++)
				g[i][j]=gamma[f[i][j]];
		
		//for (j=0;j<N;j++)
	//		aux[P-1][j]=g[P-1][j];
        
        memcpy(aux[m_P-1], g[m_P-1], sizeof(double)*m_N);
        
		for(i=m_P-2;i>=0;i--)
			for (j=0;j<m_N;j++)
				aux[i][j]=g[i][j]+aux[i+1][j];

		for(i=m_P-1;i>=0;i--)
			memcpy(g[i],aux[i],sizeof(double)*m_N);
	
		
		for(j=0;j<m_N;j++)
			g[m_P-1][j]=0;
	
		
		for(i=0;i<m_M-1;i++)
			for(j=0;j<m_N;j++)
				g[i][j]=(double)1/g[i][j]; //en este for anidado va a dar error porque la ultima fila de g es de zeros
		//como debe, y no estamos seguros d que se controle. Hay que ver si el N-1 esta bien puesto en el primer
		//for
		
		cumsum(g,m_P,m_N,g);

		//SEGUNDA PARTE
		for (j=0;j<m_N;j++)
		{
			for(i=0;i<m_M;i++)
			{
				aux_i= r[i][j] % m_P;
				aux_j= r[i][j] / m_P;
				r2[i][j]=g[aux_i][aux_j];
			}
		}
		
		for (i=0;i<m_M;i++)
		{
			sr2[i] = sum(r2[i],m_N);
			newgamma[i] = w[i] / sr2[i];
			dgamma[i] = newgamma[i] - gamma[i];
		}

		memcpy(gamma, newgamma, sizeof(double)*m_M);

        m_norm=norm(dgamma, m_M);
	}

    
	double sumGamma=sum(gamma,m_M);

    for (int i=0;i<m_M;i++) m_weights[i]=gamma[i]/sumGamma;
    //PrintArray(m_weights, m_M, "weights. ");
    
    CalculateMostProbableSolution();
    
    return true;
}

/*
 * Samples the plackett-luce model learnt.
 */
inline int CPlackettLuceModel::Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)
{
	int i,j,k , remaining,pos, samples;
	double randVal,acumul,acumul_aux;
    
    for (samples=0;samples<num_samples ;samples++){
        
        for (i=0;i<m_problem_size;i++) { m_aux1[i]=-1;m_located[i]=0;}
        for (i=0;i<m_problem_size;i++)
        {

            //get ready data
            remaining=m_problem_size-i;
            acumul_aux=0;
            pos=0;
            for (j=0;j<m_problem_size;j++)
            {
                if (m_located[j]==0)//permutation[j]==-1)
                {
                    distribution[pos]=m_weights[j];
                    acumul_aux+=m_weights[j];
                    objects[pos]=j;
                    pos++;
                }
            }

            randVal=rand()*acumul_aux;
            //randVal=(double)rand()/(((double)RAND_MAX/acumul_aux));
            acumul=distribution[0];
            for(k=1;(k<remaining && acumul<randVal);k++)
            {
                acumul += distribution[k];
            }
            k--;
        
            m_sampling_permutation[i]=objects[k];
            m_located[objects[k]]=1;
        }
       // if (inverse)
       //     population->AddToPopulation(m_sampling_permutation, samples, problem->EvaluateInv(m_sampling_permutation));
       // else
            population->AddToPopulation(m_sampling_permutation, samples, problem->Evaluate(m_sampling_permutation));
    }
    return num_samples;
}

#/*
 * Samples the plackett-luce model learnt.
 */
void CPlackettLuceModel::Sample(int *permutation)
{
	int i, j, k, remaining, pos, samples;
	double randVal, acumul, acumul_aux;

	for (i = 0; i<m_problem_size; i++) m_located[i] = 0;
	for (i = 0; i<m_problem_size; i++)
	{

		//get ready data
		remaining = m_problem_size - i;
		acumul_aux = 0;
		pos = 0;
		for (j = 0; j<m_problem_size; j++)
		{
			if (m_located[j] == 0)//permutation[j]==-1)
			{
				distribution[pos] = m_weights[j];
				acumul_aux += m_weights[j];
				objects[pos] = j;
				pos++;
			}
		}

		//randVal = rand()*acumul_aux;
		randVal=(double)rand()/(((double)RAND_MAX/acumul_aux));
		acumul = distribution[0];
		for (k = 1; (k<remaining && acumul<randVal); k++)
		{
			acumul += distribution[k];
		}
		k--;

		permutation[i] = objects[k];
		m_located[objects[k]] = 1;
	}

}

/*
 * This method updates the best solution of the optimization processes if the most probable solutions
 * are better than those solutions obtained by the population.
 */
void CPlackettLuceModel::PostProcesses()
{
 /*   if (BEST->Value()<m_mostProbable->Value())
    {
        BEST->SetGenes(m_mostProbable->Genes());
        BEST->SetValue(m_mostProbable->Value());
        CONVERGENCE_ITERATION=GEN_NUM;
        CONVERGENCE_EVALUATIONS=EVALUATIONS;
        
    }*/
}

/*
 * Returns the probability of the 'indiv' individual for Plackett-Luce model.
 */
float CPlackettLuceModel::Probability(int * indiv)
{
    double probability=1;
    double aux,aux2;
    int i,j;

        for (i=0;i<m_problem_size;i++)
        {
            aux=m_weights[indiv[i]];
            aux2=0;
            for (j=i;j<m_problem_size;j++)
                aux2=aux2+m_weights[indiv[j]];
            probability=probability*(aux/aux2);
        }
    
	return probability;
}

void CPlackettLuceModel::GibbsWeigths(float t) {
	float * abilities = new float[m_problem_size];
	float sum = 0;
	for (int i = 0; i < m_problem_size; i++) {
		abilities[i] = log(m_weights[i]);
		sum += exp(abilities[i] / t);
	}
	for (int i = 0; i < m_problem_size; i++) {
		m_weights[i] = (exp(abilities[i] / t))/(sum);
	}
	delete[] abilities;
}


/*
 * Returns the most probable solution given the vector of weights
 */
void CPlackettLuceModel::CalculateMostProbableSolution()
{
    RandomKeys(m_aux1,m_weights,m_problem_size);
    
    //Invert
    for(int i=0; i<m_problem_size; i++) m_mostProbable->Genes()[m_aux1[i]]=i;
    
    m_mostProbable->SetValue(INT_MIN);
}

void CPlackettLuceModel::printWeigths() {
	double sum = 0;
	for (int i = 0; i < m_problem_size; i++) {
		cout << "w" << i << ":" << m_weights[i] << endl;
		sum += m_weights[i];
	}
	cout << "Sum: " << sum << endl;
}

/*
 * Returns the most probable solution given the vector of weights
 */
CIndividual * CPlackettLuceModel::GetMostProbableSolution()
{
    return m_mostProbable;
}



