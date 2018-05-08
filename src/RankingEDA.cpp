/*
 *  RankingEDA.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 9/15/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "RankingEDA.h"
#include "MallowsModel.h"
#include "Cayley.h"
#include "Kendall.h"
#include "PlackettLuce.h"
#include "GeneralizedMallowsModel.h"

/*
 * The constructor.
 */
RankingEDA::RankingEDA(PBP * problem, int problem_size, long int max_evaluations, char* model_type, char* metric_type, int inverse)
{
    //1. standard initializations
    m_problem = problem;
	m_problem_size = problem_size;
    m_max_evaluations = max_evaluations;
    m_evaluations = 0;
    m_convergence_evaluations = 0;
    m_best = new CIndividual(problem_size);
	//long temp1 = MIN_LONG_INTEGER;
	//cout << "temp0" << temp1 << endl;
	m_best->SetValue(LONG_MIN);
	//cout << "initial m_best value: " << m_best->Value() << endl;
    m_pop_size = m_problem_size * 10;
    m_sel_size = m_problem_size;
    m_offspring_size = m_pop_size;
    memcpy(m_metric_type, metric_type, sizeof(char)*10);
    memcpy(m_model_type, model_type, sizeof(char)*10);
    
    m_inverse = inverse;
    
    //2. initialize the population
    m_population = new CPopulation(m_pop_size, m_offspring_size, m_problem_size);
    int *genes = new int[m_problem_size];
    for(int i=0; i<m_pop_size; i++){
		//Create random individual
		GenerateRandomPermutation(genes,m_problem_size);
        if (m_inverse)
            m_population->SetToPopulation(genes, i, m_problem->EvaluateInv(genes));
        else
            m_population->SetToPopulation(genes, i, m_problem->Evaluate(genes));
        m_evaluations++;
    }
    
    m_population->SortPopulation(0);

	m_best->SetGenes(m_population->m_individuals[0]->Genes());
	m_best->SetValue(m_population->m_individuals[0]->Value());

    //cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<endl;
    delete [] genes;
    
    
    //3. Build model structures
    if (((string)model_type)=="M"){
        m_model=new CMallowsModel(m_problem_size, m_sel_size, metric_type);
    }else if (((string)model_type)=="GM"){
        m_model=new CGeneralizedMallowsModel(m_problem_size, m_sel_size, metric_type);
    }else if (((string)model_type)=="PL"){
        m_model=new CPlackettLuceModel(m_problem_size, m_sel_size);
    }
 
	m_resultsPath = "results/progress-GM.csv";
    
}

/*
 * The destructor. It frees the memory allocated..
 */
RankingEDA::~RankingEDA()
{
    delete m_best;
    delete m_population;
    delete m_model;
}


/*
 * Experiment for the Metric suitability for the different problems.
 */
void RankingEDA::MetricSuitability_Experiment(char * results_file){

    int tests=1000;
    int samples=1000;
    int * consensus= new int[m_problem_size];
    double thetas_K[10]={1.48,1.8,2.08,2.33,2.6,2.9,3.28,3.75,4.5,10};
    double thetas_C[10]={2.85,3.25,3.56,3.85,4.15,4.47,4.83,5.3,6.1,12};
    double thetas_U[10]={3.3,3.75,4.1,4.4,4.7,5.02,5.4,5.9,6.7,12};
    double theta=0;
    int distance_max_K=(m_problem_size-1)*m_problem_size/2;
    int distance_max_C=m_problem_size-1;
    int distance_max_U=m_problem_size-1;
    int distance=0;
    double fitness_variation=0;
    int fitness_differential=0;
    int distance_max;
    
    int * sample= new int[m_problem_size];
    double total_tests=0;
   // ofstream output_file;
   // output_file.open(results_file);
    for (int j=0;j<10;j++){
        
        if (((string)m_metric_type)=="K")
            theta=thetas_K[j];
        else if (((string)m_metric_type)=="C")
            theta=thetas_C[j];
        else
            theta=thetas_U[j];
        total_tests=0;
        for (int i=0;i<tests;i++){
            fitness_variation=0;
            GenerateRandomPermutation(consensus, m_problem_size);
            m_model->Learn(consensus, theta);
        
            for (int z=0;z<samples;z++){
                m_model->Sample(sample);
                
                //distance of the samples solution with respect to the consensus ranking.
                if (((string)m_metric_type)=="K"){
                    distance = Kendall(consensus,sample,m_problem_size);
                    distance_max=distance_max_K;
                }
                else if (((string)m_metric_type)=="C"){
                    distance = Cayley(consensus,sample,m_problem_size);
                    distance_max=distance_max_C;
                }
                else{
                    distance = Ulam(consensus,sample,m_problem_size);
                    distance_max=distance_max_U;
                }
                
                //fitness differential
                fitness_differential=abs(m_problem->EvaluateInv(sample)-m_problem->EvaluateInv(consensus));
                fitness_variation+=fitness_differential*(1-distance/distance_max);
            }
            fitness_variation=fitness_variation/samples;
            total_tests+=fitness_variation;
            //cout<<fitness_variation<<endl;
        }
        total_tests=total_tests/tests;
        //cout<<total_tests<<endl;
    }
    //cout<<"--------------------------------"<<endl;
 //   output_file.close();
    delete [] consensus;
    delete [] sample;
    
}

/*
 * Running function
 */




int RankingEDA::Run(){
	string results1 = "FileName \t Solution \tFitness \t err \t FEs \n";
	string resultsCsv = "gen,fes,bestFit,avgFit,bestFound,change,changeGen,bestPerChange\n";
	/*long temp = LONG_MIN;
	cout << "temp1" << temp << endl;
	m_best->SetValue(temp);*/

	//cout << "m_evaluations: " << m_evaluations << "m_best: " << m_best->Value() <<  endl;

   // cout<<"Running..."<<endl;
    //variable initializations.
    int i, iChange = 1;
    long int newScore = m_population->m_individuals[0]->Value();
	long int lastScore;
    int * genes= new int[m_problem_size];

    int iterations = 1, changeGen = 1;
    float rate=0.001, average;

	CIndividual *bestChange = new CIndividual(m_problem_size);
	bestChange->SetGenes(m_population->m_individuals[0]->Genes());
	bestChange->SetValue(m_population->m_individuals[0]->Value());


    //EDA iteration. Stopping criterion is the maximum number of evaluations performed
    while (m_evaluations<m_max_evaluations){
		cout << "m_max_evaluations: "<< m_max_evaluations << endl;
		cout << "m_evaluations: " << m_evaluations << endl;
		bool hasChangedOccured = m_problem->changeIdentityPermutation(m_evaluations, m_max_evaluations);
		if (hasChangedOccured) {
			// Re-evaluate population fitness
			for (int i = 0; i<m_pop_size; i++) {
				m_population->SetToPopulation(m_population->m_individuals[i]->Genes(), i, m_problem->Evaluate(m_population->m_individuals[i]->Genes()));

			}
			iChange++;
			changeGen = 1;

			m_population->SortPopulation(1);
			
//			m_best->SetGenes(m_population->m_individuals[0]->Genes());
//			m_best->SetValue(m_population->m_individuals[0]->Value());

			bestChange->SetGenes(m_population->m_individuals[0]->Genes());
			bestChange->SetValue(m_population->m_individuals[0]->Value());

		}
		//average of the population
		average = m_population->AverageFitnessPopulation(m_pop_size);

		cout << "POPULATION:" << endl;
		m_population->Print();

		//results1 += iterations + "," + m_evaluations + "," +newScore + "," + average + "," + m_best->Value() + "," iChange + "," + changeGen + "\n";
		if (((string)m_model_type) == "PL") {
			resultsCsv += to_string(static_cast<long long>(iterations)) + "," + to_string(static_cast<long long>(m_evaluations)) + "," +
				std::to_string(static_cast<long double>((-1) * newScore)) + "," +
				std::to_string(static_cast<long double>((-1) * average)) + "," +
				std::to_string(static_cast<long double>((-1) * m_best->Value())) + "," +
				to_string(static_cast<long long>(iChange)) + "," +
				to_string(static_cast<long long>(changeGen)) + "," +
				std::to_string(static_cast<long double>((-1) * bestChange->Value())) + "\n";
		}

        //learn model
		cout << "LEARNING" << endl;
        m_model->Learn(m_population, m_sel_size);

		/*if (((string)m_model_type) == "PL") {
			((CPlackettLuceModel*)m_model)->printWeigths();
		}*/

        //sample the model.
		cout << "SAMPLING" << endl;
        for (i=0; i< m_offspring_size && m_evaluations<m_max_evaluations; i++){
            //distance at which the sampled solution is from the consensus ranking.
            
            
            m_model->Sample(genes);
			/*PrintArray(genes, m_problem_size, "SAMPLE: ");
			if (isPermutation(genes, m_problem_size)) {
				cout << "TRUE" <<endl;
			}else {
				cout << "FALSE" << endl;
			}*/
           if (m_inverse)
                m_population->AddToPopulation(genes, i, m_problem->EvaluateInv(genes));
            else
                m_population->AddToPopulation(genes, i, m_problem->Evaluate(genes));
            m_evaluations++;
        }

        
        //update the model.
        m_population->SortPopulation(1);


        //update indicators
		lastScore = newScore;
        newScore = m_population->m_individuals[0]->Value();
        float modification = 0;
		if ( newScore > lastScore) {
			bestChange->SetGenes(m_population->m_individuals[0]->Genes());
			bestChange->SetValue(m_population->m_individuals[0]->Value());
		}
        if (newScore > m_best->Value()){
			//results1 = "./test-output/" + "\t" + std::to_string(static_cast<long double>(newScore)) + "\t" + "0" + "\t" + to_string(static_cast<long long>(m_evaluations)) + "\n"; // ##C++0x
			cout << m_evaluations << "; " << newScore << endl;


            m_best->SetGenes(m_population->m_individuals[0]->Genes());
            m_best->SetValue(newScore);
            m_convergence_evaluations=m_evaluations;
           modification = -rate;
        }else{
            modification = rate;
        }
    

        if (((string)m_model_type)=="M"){
            ((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound = ((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound+modification;
            if (((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound < 0.001)
                ((CMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound = 0.001;
        }else if (((string)m_model_type) == "GM") {
            ((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound = ((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound + modification;
            if (((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound < 0.001)
            ((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_lower_theta_bound = 0.001;
        }
       
		//cout << "m_evaluations: " << m_evaluations << "m_best: " << m_best->Value() << "newScore: " << newScore  <<   endl;
        iterations++;
		changeGen++;

    }

	ofstream myfile1;
	myfile1.open(m_resultsPath);
	//	myfile1 << results1;
	//	cout << results << endl;
	myfile1 << resultsCsv;
	myfile1.close();

    delete [] genes;
    return 0;
}

/*
 * Returns the number of performed evaluations.
 */
int RankingEDA::GetPerformedEvaluations(){
    return m_convergence_evaluations;
}

/*
 * Returns the fitness of the best solution obtained.
 */
long int RankingEDA::GetBestSolutionFitness(){
    return m_best->Value();
}

/*
 * Returns the best solution obtained.
 */
CIndividual * RankingEDA::GetBestSolution(){
    return m_best;
}


/*
 * This method applies a swap of the given i,j positions in the array.
 */
void RankingEDA::Swap(int * array, int i, int j)
{
	int aux=array[i];
	array[i]=array[j];
	array[j]=aux;
}
