/*
 *  LOP.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "LOP.h"
#include "Tools.h"

/*
 * Class constructor.
 */
LOP::LOP()
{
	
    
}

/*
 * Class destructor.
 */
LOP::~LOP()
{
	for (int i=0;i<m_problemsize;i++)
		delete [] m_matrix[i];
	delete [] m_matrix;
    delete [] m_aux;
}

/*
 * Read LOP instance file.
 */
int LOP::Read(string filename)
{
	char line[2048]; // variable for input value
	string data="";
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	int num=0;
	while (!indata.eof())
	{

		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		if (sline=="")
		{
			break;
		}
		if (num==0)
		{
			m_problemsize = atoi(line);
		}
		else
		{
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		num++;
	}
	indata.close();

	//BUILD MATRIX
	m_matrix = new int*[m_problemsize];
	for (int i=0;i<m_problemsize;i++)
	{
		m_matrix[i]= new int[m_problemsize];
	}
    m_aux= new int[m_problemsize];
    
	istringstream iss(data);
	int i=0;
	int j=0;
	do
	{
		string sub;
	    iss >> sub;
	    if (sub!=""){
			//save distance in distances matrix.
	    	m_matrix[i][j]= atoi(sub.c_str());
	    	if (j==(m_problemsize-1))
	    	{
	    		i++;
	    		j=0;
	    	}
	    	else
	    	{
	    		j++;
	    	}
	    }
	    else
	    {
	    	break;
	    }
	} while (iss);
  /*
    for (int i=0;i<m_problemsize;i++){
        for (int j=0;j<=i;j++){
            cout<<m_matrix[i][j]-m_matrix[j][i]<<" ";
        }
        cout<<endl;
    }
    int aux[100]={69,11,39,29,73,84,45,7,15,43,12,31,56,36,77,22,18,61,58,64,40,1,14,46,8,42,16,51,88,48,9,54,35,50,0,80,74,83,70,96,37,20,34,99,91,97,94,32,90,62,81,26,19,10,59,53,86,82,85,25,68,66,92,47,65,33,30,63,79,98,41,71,2,52,13,44,75,57,93,49,4,23,78,38,87,21,89,17,27,55,28,72,5,3,6,60,67,24,95,76};
    cout<<Evaluate(aux)<<endl;exit(1);*/
	return (m_problemsize);
}

/*
 * This function evaluates the inverted solution of the given solution for the LOP problem.
 */
double LOP::EvaluateInv(int * genes)
{
    Invert(genes,m_problemsize,m_aux); //original
	double fitness=0;
    int i,j;
	for (i=0;i<m_problemsize-1;i++)
		for (j=i+1;j<m_problemsize;j++)
			fitness+= m_matrix[m_aux[i]][m_aux[j]];

	return fitness;
}

/*
 * This function evaluates the solution for the LOP problem.
 */
double LOP::Evaluate(int * genes)
{
	double fitness=0;
    int i,j;
	for (i=0;i<m_problemsize-1;i++)
		for (j=i+1;j<m_problemsize;j++)
			fitness+= m_matrix[genes[i]][genes[j]];
	return fitness;
}

/*
 * Returns the size of the problem.
 */
int LOP::GetProblemSize()
{
    return m_problemsize;
}
