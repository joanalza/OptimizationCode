/*
 *  QAP.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#include "QAP.h"

/*
 *Class constructor.
 */
QAP::QAP()
{
	
}

/*
 * Class destructor.
 */
QAP::~QAP()
{
	for (int i=0;i<m_size;i++)
	{
		delete [] m_distance_matrix[i];
		delete [] m_flow_matrix[i];
	}
	delete [] m_flow_matrix;
	delete [] m_distance_matrix;
	delete [] m_aux;
}



int QAP::Read(string filename)
{
	char line[2048]; // variable for input value
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	int num=0;
	while (!indata.eof())
	{
		//LEER LA LINEA DEL FICHERO
		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
        //cout<<"line: "<<line<<endl;
		if (num==0)
		{
			//OBTENER EL TAMAÑO DEL PROBLEMA
			m_size = atoi(sline.c_str());
			m_distance_matrix = new int*[m_size];
			m_flow_matrix = new int*[m_size];
			for (int i=0;i<m_size;i++)
			{
				m_distance_matrix[i]= new int[m_size];
				m_flow_matrix[i] = new int[m_size];
			}
		}
		else if (1<=num && num<=m_size)
		{
			//LOAD DISTANCE MATRIX
			char * pch;
			pch = strtok (line," ");
			int distance=atoi(pch);
			m_distance_matrix[num-1][0]=distance;
			for (int i=1;i < m_size; i++)
			{
				pch = strtok (NULL, " ,.");
				distance=atoi(pch);
				m_distance_matrix[num-1][i]=distance;
			}
		}
		else if (num>m_size && num<=(2*m_size))
		{
			//LOAD FLOW MATRIX
			char * pch;
			pch = strtok (line," ");
			int weight=atoi(pch);
			m_flow_matrix[num-m_size-1][0]=weight;
			for (int i=1;i < m_size; i++)
			{
				pch = strtok (NULL, " ,.");
				weight=atoi(pch);
				m_flow_matrix[num-m_size-1][i]=weight;
			}
		}
		else
		{
			break;
		}
		num++;
	}
	indata.close();
	m_aux= new int[m_size];
	return (m_size);
}

/*
 * This function evaluates the individuals for the QAP problem.
 */
double QAP::Evaluate(int * genes)
{
	double fitness=0;
	int FactA, FactB;
	int distAB, flowAB, i ,j;
	for (i=0;i<m_size;i++)
	{
		for (j=0;j<m_size;j++)
		{
			FactA = genes[i];
			FactB = genes[j];
			
			distAB= m_distance_matrix[i][j];
			flowAB= m_flow_matrix[FactA][FactB];
			fitness= fitness+(distAB*flowAB);			
		}
	}
	
	return -fitness;
}

/*
 * This function evaluates the inverted solution of the given individual for the QAP problem.
 */
double QAP::EvaluateInv(int * genes)
{
    Invert(genes, m_size,m_aux);
	double fitness=0;
	int FactA, FactB;
	int distAB, flowAB, i ,j;
	for (i=0;i<m_size;i++)
	{
		for (j=0;j<m_size;j++)
		{
			FactA = m_aux[i];
			FactB = m_aux[j];
			
			distAB= m_distance_matrix[i][j];
			flowAB= m_flow_matrix[FactA][FactB];
			fitness= fitness+(distAB*flowAB);
		}
	}
	
	return -fitness;
}

/*
 * Returns the size of the problem.
 */
int QAP::GetProblemSize()
{
    return m_size;
}

