/*
 *  PBP.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#include "PBP.h"
#include "Tools.h"

/*
 * The constructor.
 */
PBP::PBP()
{
	nextChangeIndex = 0;
}

/*
 * The destructor. It frees the memory allocated.
 */
PBP::~PBP()
{
	delete[] identityPerm;
	delete[] identityPermutationChanges ;
	delete[] m_aux;
}


void PBP::setIdentityPermutationChanges(string dynamicfilename) {
	
	//cout << "dynProfilePath" << dynamicfilename << endl;
	//cout << "dynProfilePath2" << dynamicfilename.c_str() << endl;

	std::ifstream infile(dynamicfilename.c_str());
	std::string line;
	int countlines = 0;
	
	while (std::getline(infile, line))
	{
		
		std::istringstream iss(line);
		if (countlines == 0) { 
			//sscanf(line.c_str(), "%d", &nChanges);
			stringstream(line) >> nChanges;
			identityPermutationChanges = new double[nChanges];
			//identityPermutations = new vector<int*>();
		/*}else if (countlines == 1) {
				//sscanf(line.c_str(), "%d", &pbsize);
				stringstream(line) >> pbsize;
				int* temp_array = new int[pbsize];
				for (int i = 0; i < pbsize; i++) {
					temp_array[i] = i;
				}
				identityPerm = temp_array;
				delete[] temp_array;*/
		}else{
			std::string delimiter = ";";
			int point1 = line.find(delimiter);
			std::string part0 = line.substr(0, point1);
			std::string part1 = line.substr(point1+1,line.length());
			//sscanf(part0.c_str(), "%d", &identityPermutationChanges[countlines - 2]);
			stringstream(part0) >> identityPermutationChanges[countlines - 1];
			size_t pos = 0;
			std::string delimiter2 = ",";
			std::string token;
			int* temparray = new int[pbsize];
			int i = 0;
			while ((pos = part1.find(delimiter2)) != std::string::npos) {
				token = part1.substr(0, pos);
				//std::cout << token << std::endl;
				//sscanf(token.c_str(), "%d", &temparray[i]);
				stringstream(token) >> temparray[i];
				i++;
				part1.erase(0, pos + delimiter.length());
			}
			stringstream(part1) >> temparray[i];
			/*cout << "temparray ..." << endl;
			for (int i = 0; i < pbsize; i++) {
				cout << temparray[i] << " ";
			}*/
			//cout << "remaining: " << part1 << endl;
			identityPermutations.push_back(temparray);
		}
		//cout << "line" << line << endl;
		countlines++;
	}

	// Initialize identity permutation
	identityPerm = new int[pbsize];
	for (int i = 0; i < pbsize; i++) {
		identityPerm[i] = i;
	}

}

void PBP::setPbsize(int size) {
	this->pbsize = size;
}

bool PBP::changeIdentityPermutation(int fes, int maxfes) {
	//    	System.out.println("fes: "+fes+" - maxfes: "+maxfes);
	//    	System.out.println("this.nextChangeIndex: "+this.nextChangeIndex);
	//    	System.out.println("this.identityPermutationChanges.length: "+ this.identityPermutationChanges.length);
	bool hasChangedOccured = false;
	if (nextChangeIndex < nChanges) {
		/*cout << "nextChangeIndex" << nextChangeIndex << endl;
		cout << "fes" << fes << endl;
		cout << "maxfes" << maxfes << endl;
		cout << "identityPermutationChanges[nextChangeIndex]" << identityPermutationChanges[nextChangeIndex] << endl;
		cout << "identityPermutationChanges[nextChangeIndex] * (double)maxfes" << identityPermutationChanges[nextChangeIndex] * (double)maxfes << endl;*/
		int nextChangeFes = (int)rint(identityPermutationChanges[nextChangeIndex] * (double)maxfes);
		//cout << "fes - nextChangeFes" << fes << " - " << nextChangeFes << endl;
		if (fes >= nextChangeFes) {
			/*cout << "changing id perm ..." << endl;
			for (int i = 0; i < pbsize; i++) {
				cout << identityPermutations[nextChangeIndex][i] << " ";
			}
			cout << endl;
			cout << identityPermutations[nextChangeIndex] << endl;*/
			identityPerm = identityPermutations[nextChangeIndex];
			nextChangeIndex++;
			PrintArray(identityPerm, pbsize, "Identity: ");
			//    			System.out.println("CHANGE OF IDENTITY PERMUTATION: "+ArrayUtils.tableToString(this.identityPerm));
			hasChangedOccured = true;
		}
	}
	return hasChangedOccured;
}
