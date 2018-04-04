//
//  Generic.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Generic__
#define __perms_mallows__Generic__


#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>

class Generic;
using namespace std;
class Generic{
public:
    Generic(){
        facts_ = NULL;}
    
    ~Generic(){if ( facts_ == NULL) delete [] facts_;}
    
    int**read_sample_file(int n, int m, char*path);
    
    void print_long_double_matrix(long double**matrix, int n, int m);
    
    void print_double_matrix(double**matrix, int n, int m);
    
    void print_int_matrix(int**matrix, int n, int m);
    
    void print_int_vector(int*vec, int n);
    
    void print_double_vector(double*vec, int n);
    
    void generate_random_permutation(int*sigma, int len, int first_item_in_perm);
    
    void compose(int n, int*s1, int*s2, int*res);
    
    void invert(int n, int*sigma, int*res);
    
    void get_permu_matrix(int n,int*sigma, int**matrix);
    
    void print_permus(int n);
    
    long double factorial(int val) ;

    double count_permus_no_fixed_points(int n, int k);
    
    double count_perm_fixed_points(int k, int j);
    
    bool valid_permutation(int*sigma, int n);
    
    void seed(void);
    
    void partition_function_init(int n);
    
    int get_number_of_partitions(int n);
    
    void partition_function_destructor(int n);
    
    void init_factorials (int n) ;
    
    void delete_factorials();
    
    double get_factorial_in_table(int i){return facts_[i];}
    
private:
    
    int**partition_table;// matrix for the table of partition function (number of integer partitions of n)
    
    double *facts_;
    
    int facts_n_;
};

#endif /* defined(__perms_mallows__Generic__) */
