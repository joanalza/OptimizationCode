//
//  Ferrers_diagram.h
//  perms_mallows
//
//  Created by Ekhine Irurozki on 08/07/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#ifndef __perms_mallows__Ferrers_diagram__
#define __perms_mallows__Ferrers_diagram__

#include <iostream>
#include <vector>
#include <cmath>
#include "Generic.h"

class Ferrers_diagram;
using namespace std;
class Ferrers_diagram{
    
private:
    int partition_of_n_;
    long double hook_length_;
    long double number_of_permutations_;
    int*ferrers_shape_;
    int ferrers_shape_length_;
    //int dual_len_;
    //int*dual_;
    std::vector<Ferrers_diagram*>descendants_;
    int desc_length_;
    //Ferrers_diagram*dual_;//transpose or conjugate
    
public:
    Ferrers_diagram(int n, int*partition, int k){
        partition_of_n_ = n;
        ferrers_shape_length_ = k;
        ferrers_shape_=partition;
        hook_length_ = -1;
        number_of_permutations_ = -1;
        desc_length_ = 0;
    }
    ~Ferrers_diagram(){
        delete [] ferrers_shape_;
    }
    
    int get_num_cols(){return ferrers_shape_[0];}
    
    int get_ferrers_shape_length(){return ferrers_shape_length_;}
    int*get_ferrers_shape(){return ferrers_shape_;}
    
    long double get_num_permus(){return number_of_permutations_;}
    
    long double get_hook_len(){return hook_length_;}
    
    std::vector<Ferrers_diagram*>* get_descendants(){return &descendants_;}
    
    long double calculate_hook_length();
    
    void toStr(){
        Generic gen;
        cout<<"Partition of "<<partition_of_n_<<": ";
        gen.print_int_vector(ferrers_shape_, ferrers_shape_length_);
        for(int i=0;i < descendants_.size() ; i++) {
            cout<<"\t descend "<<i<<":\t";
            gen.print_int_vector(descendants_.at(i)->get_ferrers_shape(), descendants_.at(i)->get_ferrers_shape_length());
        }
        cout<<"\t hook length  "<<hook_length_<<"\n";
    }
    
    bool is_parent_of(Ferrers_diagram * f);
    void relate_with(Ferrers_diagram * f);
};

#endif /* defined(__perms_mallows__Ferrers_diagram__) */
