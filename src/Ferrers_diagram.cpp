//
//  Ferrers_diagram.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 08/07/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include "Ferrers_diagram.h"

bool Ferrers_diagram::is_parent_of(Ferrers_diagram * f){
    int dif = 0;
    if(f->get_ferrers_shape_length() == ferrers_shape_length_ ){
        for(int i=0 ; dif <=1 && i < ferrers_shape_length_ ; i++)
            dif += abs( ferrers_shape_[i] - f->get_ferrers_shape()[i]);
        return (dif == 1);
    } else if (f->get_ferrers_shape_length() + 1 == ferrers_shape_length_ ){
        for(int i=0 ; i < ferrers_shape_length_-1 ; i++)
            dif += abs( ferrers_shape_[i] - f->get_ferrers_shape()[i]);
        return (dif == 0 && ferrers_shape_[ferrers_shape_length_-1] == 1);
    }
    return false;
}

void Ferrers_diagram::relate_with(Ferrers_diagram * f){
    descendants_.push_back(f);
    desc_length_ ++;
}

long double Ferrers_diagram::calculate_hook_length(){// optim . is very slow, desc_length_ no reall y useful
    bool positive = true;
    if(ferrers_shape_length_ == 2 && ferrers_shape_[0] == 1 && ferrers_shape_[1]==1)
        bool trace = true;
    if(hook_length_ == -1){
        //cout<<"calculando hooks "<<partition_of_n_<<", with "<<desc_length_ <<" desc"<<endl;
        if(desc_length_ == 0) hook_length_ = 1;
        else{
            long double sum = 0;
            for ( int i = 0 ; i < desc_length_; i++){
                sum += descendants_[i]->calculate_hook_length();//descendants_.at(i)->calculate_hook_length();
                if(descendants_[i]->calculate_hook_length()<0 ) positive = false;
            }
            hook_length_ = sum;
        }
        number_of_permutations_ = hook_length_ * hook_length_;
    }
    //else{cout<<"just ret" <<endl;}
    if(hook_length_<0 && positive){
        toStr();
        for ( int i = 0 ; i < desc_length_; i++){
            cout<<"descend! "<<endl;
            descendants_[i]->toStr();
    }
    }
    return hook_length_;
}


