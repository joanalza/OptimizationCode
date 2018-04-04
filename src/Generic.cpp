//
//  Generic.cpp
//  perms_mallows
//
//  Created by Ekhine Irurozki on 20/06/13.
//  Copyright (c) 2013 Ekhine Irurozki. All rights reserved.
//

#include "Generic.h"
//#include <unistd.h>
#include <fcntl.h>
#include <algorithm>

int**Generic::read_sample_file(int n, int m, char*path){
    ifstream  file;
    int**sample=new int*[m];
    for(int i=0;i<m;i++) sample[i]=new int[n];
    file.open(path);
    if(!((file )))
    {
        cout << "Cannot read input file: "<<path<<endl;
        exit(1);
    }
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            file >> sample[i][j];
    return sample;
}

void Generic::print_long_double_matrix(long double**matrix, int n, int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)cout<<matrix[i][j]<<" ";
        cout<<endl;
    }
}
void Generic::print_int_matrix(int**matrix, int n, int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)cout<<matrix[i][j]<<" ";
        cout<<endl;
    }
}
void Generic::print_double_matrix(double**matrix, int n, int m){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)cout<<matrix[i][j]<<" ";
        cout<<endl;
    }
}
void Generic::print_int_vector(int*vec, int n){
    for(int j=0;j<n;j++)cout<<vec[j]<<" ";
    cout<<endl;
}
void Generic::print_double_vector(double*vec, int n){
    for(int j=0;j<n;j++)cout<<vec[j]<<" ";
    cout<<endl;
}

void Generic::compose(int n, int*s1, int*s2, int*res){
    for(int i = 0 ; i < n ; i++) res[ i ] = s1 [ s2 [ i ] - 1 ];
}
void Generic::invert(int n, int*sigma, int*res){
    for(int i = 0 ; i < n ; i ++) res[ sigma[ i ] - 1 ] = i + 1;
    //for(int j = 0 ; j < n_ ; j ++) samples_inv[i][ samples[ i ][ j ] - 1 ] = j + 1;
}
void Generic::get_permu_matrix(int n,int*sigma, int**matrix){
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) matrix[i][j]=0;
    for(int i=0;i<n;i++) matrix[i][sigma[i]-1]=1;
}
void Generic::print_permus(int n){
    int * permu = new int[ n ];
    for(int i = 0 ; i < n ; i++) permu[ i ] = i + 1;
    sort(permu, permu + n);
    do{
        print_int_vector(permu, n);
    } while (next_permutation(permu, permu + n));
}
long double Generic::factorial(int val) {
    if(val <= 0) return 1;
    //long  N, b, c, p; // use int for fast calculation and small range of calculation..
    long   b, c;
    long double p, N;
    N=(long double)val;
    c = (long)N - 1;
    p = 1;
    while (c > 0) {
        p = 0;
        b = c;
        while (b > 0) {
            if (b & 1) {
                p += N; // p = p + N;
            }
            // if you would like to use double choose the alternative forms instead shifts
            // the code is fast even!
            // you can use the same tips on double or 64 bit int etc.... but you must... ;-)
            //b >>= 1; // b/=2; (b = b / 2;) ( b >> 1; a.s.r. is more efficent for int or long..!)
            b/=2;
            //N <<= 1; // N += N; N = N + N; N = N * 2; (N <<=1; a.s.l. is more efficent for int or long..!)
            N += N;
        } // end of: while(b>0)
        N = p;
        c--; // c = c - 1;
    } // end of: while(c > 0)
    //printf("[%d] is the factorial! \n", p);
    return p;
}
void Generic::init_factorials (int n) {
    facts_n_ = n;
    facts_=new double[n+1];
    facts_[0]=1;
    for(int i=1;i<=n;i++) facts_[i] = facts_[i-1] * i;
}
double Generic::count_permus_no_fixed_points(int n, int k){
    if ( facts_ == NULL ) { init_factorials(n);}
    else if (facts_n_ < n ) {
        cout<<"Check n in Generic::count_permus_no_fixed_points. ";exit(0);}
    double sum=0, aux=0;
    int multi=-1;
    for(int i = 1 ; i <= k ; i++){
        //num = num + (-1)^j * factorial(n-l) * factorial(l) / (factorial(j) * factorial(l-j) );
        aux=multi*facts_[k]*facts_[n - i]/(facts_[i]*facts_[k-i]);
        sum += aux;
        multi *= -1;
    }
    return facts_[n]+sum;
}
double Generic::count_perm_fixed_points(int k, int j){
    if(j<0 || j>k) return 0;
    if(k == 0 && j == 0) return 1;
    return count_perm_fixed_points(k-1, j-1) + count_perm_fixed_points(k - 1, j)*(k - 1  - j)
    + count_perm_fixed_points(k - 1, j + 1)*(j + 1);
    //(1) For every j < 0 or j > k : f(k, j) = 0.
    //(2) f(0, 0) = 1.
    //(3) For every k > 1 and k ≥ j ≥ 0, f(k, j) = f(k − 1, j − 1) + f(k − 1, j)·(k − 1  − j) + f(k − 1, j + 1)·(j + 1)
}

bool Generic::valid_permutation(int*sigma, int n){
    bool*per=new bool[n];
    for(int i=0;i<n;i++)per[i] = false;
    for(int i=0;i<n;i++)
        if(sigma[i] > 0 && sigma[i] <= n && !per[sigma[i]-1] )
            per[sigma[i]-1] = true ;
        else return false;
    delete[]per;
    return true;
}

//void Generic::seed(void) {
//    FILE* fd, buf;
//    if ((fd = fopen("/dev/urandom", O_RDONLY)) < 0) {
//        perror("/dev/urandom");
//    }
//    fread(fd, &buf, sizeof (buf));
//    fclose(fd);
//    srand(buf);
//
//	/*int fd, buf;
//	if ((fd = open("/dev/urandom", O_RDONLY)) < 0) {
//		perror("/dev/urandom");
//	}
//	read(fd, &buf, sizeof(buf));
//	close(fd);
//	srand(buf);*/
//}

void Generic::partition_function_init(int n){
    //F(n,k)=F(n,k−1)+F(n−k,k)
    partition_table=new int*[n+1];
    for(int i=0; i< n+1; i++)
        partition_table[i]=new int[n+1];
    for(int i=0; i< n+1; i++) partition_table[0][i]=1;
    for(int i=0; i< n+1; i++) partition_table[i][0]=0;
    for(int i=1; i< n+1; i++){
        for (int j = 1 ; j < n+1 ; j ++) {
            if (i - j < 0 ) partition_table [i][j] = partition_table[i][j-1];
            else partition_table [i][j] = partition_table[i][j-1] + partition_table[i-j][j];
        }
    }
}

int Generic::get_number_of_partitions(int n){
    return partition_table[n][n];
}

void Generic::partition_function_destructor(int n){
    for(int i=0; i< n+1; i++) delete [] partition_table[i];
    delete [] partition_table;
}