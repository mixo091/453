#pragma once 
#include <string>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include "./Curve.hpp"
#include "./CurveR2.hpp"
#include "./Data/Data.hpp" 
using namespace std;
//Function to parse arguments for nearest neighbor search.
void parse_arguments_search(int argc, char **argv, 
    string *i_file, string *qr_file, int *k, int *L, int *M, 
    int *probes, string * out_file, string * algo, string* metric ,double* d);

//Calulate info of a data Set.
void calculate_file_info(int *total_curves, int *number_of_points, string *filename);

//Function to store the Curve from a file.
void store_curves(string filename, int num_of_points, Curve<double> *arr);
void store_curves(string filename, int num_of_points, Data<double> *arr);
void store_as_curves_r2(string filename, int num_of_points, R2_Curve *arr);
/* Computing and returning vector v of hash function */
template <typename T>
void normal_distribution_fun(vector<T> , float, float);
void normal_distribution_fun(double *, float, float);
/* Easy to understand */
template <typename K>
double euclidean_dist(const K &v1, const K &v2) 
{
    double dist = 0;

    for(unsigned int i = 0; i < v1.size(); i++) 
        dist += pow(v1[i] - v2[i],  2);
    
    return sqrt(dist);
}
/* Modular exponetiation algorithm is used to avoid overflow
 Used in LSH to calculate  Recall (ab) mod m = ((a mod m)(b mod m)) mod m*/
int modular_pow(int, int, int);
unsigned long positive_modulo( unsigned long, unsigned);
/* Simple coin toss function return 0 or 1 for hypercube implementation */
int coinToss();
/* Calculate hamming distance of two verices for hypercube */
int hammingDistance(int, int);

/* This is a brute force method for NN, storing in a vector
* the euclidean distance of a query and every point of the input file
stored in dataset */
template <typename T>
void BruteForceNN(vector<double> qr_v,Data<T> *dataset, int data_size, vector<double> *brute_force_v) {
    // for every vector in dataset, calculate euclidean distance 
    for(int i = 0; i < data_size; i++) {
        double eu_dist = euclidean_dist(qr_v, dataset[i].getVector());

        // push eu_dist to our brute_force_v
        brute_force_v->push_back(eu_dist);
    }
}

// parsing arguments for clustering
void cluster_parse_args(int argc, char **argv, string *input, string *config,
                string *output, string *update, string *assignment, string *complete);