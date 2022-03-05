include <iostream>
#include <random>
#include <vector>
#include "../headers/HashFunction.hpp"



using namespace std;

template <typename T>
int HashFunction<T>::hfunction(vector<T> p) {
    long double product = 0;

    // create a vector by normal distribution 
    vector<double> v = vector<double>(this->dim);
    normal_distribution_fun(v, 0.0, 1.0);

    for(int i = 0; i < this->dim; i++)
        product += v[i] * p[i];
    
    return (int) floor( ( this->t + (double) product) / this->w);

}

template <typename T>
int HashFunction<T>::hashValue(vector<T> p, int tableSize) {
    int numberOfHashFunctions = this->k;
    // this will be our Sum of h * r
    int sum = 0;

    for(int i = 0; i < numberOfHashFunctions; i++) {
        // compute h(p)
        this->h[i] = hfunction(p);

        // we need to avoid overflow, so we use the formula from the slides
        // Recall (ab) mod m = ((a mod m)(b mod m)) mod m
        int res = modular_pow(this->h[i], this->r[i], M);
        
        if(res <= 0) {
            cerr << "Modular error" << endl;
            exit(-1);
        }

        sum += res;
    }

    return sum % tableSize;
}