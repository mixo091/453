#pragma once
#pragma once

#include <iostream>
#include <string.h>
#include <chrono>
#include <random>
#include <climits>
#include <math.h>
#include <vector>
#include <assert.h>     /* assert */
#include "./Functionality.hpp"

using namespace std;

class HashFunction
{
private:
    double *v; // Vector array .
    int w; // needed for hash value .
    double t;
    int num_of_hfun;
    int dim;
    int M;
    int *r; // array  of r values for g calculation.
    int *h;  // array of h values . 
public:
    HashFunction(int _w, int _k, int _dim) : w(_w), num_of_hfun(_k), dim(_dim) 
    {
        unsigned seed = chrono::steady_clock::now().time_since_epoch().count(); 
        default_random_engine e (seed); 
        /* declaring normal distribution object 'distN' and initializing its mean and standard deviation fields. */
        /* Mean and standard deviation are distribution parameters of Normal distribution. Here, we have used mean=5, and standard deviation=2. You can take mean and standard deviation as per your choice */
        uniform_int_distribution<int> distN(-20, 20);
        r = new int[num_of_hfun];
        h = new int[num_of_hfun];
        for(int i = 0; i < num_of_hfun; i++) 
            r[i] = distN(e);
        // set t
        normal_distribution_fun(&t, 0.0, (double) w);
        // calculate M , needed for modular operation
        M = (int) (pow(2, 32 ) - 5);
    }

    ~HashFunction() {
        delete [] r; delete [] h;
    }

    int hfunction(vector<double> p) {
        long double product = 0;
        vector<double> v;

        // create a vector by normal distribution 
        unsigned seed = chrono::steady_clock::now().time_since_epoch().count(); 
        default_random_engine e (seed);
        std::default_random_engine generator;
        normal_distribution<double> distN(0.0, 1.0);

        for(int i = 0; i < dim; i++)
            {v.push_back(distN(e));}

        for(int i = 0; i < dim; i++)
            product += v[i] * p[i];
         
        return (int) floor( ( this->t + (double) product) / this->w);
        


    }

    int hashValue(vector<double> p, int tableSize, int *ID) { 
        int amodM;
        int bmodM;
        for(int i = 0; i < num_of_hfun; i++) {
            h[i] = hfunction(p);
            // /avoid overflow
            amodM = positive_modulo(h[i], M);
            bmodM = positive_modulo(r[i], M);
            int res = amodM * bmodM;

            res = positive_modulo(res, M);
            *ID += res;
        }
        int g = positive_modulo(*ID, tableSize);
        return g;
    }
};