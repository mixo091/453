#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <ctime>
#include <chrono>

// #include "../Data/Data.hpp"
#include "../hashTable/HashTable.hpp"
#include "../HashFunction.hpp"
#include "../Functionality.hpp"

#define BUCKET_DIVIDER 8

using namespace std;
using namespace std::chrono;

template <typename T>
class Lsh
{
protected:
    int numberOfHashTables;
    int ht_Size;
    int vecDimension;
    int numberOfHashFunctions;
    int w;
    
   HashTable<Data<T> *> **hash_tables; 

public:
    
    Lsh(int L, int totalVectors, int _dim, int _k, int _w, Data<T> *data)
    :   numberOfHashTables(L), 
        ht_Size(totalVectors / BUCKET_DIVIDER), 
        vecDimension(_dim), 
        numberOfHashFunctions(_k), 
        w(_w)
    {   
        /* this constructor is called for LSH */
        hash_tables = new HashTable<Data<T> *>*[numberOfHashTables];
        for(int i = 0; i < numberOfHashTables; i++) {
            hash_tables[i] = new HashTable<Data<T>*>(ht_Size, w, numberOfHashFunctions, vecDimension);
        }
        // let's insert some data from the dataset in hash table
        // with the use of hash function
        for(int i = 0; i < totalVectors; i++) {
            for(int j = 0; j < numberOfHashTables; j++) 
                hash_tables[j]->insert(&data[i]);

        }
    }

    Lsh(int L, int _size, int _dim, int _k, int _w)
    :   numberOfHashTables(L), 
        ht_Size(_size), 
        vecDimension(_dim), 
        numberOfHashFunctions(_k), 
        w(_w)
    {   
        /* this constructor is called for Hypercube */
        // cout << "Constructing Lsh...\n";

        // create hash tables
        hash_tables = new HashTable<Data<T> *>*[numberOfHashTables];
        hash_tables[0] = new HashTable<Data<T>*>(ht_Size, w, numberOfHashFunctions, vecDimension);        
    }

    void PrintHTs()
    {
        /* printing table */
        for(int j = 0; j < numberOfHashTables; j++){
            cout <<"__ HASH TABLE ["<<j<<"] __"<<endl;
            hash_tables[j]->PrintHT();
        }
    }

    void ANN(Data<T> *qr_data, int qr_lines, Data<T> *in_data, int in_dataSize, 
            int N, const string &out_file, float radius) 
    {   
        /* 
         * Aproximate nearest neighbor function that calls find_NN function
         * of HashTable and stores in a map the best k pair of euclidean distance and ids.
         * We also find R-nearest neighbors.
         * Finally we prin some stats.
         */

        // open output file
        ofstream output(out_file);
        if(!output.is_open()) {
            cout << "Error with output file\n";
            exit(-1);
        }

        // create a map of <id, eu_dist> where id is int and eu_dist is double
        map<double, int> my_map;
        map<double, int> result_map;

        // used to store for every query R-near neighbors
        std::vector<int> results_of_radius_nearest_neighbors_vec;
        std::vector<int> temp_radius_nearest_neighbors;

        for(int i = 0; i < qr_lines; i++) {

            auto start_lsh = high_resolution_clock::now();

            // get every query point
            // and hash it in every hash table
            for(int j = 0; j < numberOfHashTables; j++) {
                hash_tables[j]->find_NN(&qr_data[i], my_map, N);

                for (auto it = my_map.cbegin(); it != my_map.cend(); ++it) {
                    result_map.insert(pair<double, int>((*it).first, (*it).second));    
                }
                my_map.clear();

                // call range search
                temp_radius_nearest_neighbors = hash_tables[j]->range_search(&qr_data[i], radius);

                std::copy(temp_radius_nearest_neighbors.begin(), 
                        temp_radius_nearest_neighbors.end(), 
                        std::back_inserter(results_of_radius_nearest_neighbors_vec));
            }
        
            auto stop_lsh = high_resolution_clock::now();
            auto duration_lsh = duration_cast<microseconds>(stop_lsh - start_lsh);

            auto start_bf_search = high_resolution_clock::now();

            // vector used for brute force
            vector<double> brute_force_v;

            // Brute force method for NN
            BruteForceNN(qr_data[i].getVector(), in_data, in_dataSize, &brute_force_v);

            // sort brute force vector
            sort(brute_force_v.begin(), brute_force_v.end());

            auto stop_bf_search = high_resolution_clock::now();
            auto duration_bf = duration_cast<microseconds>(stop_bf_search - start_bf_search);

            int n_ostos = 0;
            output << "\nQuery: " << qr_data[i].getId();
            for(auto it = result_map.cbegin(); it != result_map.cend(); ++it) {
                
                output << "\nNearest neighbor-" << n_ostos + 1 << ": " << (*it).second; 
                output << "\ndistanceLsh: " << (*it).first; 
                output << "\ndistanceTrue: " << brute_force_v.at(n_ostos);
                output << "\ntLsh: " << duration_lsh.count();
                output << "\ntTrue:" << duration_bf.count();
                
                if(++n_ostos == N) {
                    output << "\nR-near neighbors :" << endl;

                    for(auto vec_it = results_of_radius_nearest_neighbors_vec.cbegin(); vec_it != results_of_radius_nearest_neighbors_vec.cend(); ++vec_it) {
                        output << *vec_it << endl;
                    }

                    if(results_of_radius_nearest_neighbors_vec.size() == 0) output << "Not found\n"; 

                    break;
                }
            }
            result_map.clear();
            my_map.clear();
            brute_force_v.clear();
            temp_radius_nearest_neighbors.clear();
            results_of_radius_nearest_neighbors_vec.clear();
        }
        // close output file
        output.close();
    }



    vector<int> ReverseAssignment(vector<double>centroid,int id ,double R){
        //Result items in Range
        std::vector<int> results_of_radius_nearest_neighbors_vec;
        //in Range in hash Table
        std::vector<int> temp_radius_nearest_neighbors;

        for(int j = 0; j < numberOfHashTables; j++) {
            temp_radius_nearest_neighbors = hash_tables[j]->clustering_range_search(centroid,id, R);

            std::copy(temp_radius_nearest_neighbors.begin(), 
                temp_radius_nearest_neighbors.end(), 
                std::back_inserter(results_of_radius_nearest_neighbors_vec)
                );
        }
        
        return results_of_radius_nearest_neighbors_vec;
    }

    virtual ~Lsh() {
        // cout << "Destructing Lsh... " << numberOfHashTables << endl;
        for (int i = 0; i < numberOfHashTables; i++) {
            delete hash_tables[i];
        }

        delete [] hash_tables;

    }
};