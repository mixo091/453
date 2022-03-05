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
        cout << "number of hash tables " << numberOfHashTables << endl;
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
        map<double, string> my_map;
        map<double, string> result_map;

        double avg_duration_lsh = 0;
        double avg_duration_brute_force = 0;
        vector<double> vector_MAF;

        for(int i = 0; i < qr_lines; i++) {

            auto start_lsh = high_resolution_clock::now();

            // get every query point
            // and hash it in every hash table
            for(int j = 0; j < numberOfHashTables; j++) {
                hash_tables[j]->find_NN(&qr_data[i], my_map, N);

                for (auto it = my_map.cbegin(); it != my_map.cend(); ++it) {
                    result_map.insert(pair<double, string>((*it).first, (*it).second));    
                }
                my_map.clear();
            }
        
            auto stop_lsh = high_resolution_clock::now();
            auto duration_lsh = duration_cast<microseconds>(stop_lsh - start_lsh);
            avg_duration_lsh +=  duration_lsh.count();

            auto start_bf_search = high_resolution_clock::now();

            // vector used for brute force
            vector<pair<double, string>> brute_force_v;

            // Brute force method for NN
            BruteForceNN(qr_data[i].getVector(), in_data, in_dataSize, &brute_force_v);

            // sort brute force vector
            sort(brute_force_v.begin(), brute_force_v.end());

            auto stop_bf_search = high_resolution_clock::now();
            auto duration_bf = duration_cast<microseconds>(stop_bf_search - start_bf_search);
            avg_duration_brute_force += duration_bf.count();

            output << "Query: " << qr_data[i].get_name();
            output<< "\nAlgorithm:  { " << "LSH }";
            output << "\nApproximate Nearest neighbor:  ";
            if(result_map.size() > 0) {
                output << result_map.cbegin()->second;
            }
            output << "\nTrue Nearest neighbor:    " ; 
            if(brute_force_v.size() > 0) output << brute_force_v.at(0).second;
            output << "\ndistanceApproximate:    "; 
            if(result_map.size() > 0) {
                output << result_map.cbegin()->first; 
            }
            output << "\ndistanceTrue:  "; 
            if(brute_force_v.size() > 0) output << brute_force_v.at(0).first;
            output << endl << endl;

            if(result_map.size() > 0 && brute_force_v.size() > 0) {
                vector_MAF.push_back(result_map.cbegin()->first / brute_force_v.at(0).first);
            }

            result_map.clear();
            my_map.clear();
            brute_force_v.clear();
        }
        

        output << "\n\ntApproximateAverage:    " << (double) (avg_duration_lsh / 1000000) / qr_lines;
        output << "\ntTrueAverage:    " << (double) (avg_duration_brute_force / 1000000) / qr_lines; 
        output << "\nMAF:    " <<  *max_element(vector_MAF.begin(), vector_MAF.end()) << endl << endl;

        vector_MAF.clear();

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
