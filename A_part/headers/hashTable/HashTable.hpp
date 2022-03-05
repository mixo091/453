#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <bits/stdc++.h>

#include <assert.h>
#include <ctime>
#include <map>

#include "../Data/Data.hpp"
#include "../HashFunction.hpp"

// for frechet
#include "../frechet/frechet.hpp"
#include "../frechet/curve.hpp"
#include "../frechet/types.hpp"
#include "../frechet/point.hpp"

using namespace std;

using namespace Frechet;
using namespace Continuous;

/* A bucket is a list of points in space 
 In this case as in LSH we keep the points many times
 we will have a bucket as a list of pointers in points (vectors) */

template <typename K>
struct BucketEntry {
    /* ID: used from theory as Query Trick to avoid calculating Euclidean
    * distance, betwween query and every point in Bucket */
    int query_trick;
    /* This is our key */
    K key;

    BucketEntry(K k, int ID) {
        key = k;
        query_trick = ID;
    }

    ~BucketEntry() {

    }

    void print() {
        Data<double> *data = (Data<double> *) key;
        data->printVector();
        cout << "ID(P) = " << query_trick << endl << endl;
    }
    int getQueryTrick() { return query_trick; }
    vector<double> getVector() { return key->v; }
    int get_mark_status (){ return key->assigned_to;}
    void mark(int clustered_to){key->assigned_to = clustered_to;}
    int getId() { return key->id; }
    string get_name() { return key->get_name(); }
};

/* This is the bucket of our hash table */
template <typename K>
struct Bucket {
    int number_of_entries;
    int id;
    list<BucketEntry<K> *> bucket_entries; 

    Bucket(int i) {
        id = i;
        number_of_entries = 0;
    }

    ~Bucket() {
        // cout << "Destructing Bucket...";
        for(typename list<BucketEntry<K>*>::const_iterator it = bucket_entries.begin(); it != bucket_entries.end(); ++it)
        {
            delete *it;
        } 
        bucket_entries.clear();
    }

    void insert(K key, int ID){
        bucket_entries.push_back(new BucketEntry<K>(key, ID));
        // increase list items
        number_of_entries++;
    }

    void displayBucket(){
        // cout<<"BUCKET ["<<id<<"] :"<<endl;
        if(!bucket_entries.empty())
            for(auto &it : bucket_entries)
                it->print();

    }
};

template <typename K>
class HashTable
{   
    /* Number of buckets */
    int buckets;
    /* Hash table size */
    int table_size;
    struct Bucket<K> **hash_table ; 
    /* Every hash table has it's own hash function */
    HashFunction *h_fun;

public:
    HashTable(int ht_size, int w, int k, int dim) 
    : table_size(ht_size) 
    {
        // cout << "Constructing hash table with size " << table_size << endl;
        // allocate the buckets.
        hash_table = new Bucket<K>*[table_size];

        for(int i = 0; i < table_size; i++)
            hash_table[i] = NULL; 
        //Creating the hash function
        h_fun = new HashFunction(w, k, dim);
    }

    HashFunction *getHashFunction() { return h_fun; }

    // query_trick = (h_1 * r_1 ) mod M + (h2 * r_2) mod M.......

    /* used for continues frechet */
    void insert(vector<double> key, Data<double> *value) {
        // this is used from theory as query trick
        int query_trick = 0;
        unsigned int index = h_fun->hashValue(key, table_size, &query_trick);
       // check size of index
        assert(index <= INT_MAX);

        if(hash_table[index] != NULL) {
           // key->PrintID();
            //cout<<"ASSIGNED TO BUCKET"<<index<<endl;
            hash_table[index]->insert(value, query_trick); 
        }
        else {
            // create new Bucket and insert key
           // key->PrintID();
            //cout<<"ASSIGNED TO BUCKET"<<index<<endl;
            hash_table[index] = new Bucket<K>(index);
            hash_table[index]->insert(value, query_trick);
        }
    }

    void insert(Data<double> *key) {
        /* insertion for lsh*/ 

        // this is used from theory as query trick
        int query_trick = 0;
        unsigned int index = h_fun->hashValue(key->getVector(), table_size, &query_trick);
       // check size of index
        assert(index <= INT_MAX);

        if(hash_table[index] != NULL) {
           // key->PrintID();
            //cout<<"ASSIGNED TO BUCKET"<<index<<endl;
            hash_table[index]->insert(key, query_trick); }
        else {
            // create new Bucket and insert key
           // key->PrintID();
            //cout<<"ASSIGNED TO BUCKET"<<index<<endl;
            hash_table[index] = new Bucket<K>(index);
            hash_table[index]->insert(key, query_trick);
        }
    }

    void insertHyperCube(Data<double> *key, int index) {
        /* insertion for hypercube */

        if(hash_table[index] != NULL) hash_table[index]->insert(key, -1);
        else {
            hash_table[index] = new Bucket<K>(index);
            hash_table[index]->insert(key, -1);
        }
    }

    void search_NN_in_radius(Data<double> *q, int current_bucket_num,
                            int probes, int *count_probes_searched,
                            int M, int *count_items_searched,
                            vector<int> &v_of_neighbors_in_radius, float radius)
    {
        /* 
        * Helper searching function to find nearest neighbors in radius.
        * Used in hypercube.
        */
        bool stop_searching = false;
        std::vector<int>::iterator vec_it;

        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {

            // traverse the list of Data at current bucket.
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
                    
                double euDist = euclidean_dist(q->getVector(), (*it)->getVector());
                
                if(euDist < radius) {
                    // std::find function call
                    vec_it = std::find (v_of_neighbors_in_radius.begin(), v_of_neighbors_in_radius.end(), (*it)->getId());
                    if (vec_it == v_of_neighbors_in_radius.end())
                        v_of_neighbors_in_radius.push_back((*it)->getId());
                }

                if(++(*count_items_searched) >= M) {
                    stop_searching = true;
                    break;
                } 
            }
            // check the bucket in which query is hashed
            if(++(*count_probes_searched) >= probes) stop_searching = true;

            if(stop_searching) return;
        }
    }



    void search_NN_in_radius_cluster(vector<double> *q,int id , int current_bucket_num,
                            int probes, int *count_probes_searched,
                            int M, int *count_items_searched,
                            vector<int> &v_of_neighbors_in_radius, double R)
    {
        /* 
        * Helper searching function to find nearest neighbors in radius.
        * Used in hypercube.
        */
        bool stop_searching = false;
        std::vector<int>::iterator vec_it;

        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {

            // traverse the list of Data at current bucket.
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
                double euDist = euclidean_dist(*q, (*it)->getVector());       

                if( ((*it)->key->Clustered_already_by== -1) || 
                    ( ((*it)->key->ClusteredInRadius == R)  && ((*it)->key->Clustered_already_by != id) ) )
                {
                    if(euDist < R)
                    {  
                        (*it)->key->Clustered_already_by = id;
                        (*it)->key->ClusteredInRadius = R;
                        vec_it = std::find (v_of_neighbors_in_radius.begin(), v_of_neighbors_in_radius.end(), (*it)->getId());
                        if (vec_it == v_of_neighbors_in_radius.end())
                            v_of_neighbors_in_radius.push_back((*it)->getId());
                    }
                }
                if(++(*count_items_searched) >= M) {
                    stop_searching = true;
                    break;
                } 
            }
            // check the bucket in which query is hashed
            if(++(*count_probes_searched) >= probes) stop_searching = true;

            if(stop_searching) return;
        }
    }

    void search_NN_in_radius(Data<double> *q, vector<int> &v_of_neighbors_in_radius, float radius)
    {
        /* 
        * Helper searching function to find nearest neighbors in radius.
        * Used in LSH.
        */
        std::vector<int>::iterator vec_it;
        int query_trick = 0;
        int current_bucket_num = h_fun->hashValue(q->getVector(), table_size, &query_trick);

        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
                    
                double euDist = euclidean_dist(q->getVector(), (*it)->getVector());                
                if(euDist < radius) {
                    // std::find function call
                    vec_it = std::find (v_of_neighbors_in_radius.begin(), v_of_neighbors_in_radius.end(), (*it)->getId());
                    if (vec_it == v_of_neighbors_in_radius.end())
                        v_of_neighbors_in_radius.push_back((*it)->getId());
                }
            }
        }
    }

    void search_NN(Data<double> *q, map<double, string> &k_nearest_map, int k)
    {
        /* 
        * helper search function used in LSH in order to find k nearest neighbors
        * and store them in a map
        */
        map<double, string> distance_map;
        int query_trick = 0;

        int current_bucket_num = h_fun->hashValue(q->getVector(), table_size, &query_trick);

        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {

            // traverse the list of Data at current bucket.
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
            
                if(query_trick == (*it)->getQueryTrick()) {
                    double euDist = euclidean_dist(q->getVector(), (*it)->getVector());
                    distance_map.insert(pair<double, string>(euDist, (*it)->get_name()));
                } else 
                    continue;           
            }

            // keep the k best from the bucket.
            int item = 0 ;
            for (auto it = distance_map.cbegin(); it != distance_map.cend(); ++it) {
                k_nearest_map.insert(pair<double, string>((*it).first,(*it).second));

                // if k reached, stop adding new neigbor
                if(++item >= k) break;
            }
            distance_map.clear();   
        }
    }

    void search_NN(Data<double> *q, int current_bucket_num, int M, int probes, int k,
                    int *count_items_searched, int *count_probes_searched,
                    map<double, string> &k_nearest_map)
    {
        /*  
        * helper function to find k nearest neighbors and store them in a map
        * used in hypercube
        */
        map<double, string> distance_map;
        bool stop_searching = false;

        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {

            // traverse the list of Data at current bucket.
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
                    
                double euDist = euclidean_dist(q->getVector(), (*it)->getVector());
                distance_map.insert(pair<double, string>( euDist, (*it)->get_name())); 

                if(++(*count_items_searched) >= M) {
                    stop_searching = true;
                    break;
                }                     
            }
            // check the bucket in which query is hashed
            if(++(*count_probes_searched) >= probes) stop_searching = true;

            // keep the k best from the bucket.
            int item = 0 ;
            for (auto it = distance_map.cbegin(); it != distance_map.cend(); ++it) {
                k_nearest_map.insert(pair<double, string>((*it).first,(*it).second));

                // if k reached, stop adding new neigbor
                if(++item >= k) break;
            }
            distance_map.clear();   

            // stop searching
            if(stop_searching) return; 
        }
    }

    void find_NN(Data<double> *q, map<double, string> &k_nearest_map,
                int k, vector<int> neighbors, int current_bucket_num, 
                int M, int probes, int *count_items_searched, int *count_probes_searched) 
    {
        /* find NN for Hypercube checking current_bucket_num at first and if needed, every other neighbor*/

        bool stop_searching = false;
        // at first search is done for current bucket
        search_NN(q, current_bucket_num, M, probes, 
                k, count_items_searched, 
                count_probes_searched, k_nearest_map);

        if(*count_items_searched >= M|| *count_probes_searched >= probes) stop_searching = true;

        // search is done for every neighbor
        for(int i = 0; i < (int)neighbors.size() && stop_searching == false; i++) {

            search_NN(q, neighbors.at(i), M, probes, k,
                    count_items_searched, count_probes_searched, 
                    k_nearest_map);

            if(*count_items_searched >= M|| *count_probes_searched >= probes) stop_searching = true;
        }
    }

    
    void find_NN(Data<double> *q, map<double, string> &k_nearest_map, int k) 
    {
        /* NN for lsh. We calculate euclidean distance with the use of query trick 
        * and we keep the best-k pair of <eu_dist, id> from every bucket  in  map
        */
        search_NN(q, k_nearest_map, k);   
        
    }

    vector<int> clustering_range_search(vector<double> centroid,int id,double R){
        // clustering range search used for lsh
        vector<int> v_of_neighbors_in_radius;
        std::vector<int>::iterator vec_it;
        int query_trick = 0;
        int current_bucket_num = h_fun->hashValue(centroid, table_size, &query_trick);
        // get hashed bucket
        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {
            // iterate through the list of bucket
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
                // calculate eu distsance
                double euDist = euclidean_dist(centroid, (*it)->getVector()); 
                if( ((*it)->key->Clustered_already_by == -1) || 
                    ( ((*it)->key->ClusteredInRadius == R)  && ((*it)->key->Clustered_already_by != id) ) )
                {
                    if(euDist < R) 
                    {              
                        (*it)->key->Clustered_already_by = id;
                        (*it)->key->ClusteredInRadius = R;
                        vec_it = std::find (v_of_neighbors_in_radius.begin(), v_of_neighbors_in_radius.end(), (*it)->getId());
                        if (vec_it == v_of_neighbors_in_radius.end())
                            v_of_neighbors_in_radius.push_back((*it)->getId());
                    }
                }
            }   
        }
         return v_of_neighbors_in_radius;
    }

    vector<int> range_search(Data<double> *data, float radius) 
    {
        /*
        * Performing range search for LSH.
        * Returnin a vector of id's with eu_dist < radius.
        */
        vector<int> rad_nearest_neighbors;
        search_NN_in_radius(data, rad_nearest_neighbors, radius);

        return rad_nearest_neighbors;
    }


     vector<int> clustering_range_search(vector<double>*  centroid,int id,double R,int bucket_num, 
                    int probes, int *count_probes_searched, 
                    int M, int *count_items_searched,vector<int> neighbors){

        vector<int> v_of_neighbors_in_radius;
        std::vector<int>::iterator vec_it;
        vector<int> rad_nearest_neighbors;
        bool stop_searching = false;
        
        search_NN_in_radius_cluster(centroid,id, bucket_num,
                        probes, count_probes_searched,
                        M, count_items_searched,
                        rad_nearest_neighbors, R);

        
        if(*count_items_searched >= M|| *count_probes_searched >= probes) stop_searching = true;

        // range search is done for every neighbor
        for(size_t i = 0; i < neighbors.size() && stop_searching == false; i++) {

            search_NN_in_radius_cluster(centroid,id, bucket_num,
                probes, count_probes_searched,
                M, count_items_searched,
                rad_nearest_neighbors, R);

            if(*count_items_searched >= M|| *count_probes_searched >= probes) stop_searching = true;
        }

        return rad_nearest_neighbors;
    }

    vector<int> range_search(Data<double> *data, int bucket_num, 
                    int probes, int *count_probes_searched, 
                    int M, int *count_items_searched, 
                    float radius, vector<int> neighbors)
    {
        /*  
        *   Performing range search, adding every vertex with distance < radius in a vector of distances.
        *   Returing a vector of id's.
        *   Used in hypeprcube.
        */
        vector<int> rad_nearest_neighbors;
        bool stop_searching = false;
        
        search_NN_in_radius(data, bucket_num,
                        probes, count_probes_searched,
                        M, count_items_searched,
                        rad_nearest_neighbors, radius);

        
        if(*count_items_searched >= M|| *count_probes_searched >= probes) stop_searching = true;

        // range search is done for every neighbor
        for(size_t i = 0; i < neighbors.size() && stop_searching == false; i++) {

            search_NN_in_radius(data, bucket_num,
                        probes, count_probes_searched,
                        M, count_items_searched,
                        rad_nearest_neighbors, radius);

            if(*count_items_searched >= M|| *count_probes_searched >= probes) stop_searching = true;
        }

        return rad_nearest_neighbors;
    }

    void search_NN(vector<double> key, Data<double> *query, pair<double, string> *nn_pair, double *duration)
    {
        /* 
        * helper search function used in LSH in order to find k nearest neighbors
        * and store them in a map. Used with continuous frehet
        */

        map<double, string> distance_map;
        int query_trick = 0;

        int current_bucket_num = h_fun->hashValue(key, table_size, &query_trick);
        cout << "bucket = " << current_bucket_num << endl;

        struct Bucket<K> *bucket = hash_table[current_bucket_num];
        if(bucket) {

            // traverse the list of Data at current bucket.
            typename std::list<BucketEntry<K> *>::iterator it;
            for (it = bucket->bucket_entries.begin(); it != bucket->bucket_entries.end(); ++it) {
            
                if(query_trick == (*it)->getQueryTrick()) {
                    // we need to initialise frechet structures here
                    // our points are represented in R
                    Point p1(1);
                    Points total_input_points(1);

                    vector<double> filtered_input_curve_vector = (*it)->getVector();

                    for(size_t i = 0; i < filtered_input_curve_vector.size(); i++) {
                        p1.set(0, filtered_input_curve_vector.at(i));
                        // create all the points from this curve with same dimension to p1
                        total_input_points.add(p1, p1.dimensions());
                    }

                    // do the same for queries
                    Points total_query_points(1);
                    Point p2(1);
                    for(size_t i = 0; i < query->getVector().size(); i++) {
                        p2.set(0, query->getVector().at(i));

                        total_query_points.add(p2, p2.dimensions());
                    } 
                
                    // let's init the curves c1 and c2 now
                    Curve input_curve(total_input_points, (*it)->get_name());
                    Curve query_curve(total_query_points, query->get_name());

                    // distance
                    Distance d = distance(input_curve, query_curve);
                    *duration += d.time_searches;

                    distance_map.insert(pair<double, string>(d.value, (*it)->get_name()));

                    p1.clear(); total_input_points.clear();
                    p2.clear(); total_query_points.clear();
                    input_curve.clear(); query_curve.clear();

                    // double euDist = euclidean_dist(q, (*it)->getVector());
                } else 
                    continue;           
            }

            // keep the best from the bucket.
            if(distance_map.size() > 0) {
                *nn_pair = pair<double,string>(distance_map.begin()->first, distance_map.begin()->second);
            
                cout << "mpike " << endl;
                cout << "\n " << nn_pair->first << " " << nn_pair->second << endl;
            }
            distance_map.clear();   
        }
    }

    // void find_NN_Continuous_Frehet(vector<double> query_vector, map<double, int> k_map, int k) {
    //     search_NN(query_vector, k_map, k);
    // }

    void PrintHT(){
        for (int i = 0; i < table_size; ++i) {
            cout<<"++++++ BUCKET ["<<i<<"] ++++++"<<endl;

            if(hash_table[i] != NULL)
                hash_table[i]->displayBucket();
        } 
    }

    ~HashTable() {
        // cout << "Destructing Hashtable..." << endl;
        Bucket<K> *temp = NULL;
        for(int i = 0; i < table_size; i++) {

            temp = hash_table[i];

            if(temp != NULL) {
                delete temp;

            }

            hash_table[i] = NULL;
        }
        delete[] hash_table;

        delete h_fun;
    }
};






















/*
#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <bits/stdc++.h>
#include <assert.h>
#include <ctime>
#include <map>

#include "../Curve.hpp"
#include "../Data/Data.hpp"
#include "./HashFunction.hpp"

using namespace std;


template <typename K>
struct BucketEntry {
    //Info of an entry inside a bucket.
    int query_trick;
    K entry;
    BucketEntry(K _entry, int id) {
        entry = _entry;
        query_trick = id;
    }
    ~BucketEntry() {
    }

    int getId() { return entry->id; }
};



template <typename K>
struct Bucket {
    int number_of_entries;
    int id;
    list<BucketEntry<K>*> bucket_entries; 

    Bucket(int i) {
        id = i;
        number_of_entries = 0;
    }

    ~Bucket() {
        // cout << "Destructing Bucket...";
        for(typename list<BucketEntry<K>*>::const_iterator it = bucket_entries.begin(); it != bucket_entries.end(); ++it)
        {
            delete *it;
        } 
        bucket_entries.clear();
    }

    void insert(K key, int _id){
        bucket_entries.push_back(new BucketEntry<K>(key, _id));
        // increase list items
        number_of_entries++;
    }

    void displayBucket(){
        // cout<<"BUCKET ["<<id<<"] :"<<endl;
        if(!bucket_entries.empty())
            for(auto &it : bucket_entries)
                it->key->PrintID();

    }
};



template <typename K>
class HashTable{   
    int buckets;    // Number of buckets .
    int table_size;   // Size of the Table.
    struct Bucket<K> **hash_table ;   // Buckets.
    HashFunction *h_fun;  //It's unique HashFunction.
    int ht_id;
    //Functionality.
    public:
        //Constructor.
        HashTable(int ht_size, int w, int k, int dim ,int _id) 
        :table_size(ht_size){
            this->ht_id = _id;
            hash_table = new Bucket<K>*[table_size]; // Allocate BUckets.
            for(int i = 0; i < table_size; i++)
                hash_table[i] = NULL; 
            h_fun = new HashFunction(w, k, dim);     // Generate HashFunction.
            cout<<"HashTable["<<ht_id<<"] is initialized."<<endl;
        }

        // Get Hash Function.
        HashFunction *getHashFunction() { return h_fun; }

        // Insert Key to the HT.
        void insertGridCurve(Grid_Curve *key) {
            int query_trick = 0;  // QueryTrick;
            unsigned int index = h_fun->hashValue(key->get_curve(), table_size, &query_trick);
            assert(index <= INT_MAX); //A little Check.
            if(hash_table[index] != NULL) {
                cout<<" Curve["<<key->get_id()<<"] assigned to Bucket["<<index<<"] "<<endl;
                hash_table[index]->insert(key, query_trick); }
            else {
                cout<<" Curve["<<key->get_id()<<"] is the 1st to be assigned to Bucket["<<index<<"] "<<endl;
                hash_table[index] = new Bucket<K>(index);
                hash_table[index]->insert(key, query_trick);
            }
        }
};

*/


