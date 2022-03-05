#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <bits/stdc++.h>
#include <assert.h>
#include <ctime>
#include <limits>
#include <map>
#include <string>
#include <random>
#include <chrono>	
#include <algorithm> 
#include <math.h>	
#include <fstream>

#include "../Curve.hpp"
#include "../Data/Data.hpp"
#include "../LSH/lsh.hpp"
#include "../hypercube/hypercube.hpp"


#define NONE -1
#define MIN_CHANGE 0
#define MAX_ITERATIONS 50
#define FIRST 0 
#define LAST MAX_ITERATIONS-1
#define MIDDLE (MAX_ITERATIONS - 1)/2
using namespace std;

template <typename K>
class Cluster
{   
    public:
    int number_of_items;
    int Cluster_Id;
    list<K*> Cluster_items; 
    vector<double> centroid;
    int CentroidId;

    // Constructor
    Cluster(int id){
        cout<<"New Cluster"<<endl;
        Cluster_Id = id;
        CentroidId = NONE ;
        number_of_items = 0;
    }

    void PrintCluster(){
        cout<<"CLUSTER["<<Cluster_Id<<"]:"<<endl;
        cout<<" INFO:"<<endl;
        cout<<"  CENTROID: "<<endl;
        cout<<"   ID:["<<CentroidId<<"]\n   VAVLUE:[";
        for(auto it = centroid.begin(); it < centroid.end(); it++) {
            cout << *it << " ";
        } 
        cout <<"]"<<endl;
        cout<<     "ITEMS:"<<number_of_items<<endl;
        cout<<"     ITEMS ID[";
        for (typename std::list<K*>::iterator it = Cluster_items.begin(); it != Cluster_items.end(); ++it){
             cout<<(*it)->id<<" ,";
        }
         
        cout <<"]"<<endl;
    }
    void InsertItem(Data<double> *Item) {
        Cluster_items.push_back(Item);
        number_of_items++;
    }
    void SetCenter(vector<double> newCenter,int id ){
        cout<<id<<endl;
        centroid = newCenter;
        CentroidId = id;
    }
    vector<double> GetCenter(){
        return centroid;
    }
    int GetId(){
        return  Cluster_Id;
    }    

    void ClearCluster(){
        //Clear assigned Items.
        Cluster_items.clear();
        number_of_items = 0;
    }
    int CalculateCentroid(){
        int ClusterChange = 0;
        if(Cluster_items.size()>0){
        vector<double> old_centroid = centroid;
        int dimension = centroid.size();
        for(int i =0 ; i < dimension ;i++){
            centroid.at(i) = 0.0;
            double sum  = 0.0;
            for (typename std::list<K*>::iterator it = Cluster_items.begin(); it != Cluster_items.end(); ++it){
                 sum = sum + (*it)->v.at(i);
            }
            centroid.at(i) = sum/number_of_items;
            ClusterChange = ClusterChange + abs(int(old_centroid.at(i)) - (sum/number_of_items));
        }
        old_centroid.clear();

        }
        return ClusterChange;
    }

};

//---Implementing Clustering Structure---//

template <typename K, typename T>
class Clustering {

    Cluster<K> **Clusters;
    int ClustersNum;
    Data<T> *DataSet;
    int TotalDataItems ; 

    // Lsh data needed
    int lsh_w;
    int lsh_k;
    int lsh_L;
    int lsh_dim;
    Lsh<T> *cluster_lsh;

    // Hypercube data
    int M_cube;
    int probes_cube;
    int w_cube;
    int dimension_cube;
    hypercube<T> *cluster_hcube;

    // R^2 Frechet Clustering .
    
    public:
        // Clustering(int k , Curve<T> *dataSet ,int TotalVectors, int L, int k_lsh, int dimension){
        //     DataSet = dataSet;
        //     TotalDataItems = TotalVectors;
        //     ClustersNum = k;
        //     Clusters = new Cluster<K>*[ClustersNum];
        //     for(int i = 0 ; i<k ; i++ )
        //         Clusters[i] = new Cluster<Curve<T>>(i);

        //     // initialise lsh structure
        //     lsh_L = L; lsh_k = k;
        //     int w = 1000;
        //     cluster_lsh = new Lsh<K>(L, TotalVectors, dimension, k, w , DataSet);

        // }

        Clustering(int k_clusters, Data<T> *dataSet, int TotalVectors, 
                    int w_lsh, int k_lsh, int L_lsh, int dimension_lsh, 
                    int M, int probes, int w_cube, int cube_dim)
        : // initialisation list
        Clusters(new Cluster<K>*[k_clusters]), ClustersNum(k_clusters), 
        DataSet(dataSet), TotalDataItems(TotalVectors),
        lsh_w(w_lsh), lsh_k(k_lsh), lsh_L(L_lsh), lsh_dim(dimension_lsh), cluster_lsh(NULL),
        M_cube(M), probes_cube(probes), w_cube(w_cube), dimension_cube(cube_dim), cluster_hcube(NULL)
        // cluster_lsh(new Lsh<T>(lsh_L, TotalVectors, dimension_lsh, lsh_k, lsh_w, DataSet))
        {

            cout<<"Total Items : "<<TotalDataItems<<endl;
            // constructor used for lsh initialisation
            for(int i = 0 ; i < ClustersNum ; i++)
            {   cout<<"A cluster is being Constructed"<<endl;
                Clusters[i] = new Cluster<K>(i);
            }
//PrintData();
                
        }

        void PrintData(){
            for (int i = 0; i < TotalDataItems; i++) 
                DataSet[i].printVector();
        }


        void PrintClusters(){
            for(int i = 0 ; i < ClustersNum ; i++)
                Clusters[i]->PrintCluster();
        }

        void ClearClusters(){
            for(int i = 0 ; i < ClustersNum ; i++)
                Clusters[i]->ClearCluster();
        }

        void GenerateCentroids(){
            //=== Generate K random Centroids from the dataset.
            srand(time(NULL));
            for (int i = 0; i < ClustersNum; i++){
                int randomCentroid = rand() % (TotalDataItems- 1) + 0;
                Clusters[i]->SetCenter(DataSet[randomCentroid].v,DataSet[randomCentroid].id);
            }
        }


        void K_meansPP()
        {    //Choose first Centroid Randomly.
            default_random_engine RandomEn(chrono::system_clock::now().time_since_epoch().count());
            int randomCentroid = rand()%TotalDataItems;
            cout<<"randomCentroid "<<randomCentroid<<endl;
            Clusters[0]->SetCenter(DataSet[randomCentroid].v,DataSet[randomCentroid].id);
            Clusters[0]->PrintCluster();
            for (int i = 1; i < ClustersNum; i++)
            {
                //cout<<"FOR CENTROID:"<<i<<endl;
                double maxDist = 0;
                vector<double> P(TotalDataItems + 1); 

                for (int j = 0 ; j<TotalDataItems ;j++)
                {   
                    double Di = 0;
                    //MinimumDistFromCentroids
                    vector<double> Distances;
                    for(int c=0; c<ClustersNum; c++){
                        if(Clusters[c]->CentroidId != NONE){
                            double dist = euclidean_dist(DataSet[j].v, DataSet[Clusters[c]->CentroidId].v);
                            Distances.push_back(dist);
                        }
                    }
                    if(Distances.size()!=0){
                        Di =  *min_element(begin(Distances), end(Distances));
                    }
                    //===============================//
                   //cout<<"Di"<<j<<" "<<Di<<endl;

                    if (Di > maxDist)
                        maxDist = Di;
                    
                }
                //cout<<"maxDist: "<<maxDist<<endl;
                P[0] = 0;
                for (int j = 1; j <TotalDataItems; j++)
                {

                    double Di = 0;
                    //MinimumDistFromCentroids
                    vector<double> Distances;
                    for(int c=0; c<ClustersNum; c++){
                        if(Clusters[c]->CentroidId != NONE){
                            double dist = euclidean_dist(DataSet[j].v, DataSet[Clusters[c]->CentroidId].v);
                            Distances.push_back(dist);
                        }
                    }
                    if(Distances.size()!=0){
                        Di =  *min_element(begin(Distances), end(Distances));
                    }
                    P[j] = pow(Di / maxDist, 2) + P[j - 1];
                }

                uniform_real_distribution<double> unif(0.0, P[TotalDataItems-1]);
                double x = unif(RandomEn);
                //cout<<"x:"<<x<<endl;

                int nextCentroid = 0;
                for (int p = 1; p < TotalDataItems; p++)
                {
                    if (x <= P[p])
                    {       
                        nextCentroid = p;
                        break;
                    }
                }
                Clusters[i]->SetCenter(DataSet[nextCentroid].v,DataSet[nextCentroid].id);
            }
       }

        double minDistBetween2centers() {
            double minDist = std::numeric_limits<double>::max();
            for(int i =0 ; i<ClustersNum ; i++){
                for(int j = i+1 ; j < ClustersNum ; j++ ){
                    double dist = euclidean_dist(DataSet[Clusters[i]->CentroidId].v, DataSet[Clusters[j]->CentroidId].v);
                    if(dist<minDist){
                        minDist = dist;
                    }
                }
            }
            return minDist ;
        }

        //=============  Loyds Clustering =============================//
        void Loyds_Clustering(const string &out_file , bool complete){
            ofstream output(out_file);
            if(!output.is_open()) {
                cout << "Error with output file\n";
                exit(-1);
            }
            int ClusterChange ;
            int TotalClusteringChange =  MIN_CHANGE + 1;

            cout<<" === Kmeans++ - Initialize Centroids === // "<<endl;
            K_meansPP();
            PrintClusters();
            int iteration = 0;
            auto start = high_resolution_clock::now();
            cout<<"ok"<<endl;
            while (iteration < MAX_ITERATIONS && TotalClusteringChange > MIN_CHANGE)
            {   
                TotalClusteringChange = 0;
                if((iteration == FIRST) || (iteration == MIDDLE) || (iteration==LAST))
                    cout<<"ITERATION "<<iteration<<endl;
                // === Clear Clusters === // 
                ClearClusters();
                //cout<<"TOTAL:"<<TotalDataItems<<endl;
                // === Reassign the DataSet to the Clusters. === //
                for (int i = 0; i < TotalDataItems; i++)
                {
                   //=== Calculate Distance from Each Cluster (Centroid). === //
                    map<double, int> L2fromClusters; 
                    for(int j = 0; j < ClustersNum; j++)
                    {   
                        double euDist = euclidean_dist(DataSet[i].v,Clusters[j]->centroid);
                        L2fromClusters.insert(pair<double, int>(euDist, j));
                    }
                    //Assign to Apropriate Cluster Based on L2 ===//
                    auto it = L2fromClusters.cbegin();
                    //cout<<(*it).second<<endl;
                    int ClusterToBeAssigned = (*it).second ;
                     if((iteration == FIRST) || (iteration==LAST )){
                         if(iteration==FIRST){
                            // FirstIteration<<"ITEM ID : "<<DataSet[i].id<<" --- Assigned to CLUSTER --------- "<<ClusterToBeAssigned<<endl;
                         }else{
                             //LastIteration<<"ITEM ID : "<<DataSet[i].id<<" --- Assigned to CLUSTER --------- "<<ClusterToBeAssigned<<endl;
                         }

                     }

                    Clusters[ClusterToBeAssigned]->InsertItem(&DataSet[i]);
                    //Clear the map to use it foe the next point 
                    L2fromClusters.clear();
                }
                // === Recalculate Centroids Based  on the Reassignment. ===//
                for(int i =0 ; i< ClustersNum; i++){
                    ClusterChange =0 ;
                    ClusterChange = Clusters[i]->CalculateCentroid();
                    TotalClusteringChange += ClusterChange;
                }
                //cout<<"ClusteringChange: "<<TotalClusteringChange<<endl;
                iteration++;
                }
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            //output..
            if(complete){
                output<<"Algorithm:Classic\n";
                for(int c = 0 ; c < ClustersNum ; c++){
                    output<<"Cluster-"<<c;
                    output<<" {";  
                    output<<",centroid:" ;
                    for(auto it =Clusters[c]->centroid.begin(); it < Clusters[c]->centroid.end(); it++) {
                         output << *it << " ";
                    } 
                    output<<", "<<endl;
                    //output<<Clusters[c]->id<<" ,";
                    for (std::list<Data<double>*>::iterator it = Clusters[c]->Cluster_items.begin(); it !=  Clusters[c]->Cluster_items.end(); ++it){
                        output<<(*it)->id<<" ,";
                    }
                    output<< "}"<<endl;
                }
                output<<"ClusteringTime : "<<duration.count()<<"microseconds"<<endl;
                
            }else{
                output<<"Algorithm:Classic\n";
                for(int c = 0 ; c < ClustersNum ; c++){
                    output<<"Cluster-"<<c;
                    output<<" { size:"<<Clusters[c]->number_of_items;
                    output<<",centroid:" ;
                    for(auto it =Clusters[c]->centroid.begin(); it < Clusters[c]->centroid.end(); it++) {
                         output << *it << " ";
                    } 
                    output<< "}"<<endl;
                }
                output<<"ClusteringTime : "<<duration.count()<<"microseconds"<<endl;
            }
            PrintClusters();
        }


        void Lsh_Clustering(const string &out_file, bool complete) {
            int TotalDataItems90percent = int(TotalDataItems*90/100);
            int assigned = 0;  
            cout<<"Lsh is running...Get Results at :"<<out_file<<endl;
            //Initialize Centroids with Kmeans++ Initialazation.
            K_meansPP();
            double initialRadius =  minDistBetween2centers()/2  ; // minDist Between2centers.
            double R = initialRadius;
            
            ofstream output(out_file);
            if(!output.is_open()) {
                cout << "Error with output file\n";
                exit(-1);
            }
            output.clear();

            //Lsh Structure with the Dataset.
            vector<int> idsInQueryBall;

            // initialize lsh
            cluster_lsh = new Lsh<T>(lsh_L, TotalDataItems, lsh_dim, lsh_k, lsh_w, DataSet);

            map<int, int> AssignmentsInRadius;

            auto start = high_resolution_clock::now();
            int iteration = 0;
            while (iteration < MAX_ITERATIONS  && assigned < TotalDataItems90percent)
            {   
                cout<<" --- ITERATION --- "<<iteration<<endl;

                //------------------------------------------------------------------//
                for(int i=0 ; i < ClustersNum ;i++){
                    if(iteration == 0 ){
                        cout<<"--- Cluster "<<i<<" --- intitial Centroid "<<Clusters[i]->CentroidId<<endl;
                    }
                    cout<<"--- Query Cluster ---"<<i<<" at Radius ---"<<R<< endl;

                    //what items the Cluster go in this Radius R //
                    idsInQueryBall=cluster_lsh->ReverseAssignment(Clusters[i]->centroid,Clusters[i]->Cluster_Id,R);

                    // if this ball gets at least one point
                    if (idsInQueryBall.size() > 0 ){
                        cout<<"------ In Range of Cluster ---------- "<<i<<"---- iteration "<<iteration<<endl;
                        for(auto it = idsInQueryBall.begin(); it < idsInQueryBall.end(); it++) {
                            
                            if(it !=  idsInQueryBall.end()){
                            
                                int id = *it;
                                // if id exists
                                if(AssignmentsInRadius.count(id)){
                                    map<int, int>::iterator it1 = AssignmentsInRadius.find(id);
                                    double dist1=0.0;
                                    double dist2=0.0;
                                    //cout<<(*it1).second<<endl;
                                    dist1 = euclidean_dist(DataSet[id-1].v,Clusters[i]->centroid);
                                    dist2 = euclidean_dist(DataSet[id-1].v,Clusters[(*it1).second]->centroid);
                                    if(dist1>dist2) {
                                        AssignmentsInRadius.erase(id);
                                        AssignmentsInRadius.insert(pair<int, int>(id, i));     
                                    }
                                    
                                }else{
                                    AssignmentsInRadius.insert(pair<int, int>(id, i)); 
                                }
                            }                  
                        }
                    }

                idsInQueryBall.clear();
            }
 
                map<int, int>::iterator it2 = AssignmentsInRadius.begin();
                while (it2 != AssignmentsInRadius.end())
                {       
                    int cluster_num = (*it2).second;
                    Clusters[cluster_num]->InsertItem(&DataSet[(*it2).first]);
                    assigned++;

                    it2++;
                }
                AssignmentsInRadius.clear();

                iteration++;
                R=R*2;
            }

            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            cout<<"assigned"<<assigned<<endl;
            if(complete){
                output<<"Algorithm:Range Search LSH\n";
                for(int c = 0 ; c < ClustersNum ; c++){
                    output<<"Cluster-"<<c;
                    output<<" {";  
                    output<<",centroid:" ;
                    for(auto it =Clusters[c]->centroid.begin(); it < Clusters[c]->centroid.end(); it++) {
                         output << *it << " ";
                    } 
                    output<<", "<<endl;
                    for (std::list<Data<double>*>::iterator it = Clusters[c]->Cluster_items.begin(); it !=  Clusters[c]->Cluster_items.end(); ++it){
                        output<<(*it)->id<<" ,";
                    }
                    output<< "}"<<endl;
                    
                }
                output<<"ClusteringTime : "<<duration.count()<<"microseconds"<<endl;
                
            }else if(complete == false){
                output<<"Algorithm:Range Search LSH\n";
                for(int c = 0 ; c < ClustersNum ; c++){
                    output<<"Cluster-"<<c;
                    output<<" { size:"<<Clusters[c]->number_of_items;
                    output<<",centroid:" ;
                    for(auto it =Clusters[c]->centroid.begin(); it < Clusters[c]->centroid.end(); it++) {
                         output << *it << " ";
                    } 
                }
                 output<<"ClusteringTime : "<<duration.count()<<"microseconds"<<endl;
            }
            PrintClusters();
            return;
        }

    
    // ========== Hypercube Clustering =========== //
    void hypercube_clustering(int dimension, const string &out_file, bool complete) {
        // int w = 1000;  
        // int M = 500;
        // int k = 13 ;
        // int probes = 50;
        // int TotalClusteringChange =  MIN_CHANGE + 1;
        int TotalDataItems90percent = int(TotalDataItems*90/100);
        int assigned = 0;  

        //Initialize Centroids with Kmeans++ Initialazation.
        K_meansPP();
        double initialRadius =  minDistBetween2centers()/2  ; // minDist Between2centers. 
        double R = initialRadius;
        
        ofstream output(out_file);
        if(!output.is_open()) {
            cout << "Error with output file\n";
            exit(-1);
        }
        output.clear();

        // this is the result vectord of id;s
        vector<int> idsInQueryBall;

        // initialisation of hypercube
        hypercube<double> cube = hypercube<double>(probes_cube, M_cube, w_cube, dimension_cube, 
                                                        lsh_dim, TotalDataItems, DataSet);
        
        // temp vector to store the assigned cluster for each data point
        map<int, int> AssignmentsInRadius;

        auto start = high_resolution_clock::now();
        
        int iteration = 0;
        while (iteration < MAX_ITERATIONS  && assigned < TotalDataItems90percent)
        {   
            cout<<" --- ITERATION --- "<<iteration<<endl;
            for(int i=0 ; i < ClustersNum ;i++){
                if(iteration == 0 ){
                    cout<<"--- Cluster "<<i<<" --- intitial Centroid "<<Clusters[i]->CentroidId<<endl;
                }
                cout<<"--- Query Cluster ---"<<i<<" at Radius ---"<<R<< endl;

                // what items the Cluster go in this Radius R 
                idsInQueryBall=cube.ReverseAssignment(Clusters[i]->centroid,Clusters[i]->Cluster_Id,R, probes_cube, M_cube);

                if (idsInQueryBall.size() > 0 ) {
                    cout<<"------ In Range of Cluster ---------- "<<i<<"---- iteration "<<iteration<<endl;
                    
                    for(auto it = idsInQueryBall.begin(); it < idsInQueryBall.end(); it++) {
                        if(it !=  idsInQueryBall.end()) {
                            
                            int id = *it;
                            if(AssignmentsInRadius.count(id)) {
                                //cout<< id << "is present! "<<endl;
                                map<int, int>::iterator it1 = AssignmentsInRadius.find(id);

                                double dist1=0.0;
                                double dist2=0.0;
                                //cout<<(*it1).second<<endl;
                                dist1 = euclidean_dist(DataSet[id-1].v,Clusters[i]->centroid);
                                dist2 = euclidean_dist(DataSet[id-1].v,Clusters[(*it1).second]->centroid);
                                //cout<<"dist1:"<<dist1<<" dist2:"<<dist2<<endl;
                                if(dist1>dist2){
                                    AssignmentsInRadius.erase(id);
                                    AssignmentsInRadius.insert(pair<int, int>(id, i));     
                                }
                                
                            }else{
                                AssignmentsInRadius.insert(pair<int, int>(id, i)); 
                            }
                        }                  
                    }
                }
                idsInQueryBall.clear();
            }
 
            map<int, int>::iterator it2 = AssignmentsInRadius.begin();
            while (it2 != AssignmentsInRadius.end())
            {       //cout<<"Iteration_ "<<iteration<<" --- Assignment in Radius :"<<endl;
                // int item_id = (*it2).first;s
                int cluster_num = (*it2).second;
                        //cout<<"Item_id_"<<item_id <<" && cluster_num_"<< cluster_num<<endl; 
                Clusters[cluster_num]->InsertItem(&DataSet[(*it2).first]);
                assigned++;

                it2++;
            }
            AssignmentsInRadius.clear();

            iteration++;
            R=R*2;
                
            }
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            cout<<"assigned"<<assigned<<endl;
        
            if(complete){
                output<<"Algorithm:Range Search Hypercube\n";
                for(int c = 0 ; c < ClustersNum ; c++){
                    output<<"Cluster-"<<c;
                    output<<" {";  
                    output<<",centroid:" ;
                    for(auto it =Clusters[c]->centroid.begin(); it < Clusters[c]->centroid.end(); it++) {
                         output << *it << " ";
                    } 
                    output<<", "<<endl;
                    for (std::list<Data<double>*>::iterator it = Clusters[c]->Cluster_items.begin(); it !=  Clusters[c]->Cluster_items.end(); ++it){
                        output<<(*it)->id<<" ,";
                    }
                    output<< "}"<<endl;
                    
                }
                output<<"ClusteringTime : "<<duration.count()<<"microseconds"<<endl;
                
            }else if(complete == false){
                output<<"Algorithm:Range Search Hypercube\n";
                for(int c = 0 ; c < ClustersNum ; c++){
                    output<<"Cluster-"<<c;
                    output<<" { size:"<<Clusters[c]->number_of_items;
                    output<<",centroid:" ;
                    for(auto it =Clusters[c]->centroid.begin(); it < Clusters[c]->centroid.end(); it++) {
                         output << *it << " ";
                    } 
                    output<< "}"<<endl;
                }
                 output<<"ClusteringTime : "<<duration.count()<<"microseconds"<<endl;
            } 
            PrintClusters();
        
        
        return;
    }



};