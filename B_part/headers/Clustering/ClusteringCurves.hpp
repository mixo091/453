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
#include "../CurveR2.hpp"
#include "../Point.hpp"
#include "../DiscreteFrechet.hpp"
#include "MeanCurve.hpp"


#define NONE -1
#define MIN_CHANGE 0
#define MAX_ITERATIONS 25
#define FIRST 0 
#define LAST MAX_ITERATIONS-1
#define MIDDLE (MAX_ITERATIONS - 1)/2

// Cluster Structure //

template <typename K>
class Curve_Cluster
{   
    public:
    int number_of_items;
    int Cluster_Id;
    list<K*> Cluster_items; 
    vector<Point> centroid;
    int CentroidId;

    // Constructor
    Curve_Cluster(int id){
        Cluster_Id = id;
        CentroidId = NONE ;
        number_of_items = 0;
    }

    void PrintCluster(){
        cout<<"CLUSTER["<<Cluster_Id<<"]:"<<endl;
        cout<<" INFO:"<<endl;
        cout<<"  CENTROID: "<<endl;
        cout<<"   ID:["<<CentroidId<<"]\n   VAVLUE:[";
            for( int i =0 ; i < centroid.size() ; i++){
                centroid.at(i).print_point();
                cout<<" , ";            }
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
    void SetCenter(vector<Point> newCenter,int id ){
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




/*

        void Calculate_new_center(){
            int ClusterChange = 0;
            if(Cluster_items.size()>2){
                vector<double> old_centroid = centroid;
                int number_of_items = Cluster_items.size();
                vector<Point> new_center;
                vector<R2_Curve*> items;
                for( int i = 0 ; i <number_of_items ; i++){
                    items.push_back(Cluster_items[i]);
                }
                Tree * tree = new Tree(items);

                new_center = 
                
                
                
                
            }


            return;
        }

        */

};

template <typename K>
class Clustering_Curves{

    // Clustering Info //
    Curve_Cluster<K> **Clusters;
    int ClustersNum;
    int TotalDataItems; 
    int L_grid;
    R2_Curve* DataSet;
    double delta = 4.0;
    string out_file ;




    public:

        // Construct Structure. //
        Clustering_Curves(int K_Cluster , int Total_Curves , int Num_Of_Points , int L_Grid , R2_Curve* CurveDataSet,string  output_file)
        : Clusters(new Curve_Cluster<K>*[K_Cluster]), ClustersNum(K_Cluster), DataSet(CurveDataSet), TotalDataItems(Total_Curves),L_grid(L_Grid){
            // Construct Clusters .//
            out_file = output_file;
            for(int i = 0 ; i < K_Cluster ; i++)
                Clusters[i] = new Curve_Cluster<K>(i);
        }

        // Print Input Data Set //
        void PrintInputDataSet(){
            cout<<"\nClustering Input Curves {";
            for(int i = 0; i < TotalDataItems; i++) {
               cout<<DataSet[i].get_name()<<" , ";}
            cout<<" }"<<endl;
        }
        

        void KmeanPP(){
            // Choose First Centroid Randomly .//
            default_random_engine RandomEn(chrono::system_clock::now().time_since_epoch().count());
            int randomCentroid = rand()%TotalDataItems;
            Clusters[0]->SetCenter(DataSet[randomCentroid].v,DataSet[randomCentroid].id);
            Clusters[0]->PrintCluster();
            for (int i = 1; i < ClustersNum; i++)
            {
        
                double maxDist = 0;
                vector<double> P(TotalDataItems + 1); 

                for (int j = 0 ; j<TotalDataItems ;j++)
                {   
                    double Di = 0;
                    //MinimumDistFromCentroids
                    vector<double> Distances;
                    for(int c=0; c<ClustersNum; c++){
                        if(Clusters[c]->CentroidId != NONE){
                            double dist = Discrete_Frechet(DataSet[j].v, DataSet[Clusters[c]->CentroidId].v);
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
                            double dist = Discrete_Frechet(DataSet[j].v, DataSet[Clusters[c]->CentroidId].v);
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

        void PrintClusters(){
            for ( int i =0 ; i < ClustersNum ; i++){
                Clusters[i]->PrintCluster();
            }
        }

        //=============  Loyds Clustering for Curves =============================//
        void Loyds_Clustering_Curves(const string &out_file, bool complete){
            ofstream output(out_file);
            if(!output.is_open()) {
                cout << "Error with output file\n";
                exit(-1);
            }
            int ClusterChange ;
            int TotalClusteringChange =  MIN_CHANGE + 1;

            // === Kmeans++ - Initialize Centroids === //
            K_meansPP();
            PrintClusters();
            int iteration = 0;
            auto start = high_resolution_clock::now();
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
                    map<double, int> DistancefromClusters; 
                    for(int j = 0; j < ClustersNum; j++)
                    {   
                        double Discrete_Frechet= Discrete_Frechet(DataSet[i].v,Clusters[j]->centroid);
                        DistancefromClusters.insert(pair<double, int>(euDist, j));
                    }
                    //Assign to Apropriate Cluster Based on L2 ===//
                    auto it = DistancefromClusters.cbegin();
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
                    DistancefromClusters.clear();
                }
                // === Recalculate Centroids Based  on the Reassignment. ===//
                for(int i =0 ; i< ClustersNum; i++){
                    ClusterChange =0 ;
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

};


