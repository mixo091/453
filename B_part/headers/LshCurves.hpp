#pragma once
#include <iostream>
#include <vector>
#include <limits>
#include <map>
#include <ctime>
#include <chrono>
#include <set>
#include <cmath>
#include <tuple>
#include <random>
#include "./Curve.hpp"
#include "./CurveR2.hpp"
#include "./GridCurve.hpp"
#include "./Point.hpp"
#include "./HashFunction.hpp"
#include "./DiscreteFrechet.hpp"
#define BUCKET_DIVIDER 16 


class _Bucket{
    public:
    int BucketId;
    vector<Grid_Curve*> BucketEntries; 
    _Bucket(int _id){
      this->BucketId = _id;
      this->BucketEntries.clear();
    }   
    ~_Bucket(){
      this->BucketEntries.clear();
    }
};

class _HashTable{
    public:
      int Id;
      int Size;
      HashFunction * hashFunction; 
      vector<Grid_Curve>* snapped_data;
      vector<_Bucket> Buckets;

      _HashTable(int size, int w, int k, int dim ,int _id,vector<Grid_Curve>* _snapped_data):Size(size),Id(_id){
            hashFunction = new HashFunction(w, k, dim);
            for(int i=0 ; i<Size; i++){
              _Bucket bucket(i);
              Buckets.push_back(bucket);
            } 
            snapped_data = _snapped_data;
            cout<<"HashTable["<<Id<<"] is initialized."<<endl;
      }
      void Print_snapped_data(){
        cout<<"Snapped Data"<<endl;
        for (int i =0 ; i < snapped_data->size(); i++){
            cout <<snapped_data->at(i).get_name()<<endl;
        }


      }
      HashFunction *GetHashFunction() { return hashFunction;}

      void HashCurve(Grid_Curve* curve){
        // Calculate Index. //
        int query_trick =0 ;
        int index = hashFunction->hashValue(curve->get_curve(), Size, &query_trick);
        assert(index <= INT_MAX);
        cout<<"Curve["<<curve->get_id()<<"] assigned to Bucket["<<index<<"] of HashTable["<<Id<<"]"<<endl;
        // Drop a Identifier Of The Curve To The Bucket. //
        Buckets.at(index).BucketEntries.push_back(&snapped_data->at(curve->get_id()));
      }

      void Print(){
        cout<<"HashTable ["<<Id<<"]"<<endl;
        for(int b=0 ; b<Size ; b++){
          cout<<" Bucket["<<b<<"]"<<endl;
          if(Buckets.at(b).BucketEntries.size()>0){
            cout<<"     [";
            for(int entry=0; entry < Buckets.at(b).BucketEntries.size() ; entry++){
              cout<<Buckets.at(b).BucketEntries.at(entry)->get_id()<<" , ";
              //cout<<Buckets.at(b).BucketEntries1.at(entry).get_id()<<" originalCurve["<<Buckets.at(b).BucketEntries1.at(entry).original_curve->get_id()<<"] , ";
            }
            cout<<"]"<<endl;
          }
        }



      }

      void FindNearestNeighbor(Grid_Curve* curve, map<double, int> &k_nearest_map, int k) 
      {
        map<double, int> distance_map;
        int query_trick = 0;
        curve->original_curve->print_r2_curve();
        int queryBucket = hashFunction->hashValue(curve->get_curve(), Size, &query_trick);
        //assert(index <= INT_MAX);
        cout<<"Query "<< curve->get_id() << " falls to Bucket "<<queryBucket<<endl;

        if(Buckets.at(queryBucket).BucketEntries.size()>0){
          cout <<" ok"<<endl;

          //Lest find the nearest neighbors in this bucket.
          for(int i = 0 ; i < Buckets.at(queryBucket).BucketEntries.size() ; i++){
            Buckets.at(queryBucket).BucketEntries.at(i)->original_curve->PrintID();
            curve->original_curve->PrintID();
            double DFT_distance = Discrete_Frechet(curve->original_curve->getVector(), Buckets.at(queryBucket).BucketEntries.at(i)->original_curve->getVector());
            cout <<"frfr"  <<DFT_distance<<endl;
            distance_map.insert(pair<double, int>(DFT_distance,Buckets.at(queryBucket).BucketEntries.at(i)->get_id()));
          }
          //Now keep the Best .//
          int item = 0 ;
            for (auto it = distance_map.cbegin(); it != distance_map.cend(); ++it) {
                k_nearest_map.insert(pair<double, int>((*it).first,(*it).second));

                // if k reached, stop adding new neigbor
                if(++item >= k) break;
            }




        }
        return;
      }





};



// Grid structure needed by Lsh_Curve for generating snapped_Vector from Curve. //
class Grid{
  public:
    int number_of_curves;
    int number_of_queries;
    int num_of_points;
    default_random_engine generator;
    uniform_real_distribution<double> distribution;
    pair<double,double> t;
    double delta;
    vector<Point> grid;
    double padding_value = numeric_limits<double>::max();
    public:
    int grid_id;
    vector<R2_Curve*>  input;
    vector<R2_Curve*>  queries;
    vector<Grid_Curve> snapped_data;
    vector<Grid_Curve> snapped_queries;
    //  Grid Constructor. //
    Grid(int _number_of_curves,int _num_of_points,double _delta ,int _number_of_queries , int _grid_id, R2_Curve* _input,R2_Curve* _query_data_set):
        number_of_curves(_number_of_curves), number_of_queries(_number_of_queries),num_of_points(_num_of_points),delta(_delta){
        // Every Grid parameritized by t at each Dimension. //
        grid_id = _grid_id;
        t.first = distribution(generator);
        t.second = distribution(generator);
        for(int i = 0 ; i < number_of_curves ; i++)
          {input.push_back(&_input[i]);}
        for(int j = 0 ; j < number_of_curves ; j++)
          {queries.push_back(&_query_data_set[j]);}
        cout<<" Grid["<<grid_id<<"] is initialized."<<endl;

    }


    Grid_Curve Curve_to_GridCurve(R2_Curve curve, bool query){
      // Gets a curve in R-2 space and returns it's snapped grid_curve. //
      vector<Point> curve_points = curve.getVector();
      vector<Point> snapped_points ;
      vector<double> result_vector ;
      // Snap each Dimension. //
      for(int i =0; i < num_of_points ; i++){
        Point p = curve_points.at(i);
        double x_i = p.get_x();
        double y_i = p.get_y();
        
        
          //cout<<"x_i"<<x_i<<"y_i"<<y_i<<endl;
        double snapped_x = floor((x_i - t.first)/this->delta + (double)1/2) * this->delta +t.first;
        double snapped_y =  floor((y_i - t.second)/this->delta + (double)1/2 )* this->delta +t.first;
        Point snapped_point(snapped_x ,snapped_y);

        // Check if it is a duplicate.If not add it. //
        bool duplicate = false;
        for(int i=0; i<snapped_points.size(); i++){
            if((snapped_points.at(i).get_x() == snapped_point.get_x()) && (snapped_points.at(i).get_y() == snapped_point.get_y())){
                duplicate = true;
                break;
            }
        }
        if(duplicate == false){
          snapped_points.push_back(snapped_point);}

      }     
      // Add snapped points to a 1D vector and then do the needed padding. //
      for(int i=0; i<snapped_points.size(); i++){
        result_vector.push_back(snapped_points.at(i).get_x());
        result_vector.push_back(snapped_points.at(i).get_y());     
      }
      for(int i = result_vector.size(); i<num_of_points*2; i++){
        result_vector.push_back(padding_value);
      }

      //for ( int k = 0 ; k < result_vector.size() ; k++){
            //cout<<result_vector.at(k)<<"    ";

     // }
      //cout<<endl;

      
      // Generate and Return Grid_Curve. // 

        if( query == true ){
          cout<<"dggr"<<endl;

          Grid_Curve snappedCurve1(curve.get_name(), curve.get_id(), queries.at(curve.get_id()), result_vector);
        return snappedCurve1;}
        Grid_Curve snappedCurve(curve.get_name(), curve.get_id(), input.at(curve.get_id()), result_vector);
        return snappedCurve;


    }

    // Snapp Input Curve. //
    void SnappInputData(R2_Curve input_curve){
        Grid_Curve snapped_curve = Curve_to_GridCurve(input_curve,false);
        this->snapped_data.push_back(snapped_curve);               
    }
    // Snap Query Curve. //
    void SnappQueryData(R2_Curve query_curve){
        query_curve.print_r2_curve();
        Grid_Curve snapped_curve = Curve_to_GridCurve(query_curve,true);
        this->snapped_queries.push_back(snapped_curve);            
    }
    // Get Snapped Input. //
    vector<Grid_Curve> GetSnappedInput(){return this->snapped_data;}
    // Get Snapped Queries. //
    vector<Grid_Curve> GetSnappedQueries(){return this->snapped_queries;}

    // Print Snapped Data. //
    void PrintSnappedDataSets(){
      cout<<" Snapped Data - Grid["<<this->grid_id<<"] "<<endl;
      cout<<"   [Input_Data_Set]:  "<<endl;
      for(int i = 0 ; i < snapped_data.size() ; i++){
          snapped_data.at(i).print_grid_curve();
          cout<<"  " <<snapped_data.at(i).original_curve->get_id()<<endl;
      }
      cout<<"   [Query_Data_Set]:  "<<endl;
      for(int i = 0 ; i < snapped_queries.size() ; i++){
          snapped_queries.at(i).print_grid_curve();
          cout<<"  " <<snapped_queries.at(i).original_curve->get_id()<<endl;
      }
    }

    // Print Original Curves. //
    void PrintOriginalData(){
      cout<<"Grid {"<<this->grid_id<<"}"<<endl;
        for(int i = 0 ; i < number_of_curves ; i++)
          {input.at(i)->print_r2_curve();}
    }

    // Destructor. //
    ~Grid(){
      // Destructor.
    }


};




// Structure For Curve Lsh using Discrete Fretech. //

class Curve_Lsh{
    int L;
    int w;
    int number_of_hashfunctions;
    int delta;
    int number_of_curves;
    int  number_of_queries;
    int hash_table_size;
    int num_of_points_per_curve;
    R2_Curve* input_data_set;
    R2_Curve* query_data_set;
    _HashTable** HashTables; 
    Grid** Grids ;

    public:
        Curve_Lsh(int L, int total_curves,int _number_of_queries, int number_of_points, int _k, int _w,double _delta ,R2_Curve *data_set ,R2_Curve *_query_data_set, const string &out_file)
            :L(L), hash_table_size(total_curves / BUCKET_DIVIDER), num_of_points_per_curve(number_of_points), number_of_hashfunctions(_k), 
            w(_w),delta(_delta),input_data_set(data_set),query_data_set(_query_data_set),number_of_curves(total_curves),number_of_queries(_number_of_queries){
      
              // Initialize L Grids. // 
              Grids = new Grid*[L] ;
              for(int i = 0 ; i < L ; i++)
                {Grids[i] = new Grid(total_curves, num_of_points_per_curve,delta ,number_of_queries ,i,input_data_set,query_data_set);}

              // Snapp Input Data_Set . // 
              for(int c = 0 ; c < number_of_curves ; c++){
                for(int i = 0; i < L ; i++){
                  Grids[i]->SnappInputData(input_data_set[c]);
                }
              }


              // Snapp Query Data_Set. // 
              for(int c = 0 ; c < number_of_queries ; c++){
                 for(int i = 0; i < L ; i++){
                  Grids[i]->SnappQueryData(query_data_set[c]);
                }
              }

              // Α Little Check of the SnappedData. //
              /***cout<<" ~~ Snapped Data ~~"<<endl;***/
              //Grids[0]->PrintSnappedDataSets(); 
              cout<<"Check[3] : Data Snapped Successfully."<<endl;
              // Initialize L HashTables. //
              HashTables = new _HashTable*[L];
              for(int i = 0; i < L ; i++) 
                {HashTables[i] = new _HashTable(hash_table_size, w, number_of_hashfunctions, num_of_points_per_curve,i,&Grids[i]->snapped_data);}
              cout<<"Check[4] : HashTables Initialized Successfully."<<endl;
              HashTables[0]->Print_snapped_data();
              // Hash Snapped Data. //
              for(int c = 0 ; c < number_of_curves ; c++){
                for(int i =0 ; i < L ; i++){
                    HashTables[i]->HashCurve(&(Grids[i]->snapped_data.at(c)));
                }
              }
              // Print Hash Tables. //
                PrintHashTables();


              // Open Output File. //
              ofstream output(out_file);
              if(!output.is_open()) {
                cout << "Error with output file\n";
                exit(-1);
              }


              
              // Lets Finally Search .//

              // create a map of <id, eu_dist> where id is int and eu_dist is double
              map<double, int> my_map;
              map<double, int> result_map;
              
              for(int query=0 ; query < number_of_queries ; query++){
               // auto start_lsh = high_resolution_clock::now();
                // Get every query curve and hash it. //
                for(int j=0 ; j < L ; j++){
                  HashTables[j]->FindNearestNeighbor( &(Grids[j]->snapped_queries.at(query)), my_map,  1);
                  // Add to Results. //
                  cout<<"Size "<<my_map.size()<<endl;
                  for (auto it = my_map.cbegin(); it != my_map.cend(); ++it) {
                      cout<<"ahytjyjtjuuj4j45u5j4j457j57j57j7555555jjjjjjjjjjjjjjjjjjjjjjjjjdd res "<< (*it).first<<"  "<<(*it).second  <<endl;
                      result_map.insert(pair<double, int>((*it).first, (*it).second));    
                  }
                  my_map.clear();
                
                }



                double min_true_distance = Discrete_Frechet(query_data_set[query].getVector(),input_data_set[0].getVector());
                string neighbhor_name = input_data_set[0].name;
                int neighbor_id = 0;
                for(int n = 0 ; n < number_of_curves ; n++){
                  if(min_true_distance > Discrete_Frechet(query_data_set[query].getVector(),input_data_set[n].getVector())){
                      min_true_distance = Discrete_Frechet(query_data_set[query].getVector(),input_data_set[n].getVector());
                      neighbhor_name = input_data_set[n].name;
                      neighbor_id = n;
                  }
                }

                output << "\nQuery: " <<" [ name:"<< this->query_data_set[query].get_name()<<" , id:"<< this->query_data_set[query].get_id()<<" ] ";
                output<< "\n Algorithm : Discrete Frechet";
                for(auto it = result_map.cbegin(); it != result_map.cend(); ++it) {
                  output<<"\nApproximate Nearest neighbor: "<< this->input_data_set[(*it).second].get_name();
                  output<<"\nApproximate Distance: "<< (*it).first<<endl;
                  output<<"\nTrue Nearest neighbor: "<< neighbhor_name <<endl;
                  output<<"\nTrue Distance: "<< min_true_distance<<endl;
                  break;
              
                }
                result_map.clear();
                my_map.clear();


              }

            output.close();

            }


    
        // Function to Print The HashTables. // 
        void PrintHashTables(){
          for(int i =0 ; i < L ; i++){
            HashTables[i]->Print();
          }
        }








        
      // Constructor for Clustering .//
      Curve_Lsh(int L, int total_curves, int number_of_points, int _k, int _w,double _delta ,R2_Curve *data_set , const string &out_file):
          L(L), hash_table_size(total_curves / BUCKET_DIVIDER), num_of_points_per_curve(number_of_points), number_of_hashfunctions(_k), 
          w(_w),delta(_delta),input_data_set(data_set){

              // Initialize L Grids. // 
              Grids = new Grid*[L] ;
              for(int i = 0 ; i < L ; i++)
                {Grids[i] = new Grid(total_curves, num_of_points_per_curve,delta ,number_of_queries ,i,input_data_set,query_data_set);}

              // Snapp Input Data_Set . // 
              for(int c = 0 ; c < number_of_curves ; c++){
                for(int i = 0; i < L ; i++){
                  Grids[i]->SnappInputData(input_data_set[c]);
                }
              }

              cout<<"ggegwr"<<endl;
              // Α Little Check of the SnappedData. //
              /***cout<<" ~~ Snapped Data ~~"<<endl;***/
              //Grids[0]->PrintSnappedDataSets(); 
              cout<<"Check[3] : Data Snapped Successfully."<<endl;
              // Initialize L HashTables. //
              HashTables = new _HashTable*[L];
              for(int i = 0; i < L ; i++) 
                {HashTables[i] = new _HashTable(hash_table_size, w, number_of_hashfunctions, num_of_points_per_curve,i,&Grids[i]->snapped_data);}
              cout<<"Check[4] : HashTables Initialized Successfully."<<endl;
              HashTables[0]->Print_snapped_data();
              // Hash Snapped Data. //
              for(int c = 0 ; c < number_of_curves ; c++){
                for(int i =0 ; i < L ; i++){
                    HashTables[i]->HashCurve(&(Grids[i]->snapped_data.at(c)));
                }
              }
              // Print Hash Tables. //
              PrintHashTables();
      



        }



             

  };

          

      

