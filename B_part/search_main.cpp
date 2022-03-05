#include <iostream>
#include <vector>
#include <string.h>
#include <sstream>
#include <fstream>
#include <string>
#include  "./headers/Functionality.hpp"
#include  "./headers/Curve.hpp"
#include  "./headers/CurveR2.hpp"
#include  "./headers/LshCurves.hpp"

#include "./headers/Data/Data.hpp"
#include "./headers/LSH/lsh.hpp"

using namespace std;
int main(int argc,char *argv[]){
    string input_file, output_file, query_file;
    int k,L,M,probes;
    int w = 3000;
    int N = 2;
    L = 5;
    k =5;
    string algorithm , metric;
    double delta = 1.0 ;
    //parse_arguments_search(argc , argv , &input_file , &query_file , &k , &L , &M , 
    //&probes , &output_file , &algorithm , &metric , &delta);


    // Files. //
    input_file = "./TestSets/FinalSets/nasd_input.csv";
    query_file = "./TestSets/FinalSets/nasd_query.csv";
    output_file = "./Results/result.txt";

    // <LSH or Hypercube or Frechet> .//
    string command = "Frechet";
    //<discrete or continuous | only for –algorithm Frechet>.//
    string _metric = "discrete";

    // --- Get the Info Of the Input File && Store it's Data. ---//
    int total_input_curves, number_of_points;
    calculate_file_info(&total_input_curves, &number_of_points, &input_file);
    cout<<"---- Check[1] ----"<<endl<<"  Input_File::"<<endl<<"  Number of curves :"<<total_input_curves<<endl<<"  Number of points : "<<number_of_points<<endl;
    Curve<double> dataset[total_input_curves];
    Data <double> DataSet[total_input_curves];
    store_curves(input_file, number_of_points, dataset); //Store As a Vector.
    store_curves(input_file, number_of_points, DataSet); //Store As a Vector.
    R2_Curve curve_dataset[total_input_curves];
    store_as_curves_r2(input_file, number_of_points,curve_dataset); //Store As a Curve.

    // --- Get the Info Of the Query File && Store it's Data. ---//
    int total_query_curves, number_of_points_q;
    calculate_file_info(&total_query_curves, &number_of_points_q, &query_file);
    cout<<"---- Check[2] ----"<<endl<<"  Query_File::"<<endl<<"  Number of curves :"<<total_query_curves<<endl<<"  Number of points : "<<number_of_points_q<<endl;
    Curve<double> query_dataset[total_query_curves];
    Data <double> QuerySet[total_input_curves];
    store_curves(query_file, number_of_points_q, query_dataset); //Store As a Vector.
    store_curves(query_file, number_of_points_q, QuerySet); //Store As a Vector.
    R2_Curve query_curve_dataset[total_query_curves];
    store_as_curves_r2(query_file, number_of_points_q,query_curve_dataset); //Store As a Curve.

    if( command.compare("LSH") == 0 ){
        // Represent Curves as Vectors. //
            //Initialize Plain LSH structure . //
            cout<<"Algorithm:Lsh"<<endl;
            Lsh<double> lsh = Lsh<double>(L, total_input_curves, number_of_points, k, w,DataSet);
            //lsh.ANN(QuerySet, total_query_curves, DataSet, total_input_curves, N, output_file, 300);
    }else if ( command.compare("Hypercube") == 0){
        // Represent Curves as Vectors. //
            //Initialize Plain HyperCube structure . //
            cout<<"Algorithm:Hypercube"<<endl;

    }else if ( command.compare("Frechet") == 0 ){
    // Represent Curves on the R^2 .//
        //--- Initialize Lsh for Curves. ---//
        Curve_Lsh Lsh_curves = Curve_Lsh(L, total_input_curves,total_query_curves,number_of_points, k, w, delta, curve_dataset,query_curve_dataset,output_file);
        //Lsh_curves.print_lsh_curves_info();
    }

    return 1;}