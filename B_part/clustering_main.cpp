#include <iostream>
#include <vector>
#include <string.h>
#include <sstream>
#include <fstream>
#include <string>
// Headers //
#include  "./headers/Functionality.hpp"
// #include  "./headers/Curve.hpp"
#include  "./headers/CurveR2.hpp"
// #include  "./headers/LshCurves.hpp"
//#include "./headers/Clustering/ClusteringCurves.hpp"
#include "./headers/Clustering/Clustering.hpp"
#include "./headers/Clustering/MeanCurve.hpp"
//#include "./headers/Clustering/ClusteringUtilities.hpp"
#include "./headers/Data/Data.hpp"

using namespace std;

// Clustering Main .//
/*   Usage of : ./cluster -i <input_file> -c <configuration_file> -o <output_file> -update <Mean Frecher or 
Mean Vector> -assignment <Classic or LSH or Hypercube or LSH_Frechet>   */
int main(int argc,char *argv[]){
    // Variables used to store the args for clustering execution
    string input_file, output_file, config_file;
    string update, assignment, complete;
    bool printComplete = false;
    
    // default parameters for LSH and Cube
    int w_lsh = 1000; int w_cube = 1000;
    int k_lsh = 5, L_lsh = 5;
    // total points and dimension
    int total_input_curves, number_of_points;

    // configuration file variables
    int number_of_clusters = 4;
    int number_of_hash_tables = 0;
    int number_of_vector_hash_functions = 0; 
    int max_number_M_hypercube = 0 ;
    int number_of_hypercube_dimensions = 0;
    int number_of_probes = 0;
    input_file = "./TestSets/FinalSets/nasd_input.csv";
    output_file =  "./Results/result.txt";
    cluster_parse_args(argc, argv, 
        &input_file, &config_file, &output_file, 
        &update, &assignment, &complete); 
    printComplete = (complete == "1" ? true : false);

    cout<<"Args Info [ "<<endl;
    cout<<"    input_file  : "<<input_file<<endl;
    cout<<"    output_file : "<<output_file<<endl;
    cout<<"    config_file  : "<<config_file<<endl;
    cout<<"    assignment : "<<assignment<<endl;
    cout<<"    update : "<<update<<endl;
    cout<<"    complete : " <<printComplete<<"   ]"<<endl;
    
    // Get File Info  && Store Data.//
    calculate_file_info(&total_input_curves, &number_of_points, &input_file);
    cout<<"File Info { Total Curves : "<<total_input_curves<<" ,  Points : " <<number_of_points<<" }"<<endl;
    Data<double> dataset[total_input_curves];
    store_curves(input_file, number_of_points, dataset);
    R2_Curve curve_dataset[total_input_curves];
    store_as_curves_r2(input_file, number_of_points,curve_dataset); //Store As a Curve.

    // Init Strucutres. //
    Clustering<Data<double>, double> clustering(number_of_clusters, dataset, 
        total_input_curves, w_lsh, k_lsh, L_lsh, number_of_vector_hash_functions,
        max_number_M_hypercube, number_of_probes,w_cube, number_of_hypercube_dimensions);


    /***            Run Algorithms          ***/

    if ( assignment == "LSH"  &&  update == "MeanVector" ){

        clustering.Lsh_Clustering(output_file, printComplete);

    }else if ( assignment == "Classic"  &&  update == "MeanVector" ){

        clustering.Loyds_Clustering(output_file, printComplete);

    }else if (  assignment == "Hypercube"  &&  update == "MeanVector" ){

        clustering.hypercube_clustering(number_of_hypercube_dimensions, output_file, printComplete);

    }else if ( assignment == "LSH" && update == "MeanFrechet"){

            //Not done

    }else if ( assignment == "Classic" && update == "MeanFrechet"){

       // Clustering_Curves<R2_Curve> clustering_curves(number_of_clusters, total_input_curves, number_of_points , 5 , curve_dataset, output_file);
       // clustering_curves.Loyds_Clustering_Curves(output_file, printComplete);

    }
    
return 0 ;
}