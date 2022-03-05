#include <iostream>
#include <vector>
#include <string.h>
#include <sstream>
#include <fstream>
#include <string>

#include  "./headers/Functionality.hpp"
// #include  "./headers/Curve.hpp"
// #include  "./headers/CurveR2.hpp"
#include  "./headers/LshCurves.hpp"
#include "./headers/LSH/lsh_CFrechet.hpp"
#include "./headers/LSH/lsh.hpp"
#include "./headers/hypercube/hypercube.hpp"
#include "./headers/Data/Data.hpp"
#include "./headers/CurveR2.hpp"

#include "./headers/frechet/frechet.hpp"
#include "./headers/frechet/curve.hpp"
#include "./headers/frechet/types.hpp"
#include "./headers/frechet/point.hpp"

using namespace std;

using namespace Frechet;
using namespace Continuous;

int main(int argc,char *argv[]) {
    string input_file, output_file, query_file;
    int k = 5,L = 1, M = 10, probes = 4;
    int w = 1000;
    int total_input_curves = 0, number_of_points = 0;
    int total_query_curves = 0, query_points = 0;

    string algorithm , metric;
    double delta = 0.1;
    parse_arguments_nn(argc , argv , &input_file , &query_file , &k , &L , &M , 
                         &probes , &output_file , &algorithm , &metric , &delta);

    // Get info of input file and store it's data
    calculate_file_info(&total_input_curves, &number_of_points, &input_file);
    // get info of query file
    calculate_file_info(&total_query_curves, &query_points, &query_file);

    Data<double> dataset[total_input_curves];
    Data<double> queries[total_query_curves];

    bool stop;
    do {

        if(algorithm.compare("Frechet") == 0 ) 
        {
            if(metric.compare("continuous") == 0)
            // filtering
            {
                store_curves(input_file, number_of_points, dataset);
                // Create query dataset
                store_curves(query_file, query_points, queries);

                Data<double> filtered_curves[total_input_curves];
                filtering_curves(total_input_curves, number_of_points, dataset, filtered_curves);

                // store in lsh table
                lsh_CFrehet<double> *_lsh_frehet = new lsh_CFrehet<double>(delta,
                                    number_of_points, total_input_curves,
                                    L, k, w, filtered_curves);
                // filtering query curves
                Data<double> filtered_queries[total_query_curves];
                filtering_curves(total_query_curves, query_points, queries, filtered_queries);

                // continuous frechet    
                _lsh_frehet->perforfm_continuous_frechet(output_file, total_query_curves, query_points, filtered_queries, filtered_curves);

                delete _lsh_frehet;
    
            } else if(metric.compare("discrete") == 0) {

                R2_Curve curve_dataset[total_input_curves];
                store_as_curves_r2(input_file, number_of_points, curve_dataset); //Store As a Curve.

                R2_Curve query_curve_dataset[total_query_curves];
                store_as_curves_r2(query_file, query_points, query_curve_dataset); //Store As a Curve.

                Curve_Lsh Lsh_curves = Curve_Lsh(L, total_input_curves,total_query_curves,number_of_points, k, w, delta, curve_dataset,query_curve_dataset,output_file);

            } else {
                cout << "WARNING! Wrong metric given for the Frechet algorithm.\nPlease give correct metric:    " << endl;
                cin >> metric;
                
                continue;
            }

        } else if(algorithm.compare("LSH") == 0 ) {
            store_curves(input_file, number_of_points, dataset);
            // Create query dataset
            store_curves(query_file, query_points, queries);
            
            Lsh<double> lsh = Lsh<double>(L, total_input_curves, number_of_points, k, w,dataset);

            lsh.ANN(queries, total_query_curves, dataset, total_input_curves, 1, output_file, 0);
        
        } else if(algorithm.compare("Hypercube") == 0) {
            store_curves(input_file, number_of_points, dataset);
            // Create query dataset
            store_curves(query_file, query_points, queries);
            // create hypercube 
            hypercube<double> cube = hypercube<double>(probes, M, w, k, number_of_points, total_input_curves, dataset);

            cube.ANN(queries, total_query_curves, dataset, total_input_curves, 1, output_file, 0);

        } else {
            cout << "WARNING! Wrong metric given for the Frechet algorithm.\nPlease give correct algorithm:    " << endl;
            cin >> algorithm;
                
            continue;
        }

        string choice;
        cout << "Would you like to continue (Y/N) ?   ";

        do {
            cin >> choice;
        } while(choice != "Y" && choice != "N");

        stop = (choice == "Y");

        if(stop) {
            cout << "New query file:    "; cin >> query_file;
            cout << "\nNew algorithm:    "; cin >> algorithm;
            if(algorithm.compare("Frechet") == 0) {
                cout << "\nNew metric:    "; cin >> metric;
            }

            // get info of query file
            calculate_file_info(&total_query_curves, &query_points, &query_file);
            // Create query dataset
            store_curves(query_file, query_points, queries);
        }
    } while(stop);


    return 1;
} 

