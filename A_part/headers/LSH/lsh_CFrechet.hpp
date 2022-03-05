#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "./lsh.hpp"

template <typename T>
class lsh_CFrehet : public Lsh<T> {
    double delta;
    int number_of_points;
    int number_of_curves;
    double t;

    public:
        lsh_CFrehet(double d, int num_points, int num_curves, 
                int L, int _k, int _w,
                Data<T> *curves )
        : 
        Lsh<T>(L, num_curves / BUCKET_DIVIDER, num_points, _k, _w),
        delta(d), 
        number_of_points(num_points),
        number_of_curves(num_curves)
        {
            default_random_engine generator;
            uniform_real_distribution<double> distribution(0.0, delta);
            // initialise t for snapping
            t = distribution(generator);

            // vector of snapped points
            vector<double> snapped_points;

            // snapping, padding for input curves
            for(int i = 0; i < number_of_curves; i++) {                
                vector<double> filtered_curve_vector = curves[i].getVector();
                
                // string name = curves[i].get_name();
                // string id = curves[i].get(id);
                snapped_points = snapping(filtered_curve_vector);
                padding(&snapped_points, number_of_points);

                // cout << "curves after snapping and padding size = : "  << snapped_points.size() << endl;
                // our key now is snapped_points vector
                // insert it to hash table
                this->hash_tables[0]->insert(snapped_points, &curves[i]);

                snapped_points.clear();  
            }    
        }

        void padding(vector<T>* final_vector, int total_points) {
            random_device rd;
            default_random_engine eng(rd());
            uniform_real_distribution<double> distr(0.0, 1000);

            double padding_value = distr(eng);

            for(size_t i = final_vector->size(); (int) i < total_points; i++)
                final_vector->push_back(padding_value);
        }

        vector<T> snapping(vector<T> filtered_curve) {
            vector<T> snapped_points;
            vector<T> to_delete;

            for(size_t j = 0; j < filtered_curve.size(); j++) {
                    double y_i = filtered_curve.at(j);
                    double new_value = floor( (y_i + this->t)/((this->delta))) * this->delta;

                    snapped_points.push_back(new_value);
                }
                double max = 0.0;
                double min = 0.0;

                // find minima-maxima
                for(size_t m = 1; m < snapped_points.size(); m++) {
                    
                    if(m != snapped_points.size() - 1)
                    {   // avoid overflow
                        if(snapped_points.at(m+1) >= snapped_points.at(m-1)) {
                            max = snapped_points.at(m+1);
                            min = snapped_points.at(m-1);
                        } else {
                            max = snapped_points.at(m-1);
                            min = snapped_points.at(m+1);
                        }
                        // check if current point is between this space
                        if(snapped_points.at(m) >= max && snapped_points.at(m) <= min) {
                            // min_maxima_vector.push_back(min); min_maxima_vector.push_back(max);
                            to_delete.push_back(snapped_points.at(m));
                        } 
                    }
                }

                // delete "to_delete" values from snapped_points vector
                for(size_t m = 0; m < to_delete.size(); m++) {
                    // cout << "to delete " << to_delete.at(m) << endl;
                    std::vector<double>::iterator position = std::find(snapped_points.begin(), snapped_points.end(), to_delete.at(m));
                    if (position != snapped_points.end()) // == myVector.end() means the element was not found
                        snapped_points.erase(position);
                }

            return snapped_points;
        }

        void perforfm_continuous_frechet(const string& out_file, 
            int queries_number, int qrcurve_total_points, Data<double> *queries, 
            Data<double> *input_curves) 
            {
             // open output file
             ofstream output(out_file);
             if(!output.is_open()) {
                 cout << "Error with output file\n";
                 exit(-1);
            }
            vector<double> snapped_points;
            // used for nearest neighbor
            std::pair<double,string> nn_dist_id_pair(-1,"");

            // for the output
            map<double, string> result_map;
            
            double continuous_frechet_avg_time = 0;
            double avg_durationBF = 0;
            vector<double> vector_MAF;

            // snapping, padding for query curves
            for(int i = 0; i < queries_number; i++) {
                vector<double> filtered_query_curve_vector = queries[i].getVector();

                snapped_points = snapping(filtered_query_curve_vector);
                padding(&snapped_points, qrcurve_total_points);
             
                // hash snapped_points vector to our hash_table
                this->hash_tables[0]->search_NN(snapped_points, &queries[i], &nn_dist_id_pair, &continuous_frechet_avg_time);

                if(nn_dist_id_pair.first != -1 && nn_dist_id_pair.second != "") {
                    result_map.insert(nn_dist_id_pair);

                    cout << "mpike sto result" << endl;
                }

                // reset to -1
                nn_dist_id_pair.first = -1;
                nn_dist_id_pair.second = "";
                
                snapped_points.clear();
                   
                // brute force method
                const auto startBF = high_resolution_clock::now();
                // vector used for brute force
                vector<pair<double, string>> brute_force_v;
                brute_force_v = brute_force_continuous_frechet(queries[i].getVector(), queries[i].get_name(), input_curves, number_of_curves);
                sort(brute_force_v.begin(), brute_force_v.end());
    
                const auto endBF = high_resolution_clock::now();
                auto duration_bf = duration_cast<microseconds>(endBF - startBF);

                avg_durationBF += duration_bf.count();

                // build output
                
                output << "Query: " << queries[i].get_name();
                output<< "\nAlgorithm:  { " << "LSH_Frechet_Continuous }";
                output << "\nApproximate Nearest neighbor:  ";
                if(result_map.size() > 0) output << result_map.cbegin()->second;
                output << "\nTrue Nearest neighbor:    " ; 
                if(brute_force_v.size() > 0) output << brute_force_v.at(0).second;
                output << "\ndistanceApproximate:    "; 
                if(result_map.size() > 0) output << result_map.cbegin()->first;
                output << "\ndistanceTrue:  "; 
                if(brute_force_v.size() > 0) output << brute_force_v.at(0).first;
                output << endl << endl;
                
                 if(result_map.size() > 0 && brute_force_v.size() > 0) {
                    vector_MAF.push_back(result_map.cbegin()->first / brute_force_v.at(0).first);
                }

                result_map.clear();
            }
            output << "\n\ntApproximateAverage:    " << (double) continuous_frechet_avg_time / queries_number;
            output << "\ntTrueAverage:    " << (double) (avg_durationBF / 1000000) / queries_number; 
            output << "\nMAF:    " << *max_element(vector_MAF.begin(), vector_MAF.end()) << endl;

            vector_MAF.clear();

            output.close();
        } 

        ~lsh_CFrehet() {

        }
};