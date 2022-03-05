#include <string>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <ctime>
#include <ctype.h>
#include <random>
#include <chrono>
#include <math.h>
#include <map>
#include <fstream>

#include  "../headers/Functionality.hpp"
#include "../headers/Data/Data.hpp"
#include "../headers/Point.hpp"

#include "../headers/frechet/frechet.hpp"
#include "../headers/frechet/curve.hpp"
#include "../headers/frechet/types.hpp"
#include "../headers/frechet/point.hpp"

using namespace std;

using namespace Frechet;
using namespace Continuous;
//Function to parse arguments for nearest neighbor search.
void parse_arguments_nn(int argc, char **argv, 
    string *i_file, string *qr_file, int *k, int *L, int *M, 
    int *probes, string * out_file, string * algo, string* metric ,double* d){
        //Lets parse the args.
        if( argc<16 || argc>21 ){
            cerr << "Wrong arguments!" << endl;
        }

        //**** I Will code it later *****//
        if(strcmp(argv[1], "-i") != 0) {    
        cerr << "Please give the input file." << endl;
        exit(-1);
        }
        *i_file = argv[2];

        if(strcmp(argv[3], "-q") != 0) {    
            cerr << "Please give query file." << endl;
            exit(-1);
        } 
        *qr_file = argv[4];

        if(strcmp(argv[5], "-k") != 0) {    
            cerr << "Please give k" << endl;
            exit(-1);
        } 

        *k = stoi(argv[6]);
        
        if(strcmp(argv[7], "-L") != 0) {    
            cerr << "Please give L." << endl;
            exit(-1);
        } 
        
        *L = stoi(argv[8]);

        if(strcmp(argv[9], "-M") != 0) {    
            cerr << "Please give M." << endl;
            exit(-1);
        } 
        
        *M = stoi(argv[10]);
        
        if(strcmp(argv[11], "-probes") != 0) {    
            cerr << "Please give -complete value 0 or 1." << endl;
            exit(-1);
        } 
        
        *probes = stoi(argv[12]);

        if(strcmp(argv[13], "-o") != 0) {    
        cerr << "Please give the output file." << endl;
        exit(-1);
        }
        *out_file = argv[14];

        if(strcmp(argv[15], "-algorithm") != 0) {    
        cerr << "Please give the correct algorithm." << endl;
        exit(-1);
        }
        *algo = argv[16];

        for(int i = 17; i < argc; i++) {
            if(strcmp(argv[i],"-metric") == 0) {
                *metric = argv[i + 1];
            } else if(strcmp(argv[i],"-delta") == 0) {
                *d = stod(argv[i+1]);
            } 
        }

        return;
    }

void cluster_parse_args(int argc, char **argv, string *input, string *config,
                string *output, string *update, string *assignment, string *complete)
{
     if(argc < 11) {
        cerr << "Usage of : " << argv[0] << " -i <input_file> -c <configuration_file> -o <output_file> -update <Mean Frecher or Mean Vector> -assignment <Classic or LSH or Hypercube or LSH_Frechet>"  << endl;
        exit(-1);
    }

    if(strcmp(argv[1], "-i") != 0) {    
        cerr << "Please give the input file." << endl;
        exit(-1);
    }
    *input = argv[2];

    if(strcmp(argv[3], "-c") != 0) {    
        cerr << "Please give config file." << endl;
        exit(-1);
    } 
    *config = argv[4];

    if(strcmp(argv[5], "-o") != 0) {    
        cerr << "Please give config file." << endl;
        exit(-1);
    } 

    *output = argv[6];
    if(output->empty()) {
        cerr << "No output file given.";
        exit(-1);    
    }

    if(strcmp(argv[7], "-update") != 0) {    
        cerr << "Please give config file." << endl;
        exit(-1);
    } 
    
    *update = argv[8];

    if(strcmp(argv[9], "-assignment") != 0) {    
        cerr << "Please give config file." << endl;
        exit(-1);
    } 
    
    *assignment = argv[10];
    
    if(strcmp(argv[11], "-complete") != 0) {    
        cerr << "Please give -complete value 0 or 1." << endl;
        exit(-1);
    } 
    
    *complete = argv[12];
}

// Function to read configuration balues
int read_config_values(const string &config_file, int *number_of_clusters, 
                int *number_of_hash_tables,int *number_of_vector_hash_functions, 
                int *max_number_M_hypercube, int *number_of_hypercube_dimensions, int *number_of_probes)
{
    ifstream config(config_file);  
    if(config.is_open()) {
        string str;
        std::string line;
        
        while(getline(config, line)) 
        {
            if (line[0] == '#' || line.empty()) continue;

            auto delimiterPos = line.find(":");
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);


            if (name == "number_of_clusters") *number_of_clusters = atoi(value.c_str());
            else if (name == "number_of_hash_tables") *number_of_hash_tables = atoi(value.c_str());
            else if (name == "number_of_vector_hash_functions") *number_of_vector_hash_functions = atoi(value.c_str());
            else if (name == "max_number_M_hypercube") *max_number_M_hypercube = atoi(value.c_str());
            else if (name == "number_of_hypercube_dimensions") *number_of_hypercube_dimensions = atoi(value.c_str());
            else if (name == "number_of_probes") *number_of_probes = atoi(value.c_str());
        }
    }else{
        return -1 ;
    }

    return 1;
}

//Function to get the Number of poin
int get_num_of_points(string str) {
    stringstream ss(str);
    string temp;
    int countDim = 0;
    while(getline(ss, temp, '\t')) {
            // cout<<temp.at(0)<<endl;
            countDim++;
    }

    return countDim - 1; 
}



//Calculate needed info for the file dataset.
void calculate_file_info(int *total_curves, int *number_of_points, string *filename)
{   *total_curves = 0;
    *number_of_points = 0;
    bool flag = true;
    ifstream input_file(*filename); 
    //Check the filename.
    while(!input_file) {
        cout << "Please give us a valid filename."<<endl;        
        string new_file;
        cin >> new_file;
        input_file.open(new_file);
        //Get a new one and update filename.
        filename->clear();
        *filename = new_file;
    }
    //--------------------------------//
    string str;
    while(getline(input_file, str)) {
        //Only once : ckeck the number of points of the TimeSeries.
        if(flag) {
            flag = false;
                *number_of_points = get_num_of_points(str);
        }
        //Count total number of Curves(Time Series).
        (*total_curves)++;    
    }
    //Close the file.
    input_file.close();
}

/*
//Function to store the Curves from file.
void store_curves(string filename, int num_of_points, Curve<double> *arr) {
    //Open file.
    ifstream input_file(filename);  
    int curve_num = 0; //Count the number of curves stored.
    if(input_file.is_open()) {
        string str;
        while(getline(input_file, str)) {
            // we need a counter for vector dimensions
            int i = 0;
            istringstream ss(str);
            //We need to seperate coordinates and store them in a vector container.
            string name;
            getline(ss, name, '\t');  
            arr[curve_num].set_name(name);
            arr[curve_num].set_id(curve_num);

            string x_ij;
            while(getline(ss, x_ij, '\t')) {
                // not a whitespace char
                if(isspace(x_ij.at(0)) == 0)
                    if(i++ < num_of_points) {
                        // let's set our data
                        arr[curve_num].setVector(stod(x_ij));
                    }
            }
            curve_num++;
        }
    } else {
        cerr << "Unable to open file" << endl;
        exit(-1);
    }
    //Close the file. 
    input_file.close();
    return;
}

*/

//Function to store the Curves from file.
void store_curves(string filename, int num_of_points, Data<double> *arr) {
    //Open file.
    ifstream input_file(filename);  
    int curve_num = 0; //Count the number of curves stored.
    if(input_file.is_open()) {
        string str;
        
        while(getline(input_file, str)) {
            // we need a counter for vector dimensions
            int i = 0;
            istringstream ss(str);
            //We need to seperate coordinates and store them in a vector container.
            string name;
            getline(ss, name, '\t');  
            arr[curve_num].set_name(name);
            arr[curve_num].set_id(curve_num);

            string x_ij;
            while(getline(ss, x_ij, '\t')) {
                // not a whitespace char
                if(isspace(x_ij.at(0)) == 0)
                    if(i++ < num_of_points) {
                        // let's set our data
                        arr[curve_num].setVector(stod(x_ij));
                    }
            }
            curve_num++;
        }
    } else {
        cerr << "Unable to open file" << endl;
        exit(-1);
    }
    //Close the file. 
    input_file.close();
    return;
}


void store_as_curves_r2(string filename, int num_of_points, R2_Curve *arr){
    //Open file.
    ifstream input_file(filename);  
    int curve_num = 0; //Count the number of curves stored.
    if(input_file.is_open()) {
        string str;
        while(getline(input_file, str)) {
            // we need a counter for vector dimensions
            int i = 0;
            istringstream ss(str);
            //We need to seperate coordinates and store them in a vector container.
            string name;
            getline(ss, name, '\t');  
            arr[curve_num].set_name(name);
            arr[curve_num].set_id(curve_num);

            string y_ij;
            while(getline(ss, y_ij, '\t')) {
                // not a whitespace char
                if(isspace(y_ij.at(0)) == 0)
                    if(i++ < num_of_points) {
                        // let's set our data
                        double x_ij = (double)i;
                        arr[curve_num].setVector(my_point(x_ij,stod(y_ij)));
                    }
            }
            curve_num++;
        }
    } else {
        cerr << "Unable to open file" << endl;
        exit(-1);
    }
    //Close the file. 
    input_file.close();
    return;
}


void normal_distribution_fun(double *n, float x, float y) {
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count(); 
    default_random_engine e (seed); 
  
    /* declaring normal distribution object 'distN' and initializing its mean and standard deviation fields. */
    /* Mean and standard deviation are distribution parameters of Normal distribution. Here, we have used mean=5, and standard deviation=2. You can take mean and standard deviation as per your choice */
    normal_distribution<double> distN(x, y);
    *n  = distN(e);
}

int modular_pow(int base, int exponent, int modulus)
{
    int result = 1;
    base = base % modulus;
    while (exponent > 0)
    {
        if(exponent & 1)
            result = (result * base) & modulus;
        exponent = exponent >> 1;
        base = (base*base) % modulus;
    }
    return result;
}

long long inn_product(int *h, int *r, int dim) {
    int product = 0;
    for(int i = 0; i < dim; i++) 
        product += h[i] * r[i];

    return product;
}   

unsigned long positive_modulo( unsigned long x, unsigned y) {
    return  ( ( x % y ) + y ) % y;
}

int coinToss() {
    random_device rd;
	mt19937 gen(rd());
  
        /* declaring normal distribution object 'distN' and initializing its mean and standard deviation fields. */
        /* Mean and standard deviation are distribution parameters of Normal distribution. Here, we have used mean=5, and standard deviation=2. You can take mean and standard deviation as per your choice */
    uniform_int_distribution<int> distN(0, 1);

    return distN(gen);
}

int hammingDistance(int n1, int n2){
    int x = n1 ^ n2;
    int setBits = 0;

    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }

    return setBits;
}


/**
 * @brief 
 * Function that calculates the theoritical value of ε needed for filtering 
 * 
 * @param total_curves 
 * @param curve_points 
 * @return double 
 */
double calclulate_mean_theoritical_e(int number_of_points, Data<double> curve) {
    double a = 0.0, b = 0.0;
    double sum = 0.0;

    for(int i = 0; i < number_of_points; i++) {
        a = curve.getVector().at(i);

        if(i == number_of_points - 1) {
            break;
        } else 
            b = curve.getVector().at(i+1);
        
        sum += abs(a-b);
    }
    // return a percent of the difference
    return  1.2 *(sum / number_of_points);
}

/**
 * @brief 
 * Filtering function which keeps only the important points.
 * For any consecutive points a,b,c, if |a-b| <= e and |b-c| <= e
 * then remove b. Continue with the resulting sequence 
 * 
 * @param total_curves 
 * @param curve_points 
 */
void filtering_curves( int total_curves, int number_of_points, Data<double>* old_curves, Data<double>* new_curves) {
    // these are the three consecutibe points
    double a = 0, b = 0, c = 0; 

    for(int i = 0; i < total_curves; i++) {
        // set id and name
        new_curves[i].set_id(old_curves[i].get_id());
        new_curves[i].set_name(old_curves[i].get_name());
        // this is our e
        double e = calclulate_mean_theoritical_e(number_of_points, old_curves[i]);
        // filter every curve
        for(int j = 0; j < number_of_points; j += 3) {
            a = old_curves[i].getVector().at(j); 
            if(j == number_of_points - 1) {
                new_curves[i].setVector(a);

                break;
            } else if( j == number_of_points - 2)  {
                b = old_curves[i].getVector().at(j+1);

                new_curves[i].setVector(a); new_curves[i].setVector(b);

                break;
            }
            b = old_curves[i].getVector().at(j+1);
            c = old_curves[i].getVector().at(j+2);
            
            if( (abs(a-b) <= e ) && (abs(b-c) <= e) ) {
                new_curves[i].setVector(a);  
                new_curves[i].setVector(c);
            } else {
                new_curves[i].setVector(a);
                new_curves[i].setVector(b);
                new_curves[i].setVector(c);
            }
        }
    }
}
/**
 * @brief 
 * Brute force function for Continuous Frechet.
 * At first we initialise Frechet structure and we save all distances
 * between the query curve and all filtered curves.
 * 
 * @param query_vector_curve 
 * @param qr_name 
 * @param input_curves 
 * @param number_of_input_curves 
 * @return vector<pair<double, string>> 
 */

vector<pair<double, string>> brute_force_continuous_frechet(vector<double> query_vector_curve, string qr_name, Data<double> *input_curves, int number_of_input_curves)
 {
    // we need to initialize our frechet structures

    Point p1(1);
    Points total_input_points(1);

    Curve c1(total_input_points.dimensions());

    
    Points total_query_points(1);
    Point p2(1);

     // do this once for query curve
    for(size_t n = 0; n < query_vector_curve.size(); n++) {
        p2.set(0, query_vector_curve.at(n));
        total_query_points.add(p2, p2.dimensions());
    } 
    Curve c2(total_query_points, qr_name);

    vector<pair<double,string>> result;

    for(int i = 0; i < number_of_input_curves; i++) {
        vector<double> filtered_input_curve_vector = input_curves[i].getVector();

        c1.set_name(input_curves[i].get_name());
        
        for(size_t j = 0; j < filtered_input_curve_vector.size(); j++) {
            p1.set(0, filtered_input_curve_vector.at(j));
            // total_input_points.add(p1, p1.dimensions());
            c1.push_back(p1);
        }
        
        // distance
        Distance d = distance(c1, c2);

        result.push_back(pair<double,string>(d.value, input_curves[i].get_name()));

        total_input_points.clear(); c1.clear(); 
    }

    p2.clear(); total_query_points.clear(); c2.clear();

    return result;
}







