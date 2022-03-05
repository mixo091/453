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
#include  "../headers/Functionality.hpp"
using namespace std;
//Function to parse arguments for nearest neighbor search.
void parse_arguments_nn(int argc, char **argv, 
    string *i_file, string *qr_file, int *k, int *L, int *M, 
    int *probes, string * out_file, string * algo, string* metric ,double* d){
        //Lets parse the args.
        if( argc<19 || argc>21 ){

        }

        //**** I Will code it later *****//


        return;
    }

//Function to get the Number of poin
int get_num_of_points(string str) {
    stringstream ss(str);
    string temp;
    int countDim = 0;
    while(getline(ss, temp, '\t')) {
            cout<<temp.at(0)<<endl;
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
                        arr[curve_num].setVector(Point(x_ij,stod(y_ij)));
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









