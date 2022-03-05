#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <math.h>

using namespace std;

#include "ClusteringUtilities.hpp"

Data<double>* parseData(string filename, int dim, int totalVectors) {
    // open file for reading
    ifstream input_file(filename);  
    Data<double> *arr = new Data<double>[totalVectors];
    // we need to count the Data created
    int nVector = 0;
    if(input_file.is_open()) {
        string str;
        while(getline(input_file, str)) {
            // we need a counter for vector dimensions
            int i = 0;
            istringstream ss(str);
            // we need to seperate coordinates and store them in a vector container
            // store id
            string id;
            getline(ss, id, ' ');  
            arr[nVector].setId(stoi(id));

            string x_ij;
            while(getline(ss, x_ij, ' ')) {
                // not a whitespace char
                if(isspace(x_ij.at(0)) == 0)
                    if(i++ < dim) {
                        // let's set our data
                        arr[nVector].setVector(stod(x_ij));
                    }
            }
            nVector++;
        }
    } else {
        cerr << "Unable to open file" << endl;
        exit(-1);
    }
    // close the file 
    input_file.close();

    return arr;
}

int getVectorDim1(string str) {
    stringstream ss(str);
    string temp;
    int countDim = 0;
    while(getline(ss, temp, ' ')) {
        if(isspace(temp.at(0)) == 0)
            countDim++;
    }

    return countDim - 1; 
}

void CalculateInfo(int *tableSize, int *dim, string filename)
{
    bool flag = true;
    ifstream input_file;   
    input_file.open(filename.c_str());
    if(input_file.is_open()) {
        string str;
        while(getline(input_file, str)) {
            if(flag) {
                flag = false;
                *dim = getVectorDim1(str);
            }
            (*tableSize)++;    
        }
    } else {
        cerr << "Unable to open file" << endl;
        exit(-1);
    }
    //Close the File. 
    input_file.close();
}

//Print a set of points in d-space.
void printSet(Data<double>* set ,int size)
{
    for (int i = 0; i < size; i++) {
        set[i].printVector();
    }
} 


