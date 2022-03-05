#pragma once
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include "../Data/Data.hpp"

using namespace std;

//Function to calculate some info for the programm.
void CalculateInfo(int *, int *,string);
//Function for creating the data set given the input file.
Data<double>* parseData(string filename, int dim, int totalVectors); 
void printSet(Data<double>* set ,int size);