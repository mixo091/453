#include <iostream>
#include <vector>
#include <math.h>
#include "./Point.hpp"
using namespace std;

// Calculate Max .//
double maximum(double a, double b){	
	if (a >= b)
    	return a; 
  	return b;
}

// Calculate Min .//
double minimum_of_3(double a, double b, double c)
{
	if (a <= b && a <= c)
      	return a;
  	else if (b <= a && b <= c)
      	return b;
    else
      	return c;
}


double Euclidean_Distance(const my_point & v1,const my_point &v2){
 		double sum = 0;
		//cout<<" v1.y_1 "<<v1.y_i<<" v2.y_2 "<<v2.y_i<<endl;
		//cout <<" v1.y1"<< pow((v1.y_i - v2.y_i),2) <<endl;
 		sum = pow((v1.y_i - v2.y_i),2) + pow((v1.x_i - v2.x_i),2);
		//cout<<"sum " <<sum<<endl;
 	
 	return sqrt(sum);
 };

double Discrete_Frechet(const std::vector<my_point> & v1,const std::vector<my_point> &v2){
	int n = v1.size();
	int m = v2.size();
	cout<<" curve_size ["<<n<<","<<"]"<<endl;
	double ** Calculation;
	Calculation = (double**) malloc(n * sizeof(double*));
	if (Calculation == NULL){
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	for (int i = 0 ; i < n ; i++){
    	Calculation[i] = (double*) malloc(m * sizeof(double));	
		if (Calculation[i] == NULL){
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}
	Calculation[0][0] = Euclidean_Distance(v1.at(0), v2.at(0));
	cout<<"Calculation[0][0] : "<<Calculation[0][0] <<endl;

	for (int i = 1 ; i < m ; i++)
    	Calculation[i][0] = maximum(Calculation[i-1][0],  Euclidean_Distance(v1.at(i), v2.at(0)));
	
    for (int j = 1 ; j < n ; j++)
    	Calculation[0][j] = maximum(Calculation[0][j-1],  Euclidean_Distance(v1.at(0), v2.at(j)));

	for (int i = 1 ; i < m ; i++)
    	for (int j = 1 ; j < n ; j++)
			Calculation[i][j] = maximum(minimum_of_3(Calculation[i-1][j], Calculation[i][j-1], Calculation[i-1][j-1]), Euclidean_Distance(v1.at(i), v2.at(j)));

	double result = Calculation[m-1][n-1];
	
    for (int i = 0 ; i < m ; i++)
    	free(Calculation[i]);      
    free(Calculation);
	
	cout <<" Result Distance :"<<result<<endl;
	return result;
	
}