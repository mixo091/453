/* Data Structure which represents a point of a Curve in Space */
#pragma once
#include <iostream>

using namespace std;

class my_point{
    public:
        double x_i ; //Represents Time.//
        double y_i ;
        my_point(){
            x_i = 0.0;
            y_i = 0.0;
        }
        my_point(double x, double y):x_i(x),y_i(y){};
        void set_x(double x){ x_i = x ;};
        void set_y(double y){ y_i = y ;};
        double get_x(){return x_i ;};
        double get_y(){return y_i ;};
        void print_point(){cout << " ("<<x_i<<","<<y_i<<")";};
};