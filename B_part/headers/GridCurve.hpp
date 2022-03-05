#include <iostream>
#include <vector>
#include <string>
#include "./CurveR2.hpp"
using namespace std;

class Grid_Curve{
    private:
        vector<double> x;
    public:
        R2_Curve* original_curve;
        string name;

        int id; 
        Grid_Curve(){};
        Grid_Curve(string _name , int id , R2_Curve* fromCurve , vector<double> _x){
            this->name = _name;
            this->id = id; 
            this->original_curve = fromCurve;
            this->x = _x;
        };
        ~Grid_Curve(){};
        void set_id(int _id){this->id = _id ;};
        void set_name(string _name){ this->name = _name ;};
        int get_id(){return this->id;};
        string get_name(){return this->name ; };
        void push(double x_i){x.push_back(x_i);};
        vector<double> get_curve(){return x ; };
        void print_grid_curve(){
            cout<<" Vector_Curve [id:"<<id<<",name:"<<name<<"] :: {";
           for(int i=0; i < x.size(); i++)
                cout << x.at(i) << ' ';
            cout<<"}"<<endl;
        };
};


