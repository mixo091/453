#pragma once
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>

//Representation of a R2_Curve on R^2 space.
#include "./Point.hpp"

using namespace std;

class R2_Curve{
    long int id;
    string name;
    vector<my_point> v;
    public:


    R2_Curve() {
        id = -1;
    }
    R2_Curve(const R2_Curve &obj) : id(obj.id),name(obj.name), v(obj.v)  {}
    //R2_Curve(const long int &item_id, const  &key ) : id(item_id) { v.push_back(key); }
    ~R2_Curve() { v.clear(); }

    void set_id(const long int &item_id) { id = item_id; }
    void setVector( const my_point &d) {  v.push_back(d); }
    void set_name(string _name){ name = _name ;}

    vector<my_point> getVector() { return v; }
    long int get_id() { return id; }
    string get_name() { return name; }
    
    void print_r2_curve() {
        cout << "id = " << id << endl;
        for(size_t i = 0; i < v.size(); i++) {
            v.at(i).print_point();
            cout<<" , ";
        } 
        cout << endl;
    }
    void PrintID(){
        cout<<"item_id:"<<id<<endl;
        return ;
    }
};