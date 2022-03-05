#pragma once
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
using namespace std;
#include "./Point.hpp"
//Representation of a R2_Curve on R^2 space.

class R2_Curve{

    public:
    long int id;
    string name;
    vector<Point> v;



    R2_Curve() {
        id = -1;
    }
    R2_Curve(const R2_Curve &obj) : id(obj.id), v(obj.v) ,name(obj.name) {}
    //R2_Curve(const long int &item_id, const  &key ) : id(item_id) { v.push_back(key); }
    ~R2_Curve() { v.clear(); }
    void set_id(const long int &item_id) { id = item_id; }
    void setPoint( const Point &d) {  v.push_back(d); }
    void setVector( const Point &d) {  v.push_back(d); }
    vector<Point> getVector() { return v; }
    long int get_id() { return id; }
    string get_name() { return name; }
    void set_name(string _name){ name = _name ;}
    void print_r2_curve() {
        cout << "id = " << id << endl;
        for(int i = 0; i < v.size(); i++) {
            v.at(i).print_point();
            cout<<" , ";
        } 
        cout << endl;
    }
    int get_number_of_points(){
        return v.size();
    }
    void PrintID(){
        cout<<"item_id:"<<id<<endl;
        return ;
    }

};