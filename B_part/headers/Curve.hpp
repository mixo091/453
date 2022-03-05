#pragma once
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>

using namespace std;
//Represents a Curve_point in space.
template <typename T>
struct Curve
{
    long int id;
    string name;
    vector<T> v;

    public:
        int  Clustered_already_by = -1;
        double ClusteredInRadius ;

    Curve() {
        id = -1;
    }
    Curve(const Curve<T> &obj) : id(obj.id), v(obj.v) ,name(obj.name) {}
    Curve(const long int &item_id, const T &key ) : id(item_id) { v.push_back(key); }
    ~Curve() { v.clear(); }
    void set_id(const long int &item_id) { id = item_id; }
    void setVector( const T &d) {  v.push_back(d); }
    vector<T> getVector() { return v; }
    long int get_id() { return id; }
    string get_name() { return name; }
    void set_name(string _name){ name = _name ;}
    void print_curve() {
        cout << "id = " << id << endl;
        for(auto it = v.begin(); it < v.end(); it++) {
            cout << *it << ",";
        } 
        cout << endl;
    }
    void PrintID(){
        cout<<"item_id:"<<id<<endl;
        return ;
    }
};