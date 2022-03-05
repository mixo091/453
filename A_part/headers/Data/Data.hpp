#pragma once

#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>

using namespace std;

//Represents a Data_point in space.
template <typename T>
struct Data
{
    long int id;
    string name;
    vector<T> v;

    public:
        int  Clustered_already_by = -1;
        double ClusteredInRadius;

    Data() {
        id = -1;
    }
    Data(const Data<T> &obj) : id(obj.id), v(obj.v) {}
    Data(const long int &item_id, const T &key ) : id(item_id) { v.push_back(key); }
    Data(const long int &item_id, const vector<T> &key ) : id(item_id), v(key) {  }
    ~Data() { v.clear(); }
    
    /* setters */
    void setId(const long int &item_id) { id = item_id; }
    void set_name(const string &_name) { name  = _name;}
    void set_id(const long int &item_id) { id = item_id; }
    void setVector( const T &d) {  v.push_back(d); }

    /* getters */
    vector<T> getVector() { return v; }
    long int get_id() { return id; }
    string get_name() { return name; }

    void printVector() {
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