#ifndef GROUNDMOTION_H
#define GROUNDMOTION_H

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#define PI 3.1415927
using namespace std;

class GroundMotion
{
    //Normalized ground motion, PGA=1.

public:
    GroundMotion();
    string name;
    vector<double> t;   //time
    vector<double> a;   //acceleration
	vector<double> v;   //velocity
    void ReadGM(const string path);
    double GetTotalDuration();
    double GetDt();
	double GetPGA();
	double GetIrias();
	double GetPGV();
};

#endif // GROUNDMOTION_H
