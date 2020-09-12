#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "Stat.h"
#define PI 3.1415927
using namespace std;
class station
{
	//Normalized ground motion, PGA=1.

public:
	station();
	string name;
	double x;
	double y;
	int lithoClass;
	double C;
	double C_std;
	double phi;
	double phi_std;
	double gamma;
	double Ac;
	double slope;
	Stat *stat;
	void getLithoC_phi();
	double getAc(double slope1);
	double getAcBy(double c_input,double phi_input, double slope1,double m);
	void GetPhi_C(int n, vector<double>&C_cql, vector <double>&phi_cql);
};