#include "station.h"

station::station()
{
	stat = new Stat;
}


void station::getLithoC_phi()
{
	switch (lithoClass)
	{
	case 1:
		C = 35;
		C_std = 10;
		phi =40 ;
		phi_std=4;
		gamma = 25;
		break;
	case 2:
		C = 30;
		C_std = 9;
		phi = 35;
		phi_std = 3;
		gamma = 23;
		break;
	case 3:
		C = 25;
		C_std = 7;
		phi = 30;
		phi_std = 3;
		gamma = 20;
		break;
	case 4:
		C = 15;
		C_std = 5;
		phi = 20;
		phi_std = 2;
		gamma = 16;
		break;
	}
}

double station::getAc(double slope1)
{
	double Fs = 0.0;
	double t = 2.4;
	double m = 0.0;
	string filename = ".//check//" + to_string(slope1) + "check.txt";

	double C_temp = -1;
	double phi_temp = -1;
	while (C_temp<0)
	{
		C_temp = stat->gaussrand(C, C_std);
	}
	while (phi_temp < 0)
	{
		phi_temp = stat->gaussrand(phi, phi_std);
	}
	slope = slope1;
	double temp = C_temp / gamma / sin(slope / 180 * PI)/t;
	Fs = C_temp / gamma / sin(slope / 180 * PI)/t + tan(phi_temp / 180 * PI) / tan(slope / 180 * PI) - m*10*tan(phi_temp / 180 * PI) / gamma / tan(slope / 180 * PI);
	Ac = (Fs-1)*sin(slope / 180 * PI)*9.8; //m/s/s
	return Ac;
}

double station::getAcBy(double c_input, double phi_input,double slope1, double m)
{
	double Fs = 0.0;
	double t = 2.4;
	string filename = ".//check//" + to_string(slope1) + "check.txt";
	ofstream opt(filename.c_str(), ios::app);
	double C_temp = c_input;
	double phi_temp = phi_input;
	slope = slope1;

	Fs = C_temp / gamma / sin(slope / 180 * PI) / t + tan(phi_temp / 180 * PI) / tan(slope / 180 * PI) - m * 10 * tan(phi_temp / 180 * PI) / gamma / tan(slope / 180 * PI);

	Ac = (Fs - 1)*sin(slope / 180 * PI)*9.8; //m/s/s
	opt << C_temp << "\t" << phi_temp << "\t" << Fs <<  "\t"<<Ac << endl;
	return Ac;
}


void station::GetPhi_C(int n,vector<double>&C_cql,vector <double>&phi_cql)
{

	for (int i=0;i<n;i++)
	{
		double C_temp = -1;
		double phi_temp = -1;
		while (C_temp < 0)
		{
			C_temp = stat->gaussrand(C, C_std);
		}
		while (phi_temp < 0)
		{
			phi_temp = stat->gaussrand(phi, phi_std);
		}
		C_cql.push_back(C_temp);
		phi_cql.push_back(phi_temp);
	}
}