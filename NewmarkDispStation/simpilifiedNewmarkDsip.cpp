#include <iostream>
#include <windows.h>
#include <sstream>
#include <fstream>
#include<string> 
#include<vector>
#include <math.h>
#include "groundmotion.h"
#define PI 3.1415927
using namespace std;
double simpilifiedNewmarkDsip_Jibson07_1(double ay, double PGA)
{
	double disp = 0.0;
	disp = 0.215 + log10(pow((1 - ay / PGA), 2.341)*pow((ay / PGA), -1.438));
	return exp(disp);
}

double simpilifiedNewmarkDsip_Jibson07_2(double ay, double PGA, double M)
{
	double disp = 0.0;
	disp = -2.71 + log10(pow((1 - ay / PGA), 2.335)*pow((ay / PGA), -1.478)) + 0.424*M;
	return pow(10, disp);
}

double simpilifiedNewmarkDsip_Jibson07_3(double ay, double Ia)
{
	double disp = 0.0;
	disp = 2.401*log10(Ia) - 3.481*log10(ay) - 3.230;
	return pow(10, disp);
}
double simpilifiedNewmarkDsip_Jibson07_4(double ay, double PGA, double Ia)
{
	double disp = 0.0;
	disp = 0.561*log10(Ia) - 3.833*log10(ay / PGA) - 1.474;
	return pow(10, disp);
}

//double simpilifiedNewmarkDsip_Saygili08_1(double ay, double PGA, double Ia)
//{
//	double disp = 0.0;
//	disp = 0.561*log10(Ia) - 3.833*log10(ay / PGA) - 1.474;
//	return pow(10, disp);
//}
double simpilifiedNewmarkDsip_Hsieh11_1(double ay, double PGA, double Ia)
{
	double disp = 0.0;
	disp = 11.287*ay*log10(Ia) - 11.485*ay+1.948;
	return pow(10, disp);
}

double simpilifiedNewmarkDsip_Hsieh11_2(double ay, double PGA, double Ia)
{
	double disp = 0.0;
	disp = 0.847 *log10(Ia) - 10.62*ay+6.587*ay*log10(ay) + 1.84;
	return pow(10, disp);
}

double simpilifiedNewmarkDsip_Saygili08_1(double ay, double PGA)
{
	double disp = 0.0;
	disp = 5.52 - 4.43*(ay / PGA) - 20.39*pow((ay / PGA), 2) + 42.61*pow((ay / PGA), 3) - 28.74*pow((ay / PGA), 4) + 0.72*log(PGA);
	return exp(disp);
}

double simpilifiedNewmarkDsip_Saygili08_2(double ay, double PGA,double PGV)
{
	double disp = 0.0;
	disp = -1.56-4.58*(ay / PGA) - 20.84*pow((ay / PGA), 2) + 44.75*pow((ay / PGA), 3) - 30.50*pow((ay / PGA), 4) -0.64*log(PGA)+1.55*log(PGV);
	return exp(disp);
}

double simpilifiedNewmarkDsip_Saygili08_3(double ay, double PGA, double PGV,double Ia)
{
	double disp = 0.0;
	disp = -0.74 - 4.93*(ay / PGA) - 19.91*pow((ay / PGA), 2) + 43.75*pow((ay / PGA), 3) - 30.12*pow((ay / PGA), 4) - 1.30*log(PGA) + 1.04*log(PGV)+0.67*log(Ia);
	return exp(disp);
}


