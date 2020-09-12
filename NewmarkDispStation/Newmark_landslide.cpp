/********************************************************
* @file    : Newmark_landslide.cpp
* @brief   : Calculate the Newmark rigid block displacement
* @details : Input the ground motion to calculate the displacement	
* @author  : Chinler
* @date    : 2018-12-15
*********************************************************/

#include <iostream>
#include <windows.h>
#include <sstream>
#include <fstream>
#include<string> 
#include<vector>
#include <math.h>
#include <stack>  
#include <iomanip>  
#include<io.h>
#include <omp.h>
#include <direct.h>
#include "groundmotion.h"
#include "Stat.h"
#include "station.h"
#define PI 3.1415927
using namespace std;
double NewmarkDisp(vector<double>&time, vector<double>&acc, double ay);
void getJustCurrentDir(string path, vector<string>& files);
wchar_t *multiByteToWideChar(const string& pKey);
//double simpilifiedNewmarkDsip_Jibson07_1(double ay, double PGA);
//double simpilifiedNewmarkDsip_Jibson07_2(double ay, double PGA, double M);
//double simpilifiedNewmarkDsip_Jibson07_3(double ay, double Ia);
//double simpilifiedNewmarkDsip_Jibson07_4(double ay, double PGA, double Ia);
//
//double simpilifiedNewmarkDsip_Hsieh11_1(double ay, double PGA, double Ia);
//double simpilifiedNewmarkDsip_Hsieh11_2(double ay, double PGA, double Ia);
//double simpilifiedNewmarkDsip_Saygili08_1(double ay, double PGA);
//double simpilifiedNewmarkDsip_Saygili08_2(double ay, double PGA, double PGV);
//double simpilifiedNewmarkDsip_Saygili08_3(double ay, double PGA, double PGV, double Ia);
void getFiles(string path, vector<string>& files)
{
	//文件句柄  
	long long  hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{ 
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}
void getFiles2GM(string path, vector<string>& files)
{
	//文件句柄  
	long long  hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				string temp = fileinfo.name;
				temp.erase(temp.length()-4, temp.length());
				if ((strcmp(temp.data(), "EW")==0)|| (strcmp(temp.data(), "NS") == 0))
				{
					files.push_back(p.assign(path).append("\\").append(fileinfo.name));
				}
				
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}
void Writeinput();//根据台站位置求最近点的岩性数据
int main(int argc, char** argv)
{
	int num_omp_threads=10;
	double m = 0.0;

	num_omp_threads = atoi(argv[1]);
	m= atof(argv[2]);
	omp_set_num_threads(num_omp_threads);
	cout << "开始计算..." << endl;
	Writeinput();
	cout << "读取台站岩性数据结束" << endl;
	ifstream ipt("stations.txt");
	string tempP;
	getline(ipt, tempP);
	std::string epCenterX, epCenterY;
	ipt >> epCenterX >> epCenterY;
	getline(ipt, tempP);

	int num_stations = 0;
	ipt >> num_stations;
	string temp;
	getline(ipt, temp);
	getline(ipt,temp);//skip the second line of input file
	station * st = new station[num_stations];

	Stat st_temp;
	system("rd check");
	system("md check");
	vector<string> files;
	getJustCurrentDir(".\\groundmotion", files);
	//vector<string> files1;
	//getFiles2GM(files[0],files1);
	ofstream outputStation("SlopeTo30cm.txt");
	outputStation << "Epicentre"<<"\n";
	outputStation<< epCenterX <<"\t"<< epCenterY << "\n";
	outputStation << num_stations<<"\n";
	outputStation << "Stationname\tx\ty\tslpoe\n";
	for (int i = 0; i < num_stations; i++)
	{
		ipt >> st[i].name >> st[i].x >> st[i].y >> st[i].C >> st[i].phi >>st[i].gamma >> st[i].slope >> st[i].lithoClass;
		cout << "Start computing the " <<st[i].name << " station" << endl;
		vector<string> files1;
		getFiles2GM(files[i], files1);		
		string fileLandslideDisp= st[i].name + "_LandslideDisp.txt";
		ofstream opt1(fileLandslideDisp.c_str());
		opt1 << "Slope\tNewmarkDisp(m)\tLandslideProbability\tVeryLow\tLow\tModerate\tHigh\n";
		if (st[i].lithoClass==0)
		{
			st[i].lithoClass = 1;
		}
		if (st[i].lithoClass>0)
		{
			GroundMotion gmEW;
			GroundMotion gmNS;
			gmEW.ReadGM(files1[0]);
			gmNS.ReadGM(files1[1]);
			st[i].getLithoC_phi();
			vector<double> slope;
			vector <double> Disp;
			int n = 1000;
			vector <double> C_array;
			vector <double> phi_array;
			st[i].GetPhi_C(n, C_array, phi_array);
			for (int ii = 10; ii < 11; ii++)
			{
				double cql = ii ;
				slope.push_back(ii);
				cout << "For slope "<< ii  << endl;
				vector<double> displacement;
				double *DispTemp;
				DispTemp = new double[n];//防止并行时用vector出错
				#pragma omp parallel for
				for (int j = 0; j < n; j++)  //Monte Carlo
				{
					cout << "Monte Carlo simulation " << j << endl;
					double Ac = st[i].getAcBy(C_array[j], phi_array[j], 25,m);
				/*	cout << Ac << endl;*/
					double dispEW = NewmarkDisp(gmEW.t, gmEW.a, Ac);
					double dispNS = NewmarkDisp(gmNS.t, gmNS.a, Ac);
					if (dispEW < dispNS)
					{
						dispEW = dispNS;
					}
					DispTemp[j] = dispEW;
					//displacement.push_back(dispEW);
				}
				for (int j = 0; j < n; j++)  
				{
					displacement.push_back(DispTemp[j]);
				}
				vector <double>LandslideRatio;
				st_temp.GetLandslideRatio(displacement, LandslideRatio);
				double tempCql = st_temp.getMedian(displacement);
				opt1 << cql << "\t" << tempCql << "\t" << 0.335*(1 - exp(-0.048*pow(st_temp.getMedian(displacement)*100.0, 1.565))) << "\t"
							<< LandslideRatio[0] << "\t" << LandslideRatio[1] << "\t" << LandslideRatio[2] << "\t" << LandslideRatio[3] << "\n";
				Disp.push_back(tempCql);
			}
			for (int jj=0;jj<80;jj++)
			{
				if (Disp[jj]>0.3)
				{
					outputStation << st[i].name <<"\t"<< st[i].x << "\t" << st[i].y << "\t" << slope[jj] << "\n";
					break;
				}
				if (jj==79)
				{
					outputStation << st[i].name << "\t" << st[i].x << "\t" << st[i].y << "\t" << 90 << "\n";
				}
			}
		}

	}
	outputStation.close();
	cout << "Finish!" << endl;
	cout << "Start to generate shp file" << endl;

	return 0;
}


double NewmarkDisp(vector<double>&time, vector<double>&acc,double ay)
{
	double Displacement = 0.0;
	if (ay < 0)
	{
		Displacement = 10.0;

	}
	else {
		double v0, v1, d0, d1;
		d0 = d1 = v1 = v0 = 0.0;
		vector<double> velocity;
		double dt = time[1] - time[0];
		velocity.push_back(0.0);
		double *abs_acc;
		abs_acc = new double[time.size()];
		double *rel_acc;
		rel_acc = new double[time.size()];
		double *rel_vel;
		rel_vel = new double[time.size()];
		double *rel_dis;
		rel_dis = new double[time.size()];

		abs_acc[0] = 0;		//% total acceleration
		rel_acc[0] = 0;		//% relative acceleration
		rel_vel[0] = 0;		//% relative velocity
		rel_dis[0] = 0;		//% relative displacement

		for (unsigned int i = 1; i < time.size(); i++)
		{
			abs_acc[i] = ay;
			rel_acc[i] = acc[i] - ay;
			rel_vel[i] = rel_vel[i - 1] + 0.5*(rel_acc[i - 1] + rel_acc[i])*dt;

			if (rel_vel[i] < 0)
			{
				abs_acc[i] = 0;
				rel_vel[i] = 0;
				rel_acc[i] = 0;
			}
			rel_dis[i] = rel_dis[i - 1] + rel_vel[i - 1] * dt + (2 * rel_acc[i - 1] + rel_acc[i])*dt*dt / 6;
		}

		Displacement = rel_dis[time.size() - 1];
		delete[]abs_acc;
		delete[]rel_acc;
		delete[]rel_vel;
		delete[]rel_dis;
	}
	
	return Displacement;
}
void getJustCurrentDir(string path, vector<string>& files)
{
	//文件句柄 
	long long  hFile = 0;
	//文件信息 
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib & _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					//files.push_back(fileinfo.name);
					files.push_back(p.assign(path).append("\\").append(fileinfo.name) );
				}

			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

wchar_t *multiByteToWideChar(const string& pKey)
{
	const char* pCStrKey = pKey.c_str();
	//第一次调用返回转换后的字符串长度，用于确认为wchar_t*开辟多大的内存空间
	int pSize = MultiByteToWideChar(CP_OEMCP, 0, pCStrKey, strlen(pCStrKey) + 1, NULL, 0);
	wchar_t *pWCStrKey = new wchar_t[pSize];
	//第二次调用将单字节字符串转换成双字节字符串
	MultiByteToWideChar(CP_OEMCP, 0, pCStrKey, strlen(pCStrKey) + 1, pWCStrKey, pSize);
	return pWCStrKey;
}

