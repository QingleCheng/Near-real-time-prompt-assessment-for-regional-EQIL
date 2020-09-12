#include "groundmotion.h"

GroundMotion::GroundMotion()
{

}

void GroundMotion::ReadGM(const string path)
{
    name=path;
    string s=path;
    ifstream fin(s.c_str());
    if(!fin.good())
    {
        cout<<"Failed to read ground motion file "<<s.c_str()<<"!\n";
        exit(1);
    }
    int nop=0;  //number of points
    fin>>nop;
    t.clear();
    a.clear();
    t.resize(nop);
    a.resize(nop);
    for(int i=0;i<nop;++i)
    {
        fin>>t[i]>>a[i];
    }
    fin.close();
}

double GroundMotion::GetTotalDuration()
{
    return t.back();
}

double GroundMotion::GetDt()
{
    if(t.size()>=2)
        return t[1]-t[0];

    return 0;
}

double GroundMotion::GetPGA()
{
	double PGA = 0.0;
	for (unsigned int i=0;i<t.size();i++)
	{
		if (PGA<abs(a[i]))
		{
			PGA = abs(a[i]);
		}
	}
	return PGA;
}

double GroundMotion::GetIrias()
{
	double Ia = 0.0;
	for (unsigned int i=0;i<t.size()-1;i++)
	{
		Ia = Ia + (a[i] * a[i] + a[i + 1] * a[i + 1]) / 2.0*GetDt();
	}
	return Ia*PI/2.0/9.8;
}

double GroundMotion::GetPGV()
{
	double PGV=0.0;
	double temp=0.0;
	double dt = GetDt();
	v.push_back(0.0);
	for (unsigned int i=1;i<t.size();i++)
	{
		temp = temp + (a[i-1]+a[i])*dt/2.0;
		v.push_back(temp);
		if (PGV<abs(temp))
		{
			PGV = abs(temp);
		}
	}
	return PGV;

}