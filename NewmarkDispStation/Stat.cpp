#include "Stat.h"

Stat::Stat()
{

}

void Stat::SetSeed(int seed)
{
    _seed=seed;
    if (_seed==0)
    {
        srand((int)time(0));
    }
    else
    {
        srand(_seed);
    }
}


double Stat::random(double min, double max)
{
    double r=(double)rand() / (double)RAND_MAX;
    r=(max-min)*r+min;
    return r;
}

int Stat::random(int min, int max)
{
    return (rand() % (max-min+1))+ min;
}

int Stat::round(double r)
{
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

double Stat::CDF_normal(double x,double median,double std)
{
    if (std==0)
    {
        std=0.0000001;
    }
    x=(x-median)/std;	//convert into standard normal distribution
    int index=round(fabs(x)*100.0);
    double y=0.0;
    if(index<0)
    {
        std::cout<<"Error calculating CDF!\n";
        exit(2);
    }

    if (index>=390)
    {
        y=1.0000;
    }
    else if (index>299)	//Look up the other standard normal distribution table
    {
        index=round(index/10.0-30.0);
        y=StandardNormalDistributionTable2[index];
    }
    else	//Look up the standard normal distribution table
    {
        y=StandardNormalDistributionTable1[index];
    }

    if (x<0)
    {
        y=1.0-y;
    }
    return y;
}

double Stat::gaussrand(double median, double std)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    X = X * std + median;

    return X;

}



void Stat::GetLandslideRatio(vector<double> disp,vector<double>&Ratio)
{
	int countVeryLow,conutLow,countModerate,countHigh;
	countVeryLow = conutLow = countModerate = countHigh = 0;
	double veryLow = 0.05;
	double low = 0.15;
	double moderate = 0.3;
	for (unsigned int i=0;i<disp.size();i++)
	{
		if (disp[i]<veryLow)
		{
			countVeryLow++;
		} 
		else if(disp[i]<low && disp[i]>veryLow )
		{
			conutLow++;
		}
		else if (disp[i]<moderate && disp[i]>low )
		{
			countModerate++;
		}
		else {
			countHigh++;
		}
	}
	double size_disp = disp.size();
	Ratio.push_back(countVeryLow / size_disp);
	Ratio.push_back(conutLow / size_disp);
	Ratio.push_back(countModerate / size_disp);
	Ratio.push_back(countHigh / size_disp);
}