/********************************************************
* @file    : main.cpp
* @brief   : 	
* @details : 	
* @author  : Chinler
* @date    : 2019-2-28
*********************************************************/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>
#include <jansson.h>     // for writing json
#include <nanoflann.hpp> // for searching for nearest point

using namespace nanoflann;

struct locations {
	locations() :x(0), y(0), LithoGroup(0){}
	locations(std::string st, double a, double b,int c) :station(st), x(a), y(b), LithoGroup(c){}
	std::string station;
	double x;
	double y;
	int LithoGroup;
};


template <typename T>
struct PointCloud
{
	struct Point
	{
		Point() : stationTag(0), x(0.), y(0.) {}
		Point(int tag, T(a), T(b)) : stationTag(tag), x(a), y(b) {}
		int stationTag;
		T  x, y;
	};

	std::vector<Point>  pts;

	inline size_t kdtree_get_point_count() const { return pts.size(); }

	inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t /*size*/) const
	{
		const T d0 = p1[0] - pts[idx_p2].x;
		const T d1 = p1[1] - pts[idx_p2].y;
		return d0*d0 + d1*d1;
	}

	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) return pts[idx].x;
		else return pts[idx].y;
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};
int CountLine(std::string filename)
{
	std::ifstream ipt(filename);
	char c;
	int lineCnt = 0;
	while (ipt.get(c))
	{
		if (c == '\n')
			lineCnt++;
	}
	return lineCnt+1;
}
void Writeinput() {
	std::map<int, locations> stationLocations;
	PointCloud<float> cloud;
	//
	// first parse the station file & put each station into the cloud of points
	//
	std::string inputPoint = ".\\NationalPoint.txt";
	std::ifstream ipt(inputPoint.c_str());
	int count = 0;
	int NumLine = 1844032;//CountLine(inputPoint);节省时间
	std::string temp;
	getline(ipt, temp);
	for (int i=0;i<NumLine;i++)
	{
		std::string s;
		ipt >> s;
		std::string station(s);
		double x, y;
		int LithoGroup;
		ipt >> LithoGroup>>x >> y;
		stationLocations[count] = locations(station, x, y, LithoGroup);
		cloud.pts.resize(count + 1);
		cloud.pts[count].stationTag = count;
		cloud.pts[count].x = x;
		cloud.pts[count].y = y;
		count++;
	}
	//
	// now find nearest point in the cloud
	//

	// build the kd tree
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, PointCloud<float> >,
		PointCloud<float>,
		2
	> my_kd_tree_t;

	my_kd_tree_t   index(2, cloud, KDTreeSingleIndexAdaptorParams(10));
	index.buildIndex();
	//
	// now parse the point file for the location and
	//
	std::string inputStation = ".\\StationList.txt";
	int LineStation = CountLine(inputStation)-2;
	std::ifstream ipt1(inputStation.c_str());
	getline(ipt1, temp);
	std::string epCenterX, epCenterY;
	ipt1 >> epCenterX >> epCenterY;
	getline(ipt1, temp);
	std::ofstream opt("stations.txt");
	opt << "震中" << "\n";
	opt << epCenterX <<"\t"<< epCenterY << "\n";
	opt << LineStation << "\n";
	opt << "Name\tx\ty\tC\tphi\tgamma\tslope\tlithoClass\n";
	for (int i=0;i<LineStation;i++)
	{
		float buildingLoc[2];
		int ID;
		std::string StationName;
		ipt1 >> ID >> StationName >> buildingLoc[0] >> buildingLoc[1];
		//
		// do a knn search to find nearest point
		//
		long unsigned int num_results = 1;
		long long unsigned int ret_index;
		//unsigned int ret_index;
		float out_dist_sqr;
		nanoflann::KNNResultSet<float> resultSet(num_results);
		resultSet.init(&ret_index, &out_dist_sqr);
		index.findNeighbors(resultSet, &buildingLoc[0], nanoflann::SearchParams(10));
		//
		// create the event
		//
		int stationTag = ret_index;
		std::map<int, locations>::iterator stationIter;
		stationIter = stationLocations.find(stationTag);
		if (stationIter != stationLocations.end()) {
			std::cerr << stationIter->second.station;
			//stationName = "EQX" + stationIter->second.station + "\n";
			opt << StationName<<"\t"<< buildingLoc[0]<<"\t"<< buildingLoc[1] << "\t" << "-1\t-1\t-1\t-1"<<"\t"
				<<stationIter->second.LithoGroup <<"\n";
		}
	}
}