#ifndef ROOTPLOTTER_H
#define ROOTPLOTTER_H

#include <vector>
#include <string>
#include <unordered_map>

#include "Nucleus.h"

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>

class RootPlotter
{
public:
	RootPlotter();
	~RootPlotter();
	
	void Run(const std::string& inputname, const std::string& outputname);
	
private:
	void FillData(const Mask::Nucleus& nuc);
	void MyFill(const std::string& name, const std::string& title, int bins, float min, float max, double val);
	void MyFill(const std::string& name, const std::string& title, int binsx, float minx, float maxx, int binsy, float miny, float maxy, 
				double valx, double valy);
	void MyFill(const std::string& name, const std::string& title, double valx, double valy, int color); //TGraph

	std::unordered_map<std::string, std::shared_ptr<TObject>> m_map;

	static constexpr double s_rad2deg = 180.0/M_PI;
};

#endif