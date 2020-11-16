#ifndef PLOTTER_H
#define PLOTTER_H

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <THashTable.h>
#include "Nucleus.h"

struct GraphData {
	std::string name;
	std::string title;
	std::vector<double> xvec;
	std::vector<double> yvec;
	int color;
};

class Plotter {
public:
	Plotter();
	~Plotter();
	inline void ClearTable() { table->Clear(); };
	inline THashTable* GetTable() {
		GenerateGraphs(); 
		return table; 
	};
	void FillData(const Nucleus& nuc);

private:
	THashTable* table;

	void GenerateGraphs();
	void MyFill(const char* name, const char* title, int bins, float min, float max, double val);
	void MyFill(const char* name, const char* title, int binsx, float minx, float maxx, int binsy, float miny, float maxy, double valx, double valy);
	void MyFill(const char* name, const char* title, double valx, double valy, int color); //TGraph

	std::vector<TGraph*> garbage_collection;
	std::vector<GraphData> graphs;

	static constexpr double rad2deg = 180.0/M_PI;

};

#endif