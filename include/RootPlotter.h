#ifndef ROOTPLOTTER_H
#define ROOTPLOTTER_H

#include <vector>
#include <string>

#include "Nucleus.h"
#include "RxnType.h"
#include "MaskFile.h"

#include <THashTable.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>

struct GraphData {
	std::string name;
	std::string title;
	std::vector<double> xvec;
	std::vector<double> yvec;
	int color;
};

class RootPlotter {
public:
	RootPlotter();
	~RootPlotter();
	inline void ClearTable() { table->Clear(); };
	inline THashTable* GetTable() {
		GenerateGraphs(); 
		return table; 
	};
	void FillData(const Mask::Nucleus& nuc, double detKE = 0.0, const std::string& modifier = "");
	void FillCorrelations(const Mask::MaskFileData& data, Mask::RxnType type);
	void FillCorrelationsDetected(const Mask::MaskFileData& data, Mask::RxnType type);

private:
	THashTable* table;

	void GenerateGraphs();
	void MyFill(const std::string& name, const std::string& title, int bins, float min, float max, double val);
	void MyFill(const std::string& name, const std::string& title, int binsx, float minx, float maxx, int binsy, float miny, float maxy, double valx, double valy);
	void MyFill(const std::string& name, const std::string& title, double valx, double valy, int color); //TGraph

	std::vector<TGraph*> garbage_collection;
	std::vector<GraphData> graphs;

	static constexpr double rad2deg = 180.0/M_PI;

};

#endif