#include "Plotter.h"
#include <iostream>

namespace Mask {

Plotter::Plotter() :
	table(new THashTable())
{
}

Plotter::~Plotter() {
	for(unsigned int i=0; i<garbage_collection.size(); i++) {
		delete garbage_collection[i];
	}
	garbage_collection.clear();
	delete table;
}

void Plotter::FillData(const Nucleus& nuc) {
	std::string sym = nuc.GetIsotopicSymbol();
	std::string ke_vs_th_name = sym + "_ke_vs_theta";
	std::string ke_vs_th_title = ke_vs_th_name + ";#theta_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ke_vs_ph_name = sym + "_ke_vs_phi";
	std::string ke_vs_ph_title = ke_vs_ph_name + ";#phi_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ex_name = sym + "_ex";
	std::string ex_title = ex_name + ";E_{ex} (MeV);counts";
	
	MyFill(ke_vs_th_name.c_str(), ke_vs_th_title.c_str(), nuc.GetTheta()*rad2deg, nuc.GetKE(), 2);
	MyFill(ke_vs_ph_name.c_str(), ke_vs_ph_title.c_str(), nuc.GetPhi()*rad2deg, nuc.GetKE(), 4);
	MyFill(ex_name.c_str(),ex_title.c_str(),260,-1.0,25,nuc.GetExcitationEnergy());
	
}

void Plotter::MyFill(const char* name, const char* title, int bins, float min, float max, double val) {
	TH1F* h = (TH1F*) table->FindObject(name);
	if(h) {
		h->Fill(val);
	} else {
		h = new TH1F(name, title, bins, min, max);
		h->Fill(val);
		table->Add(h);
	}
}

void Plotter::MyFill(const char* name, const char* title, int binsx, float minx, float maxx, int binsy, float miny, float maxy, double valx, double valy) {
	TH2F* h = (TH2F*) table->FindObject(name);
	if(h) {
		h->Fill(valx, valy);
	} else {
		h = new TH2F(name, title, binsx, minx, maxx, binsy, miny, maxy);
		h->Fill(valx, valy);
		table->Add(h);
	}
}

void Plotter::MyFill(const char* name, const char* title, double valx, double valy, int color) {
	for(auto& g : graphs) {
		if(g.name == name) {
			g.xvec.push_back(valx);
			g.yvec.push_back(valy);
			return;
		}
	}

	GraphData new_g;
	new_g.name = name;
	new_g.title = title;
	new_g.xvec.push_back(valx);
	new_g.yvec.push_back(valy);
	new_g.color = color;

	graphs.push_back(new_g);
}

void Plotter::GenerateGraphs() {
	for(auto& g : graphs) {
		TGraph* graph = new TGraph(g.xvec.size(), &(g.xvec[0]), &(g.yvec[0]));
		graph->SetName(g.name.c_str());
		graph->SetTitle(g.title.c_str());
		graph->SetMarkerColor(g.color);
		table->Add(graph);
		garbage_collection.push_back(graph);
	}
}

};