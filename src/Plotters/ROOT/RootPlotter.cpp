#include "RootPlotter.h"
#include "MaskFile.h"
#include <TFile.h>

#include <iostream>

RootPlotter::RootPlotter() :
	table(new THashTable())
{
}

RootPlotter::~RootPlotter() {}

void RootPlotter::FillData(const Mask::Nucleus& nuc, const std::string& modifier) {
	std::string sym = nuc.GetIsotopicSymbol();
	std::string ke_vs_th_name = sym + modifier + "_ke_vs_theta";
	std::string ke_vs_th_title = ke_vs_th_name + ";#theta_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ke_vs_ph_name = sym + modifier + "_ke_vs_phi";
	std::string ke_vs_ph_title = ke_vs_ph_name + ";#phi_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ex_name = sym + modifier + "_ex";
	std::string ex_title = ex_name + ";E_{ex} (MeV);counts";
	std::string angdist_name = sym + modifier +"_angDist";
	std::string angdist_title = angdist_name+";cos#right(#theta_{CM}#left);counts";
	
	MyFill(ke_vs_th_name.c_str(), ke_vs_th_title.c_str(), nuc.GetTheta()*rad2deg, nuc.GetKE(), 2);
	MyFill(ke_vs_ph_name.c_str(), ke_vs_ph_title.c_str(), nuc.GetPhi()*rad2deg, nuc.GetKE(), 4);
	MyFill(ex_name.c_str(),ex_title.c_str(),260,-1.0,25,nuc.GetExcitationEnergy());
	MyFill(angdist_name.c_str(), angdist_title.c_str(),100,-1.0,1.0,std::cos(nuc.GetThetaCM()));
	
}

void RootPlotter::MyFill(const std::string& name, const std::string& title, int bins, float min, float max, double val) {
	TH1F* h = (TH1F*) table->FindObject(name.c_str());
	if(h) {
		h->Fill(val);
	} else {
		h = new TH1F(name.c_str(), title.c_str(), bins, min, max);
		h->Fill(val);
		table->Add(h);
	}
}

void RootPlotter::MyFill(const std::string& name, const std::string& title, int binsx, float minx, float maxx, int binsy, float miny, float maxy, double valx, double valy) {
	TH2F* h = (TH2F*) table->FindObject(name.c_str());
	if(h) {
		h->Fill(valx, valy);
	} else {
		h = new TH2F(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
		h->Fill(valx, valy);
		table->Add(h);
	}
}

void RootPlotter::MyFill(const std::string& name, const std::string& title, double valx, double valy, int color) {
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

void RootPlotter::GenerateGraphs() {
	for(auto& g : graphs) {
		TGraph* graph = new TGraph(g.xvec.size(), &(g.xvec[0]), &(g.yvec[0]));
		graph->SetName(g.name.c_str());
		graph->SetTitle(g.title.c_str());
		graph->SetMarkerColor(g.color);
		table->Add(graph);
		garbage_collection.push_back(graph);
	}
}

int main(int argc, char** argv) {
	if(argc != 3) {
		std::cout<<"Unable to run ROOT plotting tool with incorrect number of arguments! Expected 2 args, given: "<<argc<<" Exiting."<<std::endl;
		return 1;
	}

	std::string inputname = argv[1];
	std::string outputname = argv[2];

	Mask::MaskFile input(inputname, Mask::MaskFile::FileType::read);

	TFile* root_out = TFile::Open(outputname.c_str(), "RECREATE");

	RootPlotter plotter;

	Mask::MaskFileHeader header = input.ReadHeader();

	std::cout<<"File Header -- rxn type: "<<header.rxn_type<<" nsamples: "<<header.nsamples<<std::endl;

	Mask::MaskFileData data;
	Mask::Nucleus nucleus;

	double flush_frac = 0.05;
	int count = 0, flush_val = flush_frac*header.nsamples, flush_count = 0;
	while(true) {
		if(count == flush_val) {
			count = 0;
			flush_count++;
			std::cout<<"\rPercent of file processed: "<<flush_frac*flush_count*100<<"%"<<std::flush;
		}

		data = input.ReadData();
		if(data.eof)
			break;
		for(unsigned int i=0; i<data.Z.size(); i++) {
			nucleus.SetIsotope(data.Z[i], data.A[i]);
			nucleus.SetVectorSpherical(data.theta[i], data.phi[i], data.p[i], data.E[i]);
			plotter.FillData(nucleus);
			if(data.detect_flag[i] == true) {
				plotter.FillData(nucleus, "detected");
			}
		}
		count++;
	}
	std::cout<<std::endl;

	input.Close();

	root_out->cd();
	plotter.GetTable()->Write();
	root_out->Close();

	return 0;
}