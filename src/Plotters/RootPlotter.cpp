#include "RootPlotter.h"
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <iostream>

RootPlotter::RootPlotter()
{
	TH1::AddDirectory(kFALSE);
	//Enforce dictionary linking
	if(Mask::EnforceDictionaryLinked())
	{
		std::cout<<"Dictionary Linked"<<std::endl;
	}
}

RootPlotter::~RootPlotter() {}

void RootPlotter::Run(const std::string& inputname, const std::string& outputname)
{
	TFile* input = TFile::Open(inputname.c_str(), "READ");
	TTree* tree = (TTree*) input->Get("SimTree");
	std::vector<Mask::Nucleus>* dataHandle = new std::vector<Mask::Nucleus>();
	tree->SetBranchAddress("nuclei", &dataHandle);

	TFile* output = TFile::Open(outputname.c_str(), "RECREATE");

	double flushFrac = 0.05;
	uint64_t nentries = tree->GetEntries();
	uint64_t flushVal = flushFrac*nentries;
	uint64_t count=0;
	uint64_t flushCount = 0;

	for(uint64_t i=0; i<nentries; i++)
	{
		tree->GetEntry(i);
		for(Mask::Nucleus& nuc : *(dataHandle))
		{
			FillData(nuc);
		}
	}

	input->Close();
	delete dataHandle;

	output->cd();
	for(auto& obj : m_map)
		obj.second->Write(obj.second->GetName(), TObject::kOverwrite);
	output->Close();
}

void RootPlotter::FillData(const Mask::Nucleus& nuc)
{
	std::string modifier = "";
	if(nuc.isDetected)
		modifier = "detected";

	std::string sym = nuc.isotopicSymbol;
	std::string ke_vs_th_name = sym + modifier + "_ke_vs_theta";
	std::string ke_vs_th_title = ke_vs_th_name + ";#theta_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ke_vs_ph_name = sym + modifier + "_ke_vs_phi";
	std::string ke_vs_ph_title = ke_vs_ph_name + ";#phi_{lab} (degrees);Kinetic Energy (MeV)";
	std::string th_vs_ph_name = sym + modifier + "_theta_vs_phi";
	std::string th_vs_ph_title = th_vs_ph_name + ";#theta_{lab};#phi_{lab}";
	std::string ex_name = sym + modifier + "_ex";
	std::string ex_title = ex_name + ";E_{ex} (MeV);counts";
	std::string angdist_name = sym + modifier +"_angDist";
	std::string angdist_title = angdist_name+";cos#right(#theta_{CM}#left);counts";
	
	if(nuc.isDetected)
	{
		MyFill(ke_vs_th_name.c_str(), ke_vs_th_title.c_str(), nuc.vec4.Theta()*s_rad2deg, nuc.GetKE(), 2);
		MyFill(ke_vs_ph_name.c_str(), ke_vs_ph_title.c_str(), nuc.vec4.Phi()*s_rad2deg, nuc.GetKE(), 4);
		MyFill(th_vs_ph_name.c_str(), th_vs_ph_title.c_str(), nuc.vec4.Theta()*s_rad2deg, nuc.vec4.Phi()*s_rad2deg, 2);
		MyFill(ex_name.c_str(),ex_title.c_str(),260,-1.0,25,nuc.GetExcitationEnergy());
		MyFill(angdist_name.c_str(), angdist_title.c_str(),20,-1.0,1.0,std::cos(nuc.thetaCM));
	}
	else
	{
		MyFill(ke_vs_th_name.c_str(), ke_vs_th_title.c_str(), nuc.vec4.Theta()*s_rad2deg, nuc.detectedKE, 2);
		MyFill(ke_vs_ph_name.c_str(), ke_vs_ph_title.c_str(), nuc.vec4.Phi()*s_rad2deg, nuc.detectedKE, 4);
		MyFill(th_vs_ph_name.c_str(), th_vs_ph_title.c_str(), nuc.vec4.Theta()*s_rad2deg, nuc.vec4.Phi()*s_rad2deg, 2);
		MyFill(ex_name.c_str(),ex_title.c_str(),260,-1.0,25,nuc.GetExcitationEnergy());
		MyFill(angdist_name.c_str(), angdist_title.c_str(),20,-1.0,1.0,std::cos(nuc.thetaCM));
	}
	
}

void RootPlotter::MyFill(const std::string& name, const std::string& title, int bins, float min, float max, double val)
{
	auto iter = m_map.find(name);
	if(iter != m_map.end())
	{
		std::shared_ptr<TH1> h = std::static_pointer_cast<TH1>(iter->second);
		h->Fill(val);
	}
	else
	{
		std::shared_ptr<TH1F> h = std::make_shared<TH1F>(name.c_str(), title.c_str(), bins, min, max);
		h->Fill(val);
		m_map[name] = h;
	}
}

void RootPlotter::MyFill(const std::string& name, const std::string& title, int binsx, float minx, float maxx,
						 int binsy, float miny, float maxy, double valx, double valy)
{
	auto iter = m_map.find(name);
	if(iter != m_map.end())
	{
		std::shared_ptr<TH2> h = std::static_pointer_cast<TH2>(iter->second);
		h->Fill(valx, valy);
	}
	else
	{
		std::shared_ptr<TH2F> h = std::make_shared<TH2F>(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
		h->Fill(valx, valy);
		m_map[name] = h;
	}
}

void RootPlotter::MyFill(const std::string& name, const std::string& title, double valx, double valy, int color)
{
	auto iter = m_map.find(name);
	if(iter != m_map.end())
	{
		std::shared_ptr<TGraph> g = std::static_pointer_cast<TGraph>(iter->second);
		g->SetPoint(g->GetN(), valx, valy);
	}
	else
	{
		std::shared_ptr<TGraph> g = std::make_shared<TGraph>(1, &valx, &valy);
		g->SetName(name.c_str());
		g->SetTitle(title.c_str());
		g->SetMarkerColor(color);
		m_map[name] = g;
	}
}