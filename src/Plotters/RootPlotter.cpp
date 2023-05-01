#include "RootPlotter.h"
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <iostream>

static double FullPhi(double phi)
{
	return phi < 0.0 ? 2.0*M_PI + phi : phi; 
}

RootPlotter::RootPlotter()
{
	TH1::AddDirectory(kFALSE);
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
		count++;
		if(count == flushVal)
		{
			count = 0;
			flushCount++;
			std::cout<<"\rPercent of data processed: "<<flushCount*flushFrac*100<<"%"<<std::flush;
		}
		tree->GetEntry(i);
		for(Mask::Nucleus& nuc : *(dataHandle))
		{
			FillData(nuc);
		}
	}
	std::cout<<std::endl;
	input->Close();
	delete dataHandle;

	output->cd();
	for(auto& obj : m_map)
		obj.second->Write();
	output->Close();
}

void RootPlotter::FillData(const Mask::Nucleus& nuc)
{
	std::string mod = "detected";

	std::string sym = nuc.isotopicSymbol;
	std::string ke_vs_th_name = sym + "_ke_vs_theta";
	std::string ke_vs_th_title = ke_vs_th_name + ";#theta_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ke_vs_th_name_det = sym + mod + "_ke_vs_theta";
	std::string ke_vs_th_title_det = ke_vs_th_name + ";#theta_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ke_vs_ph_name = sym + "_ke_vs_phi";
	std::string ke_vs_ph_title = ke_vs_ph_name + ";#phi_{lab} (degrees);Kinetic Energy (MeV)";
	std::string ke_vs_ph_name_det = sym + mod + "_ke_vs_phi";
	std::string ke_vs_ph_title_det = ke_vs_ph_name + ";#phi_{lab} (degrees);Kinetic Energy (MeV)";
	std::string th_vs_ph_name = sym + "_theta_vs_phi";
	std::string th_vs_ph_title = th_vs_ph_name + ";#theta_{lab};#phi_{lab}";
	std::string th_vs_ph_name_det = sym + mod + "_theta_vs_phi";
	std::string th_vs_ph_title_det = th_vs_ph_name + ";#theta_{lab};#phi_{lab}";
	std::string ex_name = sym + "_ex";
	std::string ex_title = ex_name + ";E_{ex} (MeV);counts";
	std::string ex_name_det = sym + mod + "_ex";
	std::string ex_title_det = ex_name + ";E_{ex} (MeV);counts";
	std::string angdist_name = sym +"_angDist";
	std::string angdist_title = angdist_name+";cos#right(#theta_{CM}#left);counts";
	std::string angdist_name_det = sym + mod +"_angDist";
	std::string angdist_title_det = angdist_name+";cos#right(#theta_{CM}#left);counts";
	std::string hist_ke_th_name = sym + "_hist_ke_vs_theta";
	std::string hist_ke_th_title = hist_ke_th_name + ";#theta_{lab};Kinetic Energy (MeV)";
	std::string hist_ke_th_name_det = sym + mod + "_hist_ke_vs_theta";
	std::string hist_ke_th_title_det = hist_ke_th_name + ";#theta_{lab};Kinetic Energy (MeV)";
	
	MyFill(ke_vs_th_name, ke_vs_th_title, nuc.vec4.Theta()*s_rad2deg, nuc.GetKE(), 2);
	MyFill(ke_vs_ph_name, ke_vs_ph_title, FullPhi(nuc.vec4.Phi())*s_rad2deg, nuc.GetKE(), 4);
	MyFill(th_vs_ph_name, th_vs_ph_title, nuc.vec4.Theta()*s_rad2deg, FullPhi(nuc.vec4.Phi())*s_rad2deg, 2);
	MyFill(ex_name, ex_title, 260, -1.0, 25, nuc.GetExcitationEnergy());
	MyFill(angdist_name, angdist_title, 20, -1.0, 1.0, std::cos(nuc.thetaCM));
	MyFill(hist_ke_th_name, hist_ke_th_title, 180, 0.0, 180.0, 400, 0.0, 20.0, nuc.vec4.Theta()*s_rad2deg, nuc.GetKE());
	if(nuc.isDetected && nuc.detectedKE > 0.25)
	{
		MyFill(ke_vs_th_name_det, ke_vs_th_title_det, nuc.vec4.Theta()*s_rad2deg, nuc.detectedKE, 2);
		MyFill(ke_vs_ph_name_det, ke_vs_ph_title_det, FullPhi(nuc.vec4.Phi())*s_rad2deg, nuc.detectedKE, 4);
		MyFill(th_vs_ph_name_det, th_vs_ph_title_det, nuc.vec4.Theta()*s_rad2deg, FullPhi(nuc.vec4.Phi())*s_rad2deg, 2);
		MyFill(ex_name_det, ex_title_det, 260, -1.0, 25, nuc.GetExcitationEnergy());
		MyFill(angdist_name_det, angdist_title_det, 20, -1.0, 1.0, std::cos(nuc.thetaCM));
		MyFill(hist_ke_th_name_det, hist_ke_th_title_det, 180, 0.0, 180.0, 400, 0.0, 20.0, nuc.vec4.Theta()*s_rad2deg, nuc.detectedKE);
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