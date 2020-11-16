#include "SabreEfficiency.h"
#include "Kinematics.h"
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>

SabreEfficiency::SabreEfficiency() : 
	m_rxn_type(-1)
{
	detectors.reserve(5);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI0*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI1*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI2*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI3*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI4*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
}

SabreEfficiency::~SabreEfficiency() {}

void SabreEfficiency::CalculateEfficiency(const char* file) {
	std::cout<<"----------SABRE Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<file<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

	switch(m_rxn_type) {
		case Kinematics::ONESTEP_DECAY:
		{
			RunDecay(file);
			break;
		}
		case Kinematics::TWOSTEP:
		{
			Run2Step(file);
			break;
		}
		case Kinematics::THREESTEP:
		{
			Run3Step(file);
			break;
		}
	}
	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}

void SabreEfficiency::RunDecay(const char* file) {
	TFile* input = TFile::Open(file, "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	NucData* eject = new NucData();
	NucData* resid = new NucData();


	tree->SetBranchAddress("ejectile", &eject);
	tree->SetBranchAddress("residual", &resid);

	double nevents = tree->GetEntries();
	std::vector<double> resid_thetas, eject_thetas;
	std::vector<double> resid_phis, eject_phis;
	std::vector<double> resid_kes, eject_kes;

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(eject->KE >= ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(eject->theta, eject->phi)) {
					eject_thetas.push_back(eject->theta);
					eject_phis.push_back(eject->phi);
					eject_kes.push_back(eject->KE);
					break;
				}
			}
		}

		if(resid->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(resid->theta, resid->phi)) {
					resid_thetas.push_back(resid->theta);
					resid_phis.push_back(resid->phi);
					resid_kes.push_back(resid->KE);
					break;
				}
			}
		}

	}

	double ejecteff = ((double) eject_thetas.size())/nevents;
	double resideff = ((double) resid_thetas.size())/nevents;
	TParameter<double> eject_eff("Light Breakup Efficiency", ejecteff);
	TParameter<double> resid_eff("Heavy Breakup Efficiency", resideff);

	input->cd();
	eject_eff.Write();
	resid_eff.Write();
	input->Close();
}

void SabreEfficiency::Run2Step(const char* file) {

	TFile* input = TFile::Open(file, "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	NucData* break1 = new NucData();
	NucData* break2 = new NucData();


	tree->SetBranchAddress("breakup1", &break1);
	tree->SetBranchAddress("breakup2", &break2);

	double nevents = tree->GetEntries();
	std::vector<double> b1_thetas, b2_thetas;
	std::vector<double> b1_phis, b2_phis;
	std::vector<double> b1_kes, b2_kes;

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(break1->KE >= ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(break1->theta, break1->phi)) {
					b1_thetas.push_back(break1->theta);
					b1_phis.push_back(break1->phi);
					b1_kes.push_back(break1->KE);
					break;
				}
			}
		}

		if(break2->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(break2->theta, break2->phi)) {
					b2_thetas.push_back(break2->theta);
					b2_phis.push_back(break2->phi);
					b2_kes.push_back(break2->KE);
					break;
				}
			}
		}

	}

	double b1eff = ((double) b1_thetas.size())/nevents;
	double b2eff = ((double) b2_thetas.size())/nevents;
	TParameter<double> break1_eff("Light Breakup Efficiency", b1eff);
	TParameter<double> break2_eff("Heavy Breakup Efficiency", b2eff);

	input->cd();
	break1_eff.Write();
	break2_eff.Write();
	input->Close();

}

void SabreEfficiency::Run3Step(const char* file) {
	TFile* input = TFile::Open(file, "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	NucData* break1 = new NucData();
	NucData* break3 = new NucData();
	NucData* break4 = new NucData();


	tree->SetBranchAddress("breakup1", &break1);
	tree->SetBranchAddress("breakup3", &break3);
	tree->SetBranchAddress("breakup4", &break4);

	double nevents = tree->GetEntries();
	std::vector<double> b1_thetas, b3_thetas, b4_thetas;
	std::vector<double> b1_phis, b3_phis, b4_phis;
	std::vector<double> b1_kes, b3_kes, b4_kes;

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(break1->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(break1->theta, break1->phi)) {
					b1_thetas.push_back(break1->theta);
					b1_phis.push_back(break1->phi);
					b1_kes.push_back(break1->KE);
					break;
				}
			}
		}


		if(break3->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(break3->theta, break3->phi)) {
					b3_thetas.push_back(break3->theta);
					b3_phis.push_back(break3->phi);
					b3_kes.push_back(break3->KE);
					break;
				}
			}
		}

		if(break4->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.IsInside(break4->theta, break4->phi)) {
					b4_thetas.push_back(break4->theta);
					b4_phis.push_back(break4->phi);
					b4_kes.push_back(break4->KE);
					break;
				}
			}
		}
	}

	double b1eff = ((double) b1_thetas.size())/nevents;
	double b3eff = ((double) b3_thetas.size())/nevents;
	double b4eff = ((double) b4_thetas.size())/nevents;
	TParameter<double> break1_eff("Light Initial Breakup Efficiency", b1eff);
	TParameter<double> break3_eff("Light Final Breakup Efficiency", b3eff);
	TParameter<double> break4_eff("Heavy Final Breakup Efficiency", b4eff);

	input->cd();
	break1_eff.Write();
	break3_eff.Write();
	break4_eff.Write();
	input->Close();
}