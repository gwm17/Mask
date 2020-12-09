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

 	G3Vec coords;
 	for(int i=0; i<5; i++) {
 		for(int j=0; j<detectors[i].GetNumberOfRings(); j++) {
 			for(int k=0; k<4; k++) {
 				coords = detectors[i].GetRingTiltCoords(j, k);
 				ringxs.push_back(coords.GetX());
 				ringys.push_back(coords.GetY());
 				ringzs.push_back(coords.GetZ());
 			}
 		}
 	}

 	for(int i=0; i<5; i++) {
 		for(int j=0; j<detectors[i].GetNumberOfWedges(); j++) {
 			for(int k=0; k<4; k++) {
 				coords = detectors[i].GetWedgeTiltCoords(j, k);
 				wedgexs.push_back(coords.GetX());
 				wedgeys.push_back(coords.GetY());
 				wedgezs.push_back(coords.GetZ());
 			}
 		}
 	}
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
	std::vector<double> resid_xs, eject_xs;
	std::vector<double> resid_ys, eject_ys;
	std::vector<double> resid_zs, eject_zs;

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	G3Vec coordinates;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(eject->KE >= ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				coordinates = det.GetTrajectoryCoordinates(eject->theta, eject->phi);
				if(coordinates.GetX() != 0) {
					eject_xs.push_back(coordinates.GetX());
					eject_ys.push_back(coordinates.GetY());
					eject_zs.push_back(coordinates.GetZ());
					break;
				}
			}
		}

		if(resid->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.GetTrajectoryCoordinates(resid->theta, resid->phi).GetX() != 0) {
					resid_xs.push_back(coordinates.GetX());
					resid_ys.push_back(coordinates.GetY());
					resid_zs.push_back(coordinates.GetZ());
					break;
				}
			}
		}

	}

	double ejecteff = ((double) eject_xs.size())/nevents;
	double resideff = ((double) resid_xs.size())/nevents;
	TParameter<double> eject_eff("Light Breakup Efficiency", ejecteff);
	TParameter<double> resid_eff("Heavy Breakup Efficiency", resideff);

	TGraph2D* gde = new TGraph2D(eject_xs.size(), &(eject_xs[0]), &(eject_ys[0]), &(eject_zs[0]));
	gde->SetName("detected_eject_points");
	gde->SetMarkerStyle(2);
	gde->SetMarkerColor(2);

	TGraph2D* gr = new TGraph2D(ringxs.size(), &(ringxs[0]), &(ringys[0]), &(ringzs[0]));
	gr->SetName("ring_detector_edges");
	gr->SetTitle("SABRE Detector; x(m); y(m); z(m)");
	gr->SetMarkerStyle(2);

	TGraph2D* gw = new TGraph2D(wedgexs.size(), &(wedgexs[0]), &(wedgeys[0]), &(wedgezs[0]));
	gw->SetName("wedge_detector_edges");
	gw->SetTitle("SABRE Detector Wedges; x(m); y(m); z(m)");
	gw->SetMarkerStyle(2);

	TCanvas* canvas = new TCanvas();
	canvas->SetName("detectors_and_particles");
	canvas->cd();
	gr->Draw("AP");
	gw->Draw("same P");
	gde->Draw("same P");

	input->cd();
	eject_eff.Write();
	resid_eff.Write();
	gr->Write();
	gw->Write();
	gde->Write();
	canvas->Write();
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
				if(det.GetTrajectoryCoordinates(break1->theta, break1->phi).GetX() != 0) {
					b1_thetas.push_back(break1->theta);
					b1_phis.push_back(break1->phi);
					b1_kes.push_back(break1->KE);
					break;
				}
			}
		}

		if(break2->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.GetTrajectoryCoordinates(break2->theta, break2->phi).GetX() != 0) {
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
				if(det.GetTrajectoryCoordinates(break1->theta, break1->phi).GetX() != 0) {
					b1_thetas.push_back(break1->theta);
					b1_phis.push_back(break1->phi);
					b1_kes.push_back(break1->KE);
					break;
				}
			}
		}


		if(break3->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.GetTrajectoryCoordinates(break3->theta, break3->phi).GetX() != 0) {
					b3_thetas.push_back(break3->theta);
					b3_phis.push_back(break3->phi);
					b3_kes.push_back(break3->KE);
					break;
				}
			}
		}

		if(break4->KE > ENERGY_THRESHOLD) {
			for(auto& det : detectors) {
				if(det.GetTrajectoryCoordinates(break4->theta, break4->phi).GetX() != 0) {
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