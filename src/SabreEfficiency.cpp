#include "SabreEfficiency.h"
#include "Kinematics.h"
#include <fstream>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>

SabreEfficiency::SabreEfficiency() : 
	m_rxn_type(-1), deadlayer(DEADLAYER_THIN)
{
	detectors.reserve(5);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI0*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI1*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI2*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI3*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI4*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);

 	Mask::Vec3 coords;
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
	std::vector<int> dead_z = {14};
	std::vector<int> dead_a = {28};
	std::vector<int> dead_stoich = {1};
	deadlayer.SetElements(dead_z, dead_a, dead_stoich);
}

SabreEfficiency::~SabreEfficiency() {}

void SabreEfficiency::MyFill(THashTable* table, const char* name, const char* title, int bins, float min, float max, double val) {
	TH1F* h = (TH1F*) table->FindObject(name);
	if(h) {
		h->Fill(val);
	} else {
		h = new TH1F(name, title, bins, min, max);
		h->Fill(val);
		table->Add(h);
	}
}

void SabreEfficiency::MyFill(THashTable* table, const char* name, const char* title, int binsx, float minx, float maxx, int binsy, float miny, float maxy, double valx, double valy) {
	TH2F* h = (TH2F*) table->FindObject(name);
	if(h) {
		h->Fill(valx, valy);
	} else {
		h = new TH2F(name, title, binsx, minx, maxx, binsy, miny, maxy);
		h->Fill(valx, valy);
		table->Add(h);
	}
}

void SabreEfficiency::CalculateEfficiency(const char* file) {
	std::cout<<"----------SABRE Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<file<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

	if(!dmap.IsValid()) {
		std::cerr<<"Unable to run SABRE Efficiency without a dead channel map."<<std::endl;
		std::cerr<<"If you have no dead channels, simply make a file that's empty"<<std::endl;
		std::cerr<<"Exiting."<<std::endl;
		std::cout<<"---------------------------------------------"<<std::endl;
	}

	switch(m_rxn_type) {
		case Mask::Kinematics::ONESTEP_DECAY:
		{
			RunDecay(file);
			break;
		}
		case Mask::Kinematics::TWOSTEP:
		{
			Run2Step(file);
			break;
		}
		case Mask::Kinematics::THREESTEP:
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

	Mask::NucData* eject = new Mask::NucData();
	Mask::NucData* resid = new Mask::NucData();


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

	Mask::Vec3 coords;
	double thetaIncident, eloss;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(eject->KE >= ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(eject->theta, eject->phi);
				if(chan.first != -1 && chan.second != -1) {
					if(dmap.IsDead(j, chan.first, 0) || dmap.IsDead(j, chan.second, 1)) break; //dead channel check
					coords = det.GetTrajectoryCoordinates(eject->theta, eject->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(eject->Z, eject->A, eject->KE, M_PI - thetaIncident);
					if((eject->KE - eloss) <= ENERGY_THRESHOLD) break; //deadlayer check

					eject_xs.push_back(coords.GetX());
					eject_ys.push_back(coords.GetY());
					eject_zs.push_back(coords.GetZ());
					break;
				}
			}
		}

		if(resid->KE > ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(resid->theta, resid->phi);
				if(chan.first != -1 && chan.second != -1) {
					if(dmap.IsDead(j, chan.first, 0) || dmap.IsDead(j, chan.second, 1)) break;
					coords = det.GetTrajectoryCoordinates(resid->theta, resid->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(resid->Z, resid->A, resid->KE, M_PI - thetaIncident);
					if((resid->KE - eloss) <= ENERGY_THRESHOLD) break;

					resid_xs.push_back(coords.GetX());
					resid_ys.push_back(coords.GetY());
					resid_zs.push_back(coords.GetZ());
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
	gde->SetMarkerStyle(1);
	gde->SetMarkerColor(2);

	TGraph2D* gr = new TGraph2D(ringxs.size(), &(ringxs[0]), &(ringys[0]), &(ringzs[0]));
	gr->SetName("ring_detector_edges");
	gr->SetTitle("SABRE Detector; x(m); y(m); z(m)");
	gr->SetMarkerStyle(1);

	TGraph2D* gw = new TGraph2D(wedgexs.size(), &(wedgexs[0]), &(wedgeys[0]), &(wedgezs[0]));
	gw->SetName("wedge_detector_edges");
	gw->SetTitle("SABRE Detector Wedges; x(m); y(m); z(m)");
	gw->SetMarkerStyle(1);

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

	Mask::NucData* break1 = new Mask::NucData();
	Mask::NucData* break2 = new Mask::NucData();

	THashTable* table = new THashTable();

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

	int cnt_09=0, cnt_08=0, cnt_07=0, cnt_06=0;
	double costheta_cm =0;

	Mask::Vec3 coords;
	double thetaIncident, eloss;
	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(break1->KE >= ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(break1->theta, break1->phi);
				if(chan.first != -1 && chan.second != -1) {
					if(dmap.IsDead(j, chan.first, 0) || dmap.IsDead(j, chan.second, 1)) break;
					costheta_cm = std::cos(break1->theta_cm);
					double this_bin;
					std::string this_name;
					for(int k=-9; k<=10; k++) {
						this_bin = 0.1*k;
						this_name = "counter_histo_bin"+std::to_string(k);
						if(costheta_cm > this_bin-0.05 && costheta_cm < this_bin+0.05) {
							MyFill(table, this_name.c_str(), "ctheta",20,-1.0,1.0,costheta_cm);
							break;
						} 
					}
					/*
					if(costheta_cm > -0.95 && costheta_cm < -0.85) {
						cnt_09++;
					} else if(costheta_cm > -0.85 && costheta_cm < -0.75) {
						cnt_08++;
					} else if(costheta_cm > -0.75 && costheta_cm < -0.65) {
						cnt_07++;
					} else if(costheta_cm > -0.65 && costheta_cm < -0.55) {
						cnt_06++;
					}
					*/
					coords = det.GetTrajectoryCoordinates(break1->theta, break1->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(break1->Z, break1->A, break1->KE, M_PI - thetaIncident);
					if((break1->KE - eloss) <= ENERGY_THRESHOLD) break;
					

					b1_thetas.push_back(break1->theta);
					b1_phis.push_back(break1->phi);
					b1_kes.push_back(break1->KE);
					break;
				}
			}
		}

		if(break2->KE > ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(break2->theta, break2->phi);
				if(chan.first != -1 && chan.second != -1) {
					if(dmap.IsDead(j, chan.first, 0) || dmap.IsDead(j, chan.second, 1)) break;
					coords = det.GetTrajectoryCoordinates(break2->theta, break2->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(break2->Z, break2->A, break2->KE, M_PI - thetaIncident);
					if((break2->KE - eloss) <= ENERGY_THRESHOLD) break;

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
	double b1eff_09 = cnt_09/nevents;
	double b1eff_08 = cnt_08/nevents;
	double b1eff_07 = cnt_07/nevents;
	double b1eff_06 = cnt_06/nevents;
	TParameter<double> break1_eff("Light Breakup Efficiency", b1eff);
	TParameter<double> break1_eff_09("Light Breakup Efficiency CosTheta -0.9", b1eff_09);
	TParameter<double> break1_eff_08("Light Breakup Efficiency CosTheta -0.8", b1eff_08);
	TParameter<double> break1_eff_07("Light Breakup Efficiency CosTheta -0.7", b1eff_07);
	TParameter<double> break1_eff_06("Light Breakup Efficiency CosTheta -0.6", b1eff_06);
	TParameter<double> break2_eff("Heavy Breakup Efficiency", b2eff);

	input->cd();
	table->Write();
	break1_eff.Write();
	break1_eff_09.Write();
	break1_eff_08.Write();
	break1_eff_07.Write();
	break1_eff_06.Write();
	break2_eff.Write();
	input->Close();

}

void SabreEfficiency::Run3Step(const char* file) {
	TFile* input = TFile::Open(file, "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	Mask::NucData* break1 = new Mask::NucData();
	Mask::NucData* break3 = new Mask::NucData();
	Mask::NucData* break4 = new Mask::NucData();


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

	bool break1_flag, break3_flag, break4_flag;
	int b1b3_count=0, b1b4_count=0, b3b4_count=0, b1b3b4_count=0;

	TH2F* b3b4_thetas = new TH2F("b3b4_theta_theta","b3b4_theta_theta;#theta_{3};#theta_{4}",180,0,180,180,0,180);
	TH2F* b3b4_kes = new TH2F("b3b4_ke_ke","b3b4_ke_ke;KE_{3};KE_{4}",150,0,15,150,0,15);
	TH2F* b3b4_phis = new TH2F("b3b4_phi_phi","b3b4_phi_phi;#phi_{3};#phi_{4}",360,0,360,360,0,360);
	Mask::Vec3 coords;
	double thetaIncident, eloss;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		break1_flag = false;
		break3_flag = false;
		break4_flag = false;

		tree->GetEntry(i);

		if(break1->KE > ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(break1->theta, break1->phi);
				if(chan.first != -1 && chan.second != -1) {
					coords = det.GetTrajectoryCoordinates(break1->theta, break1->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(break1->Z, break1->A, break1->KE, M_PI - thetaIncident);
					if((break1->KE - eloss) <= ENERGY_THRESHOLD) break;

					b1_thetas.push_back(break1->theta);
					b1_phis.push_back(break1->phi);
					b1_kes.push_back(break1->KE);
					break1_flag = true;
					break;
				}
			}
		}


		if(break3->KE > ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(break3->theta, break3->phi);
				if(chan.first != -1 && chan.second != -1) {
					coords = det.GetTrajectoryCoordinates(break3->theta, break3->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(break3->Z, break3->A, break3->KE, M_PI - thetaIncident);
					if((break3->KE - eloss) <= ENERGY_THRESHOLD) break;

					b3_thetas.push_back(break3->theta);
					b3_phis.push_back(break3->phi);
					b3_kes.push_back(break3->KE);
					break3_flag = true;
					break;
				}
			}
		}

		if(break4->KE > ENERGY_THRESHOLD) {
			for(int j=0; j<5; j++) {
				auto& det = detectors[j];
				auto chan = det.GetTrajectoryRingWedge(break4->theta, break4->phi);
				if(chan.first != -1 && chan.second != -1) {
					coords = det.GetTrajectoryCoordinates(break4->theta, break4->phi);
					thetaIncident = std::acos(coords.Dot(det.GetNormTilted())/(coords.GetR()));
					eloss = deadlayer.getEnergyLossTotal(break4->Z, break4->A, break4->KE, M_PI - thetaIncident);
					if((break4->KE - eloss) <= ENERGY_THRESHOLD) break;

					b4_thetas.push_back(break4->theta);
					b4_phis.push_back(break4->phi);
					b4_kes.push_back(break4->KE);
					break4_flag = true;
					break;
				}
			}
		}

		if(break1_flag && break3_flag && break4_flag) {
			b1b3b4_count++;
			b1b3_count++;
			b1b4_count++;
			b3b4_count++;
			b3b4_thetas->Fill(b3_thetas[b3_thetas.size()-1]/DEG2RAD, b4_thetas[b4_thetas.size()-1]/DEG2RAD);
			b3b4_kes->Fill(b3_kes[b3_kes.size()-1], b4_kes[b4_kes.size()-1]);
			b3b4_phis->Fill(b3_phis[b3_phis.size()-1]/DEG2RAD, b4_phis[b4_phis.size()-1]/DEG2RAD);
		} else if(break1_flag && break3_flag) {
			b1b3_count++;
		} else if(break1_flag && break4_flag) {
			b1b4_count++;
		} else if(break3_flag && break4_flag) {
			b3b4_count++;
			b3b4_thetas->Fill(b3_thetas[b3_thetas.size()-1]/DEG2RAD, b4_thetas[b4_thetas.size()-1]/DEG2RAD);
			b3b4_kes->Fill(b3_kes[b3_kes.size()-1], b4_kes[b4_kes.size()-1]);
			b3b4_phis->Fill(b3_phis[b3_phis.size()-1]/DEG2RAD, b4_phis[b4_phis.size()-1]/DEG2RAD);
		}
	}

	double b1eff = ((double) b1_thetas.size())/nevents;
	double b3eff = ((double) b3_thetas.size())/nevents;
	double b4eff = ((double) b4_thetas.size())/nevents;
	double b1b3eff = b1b3_count/nevents;
	double b1b4eff = b1b4_count/nevents;
	double b3b4eff = b3b4_count/nevents;
	double b1b3b4eff = b1b3b4_count/nevents;
	TParameter<double> break1_eff("Light Initial Breakup Efficiency", b1eff);
	TParameter<double> break3_eff("Light Final Breakup Efficiency", b3eff);
	TParameter<double> break4_eff("Heavy Final Breakup Efficiency", b4eff);
	TParameter<double> b1b3_eff("Break1 Coincident with Break3", b1b3eff);
	TParameter<double> b1b4_eff("Break1 Coincident with Break4", b1b4eff);
	TParameter<double> b3b4_eff("Break3 Coincident with Break4", b3b4eff);
	TParameter<double> b1b3b4_eff("All Breakups Coincident", b1b3b4eff);

	input->cd();
	break1_eff.Write();
	break3_eff.Write();
	break4_eff.Write();
	b1b3_eff.Write();
	b1b4_eff.Write();
	b3b4_eff.Write();
	b1b3b4_eff.Write();
	b3b4_thetas->Write();
	b3b4_phis->Write();
	b3b4_kes->Write();
	input->Close();
}
