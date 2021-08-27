#include "SabreEfficiency.h"
#include <fstream>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>

SabreEfficiency::SabreEfficiency() : 
	DetectorEfficiency(), deadlayer(DEADLAYER_THIN), sabre_eloss(SABRE_THICKNESS)
{
	detectors.reserve(5);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI0*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI1*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI2*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI3*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI4*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);

 	
	std::vector<int> dead_z = {14};
	std::vector<int> dead_a = {28};
	std::vector<int> dead_stoich = {1};
	deadlayer.SetElements(dead_z, dead_a, dead_stoich);
	sabre_eloss.SetElements(dead_z, dead_a, dead_stoich);
}

SabreEfficiency::~SabreEfficiency() {}

void SabreEfficiency::CalculateEfficiency(const std::string& file) {
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
		case Mask::Kinematics::ONESTEP_RXN:
		{
			Run1Step(file);
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

void SabreEfficiency::DrawDetectorSystem(const std::string& filename) {
	TFile* output = TFile::Open(filename.c_str(), "RECREATE");

	std::vector<double> ringxs, ringys, ringzs;
    std::vector<double> wedgexs, wedgeys, wedgezs;
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

 	TGraph2D* gr = new TGraph2D(ringxs.size(), &(ringxs[0]), &(ringys[0]), &(ringzs[0]));
	gr->SetName("ring_detector_edges");
	gr->SetTitle("SABRE Detector; x(m); y(m); z(m)");
	gr->SetMarkerStyle(1);

	TGraph2D* gw = new TGraph2D(wedgexs.size(), &(wedgexs[0]), &(wedgeys[0]), &(wedgezs[0]));
	gw->SetName("wedge_detector_edges");
	gw->SetTitle("SABRE Detector Wedges; x(m); y(m); z(m)");
	gw->SetMarkerStyle(1);

	TCanvas* canvas = new TCanvas();
	canvas->SetName("detector_system");
	canvas->cd();
	gr->Draw("AP");
	gw->Draw("same P");

	canvas->Write();
	gr->Write();
	gw->Write();

	output->Close();
}

double SabreEfficiency::RunConsistencyCheck() {
	double theta, phi;
	double npoints = 5.0*16.0*4.0;
	int count=0;
 	for(int h=0; h<5; h++) {
 		for(int j=0; j<16; j++) {
 			for(int k=0; k<4; k ++) {
 				theta  = detectors[h].GetRingTiltCoords(j, k).GetTheta();
 				phi = detectors[h].GetRingTiltCoords(j, k).GetPhi();
 				for(int i=0; i<5; i++) {
 					auto channels = detectors[i].GetTrajectoryRingWedge(theta, phi);
 					if(channels.first != -1) {
 						count++;
 					}
 				}
 			}
 		}
 	}

 	return ((double)count)/npoints;
}

/*Returns if detected, as well as total energy deposited in SABRE*/
std::pair<bool,double> SabreEfficiency::IsSabre(Mask::NucData* nucleus) {
	if(nucleus->KE <= ENERGY_THRESHOLD) {
		return std::make_pair(false, 0.0);
	}

	Mask::Vec3 coords;
	double thetaIncident, eloss, e_deposited;
	for(int i=0; i<5; i++) {
		auto chan = detectors[i].GetTrajectoryRingWedge(nucleus->theta, nucleus->phi);
		if(chan.first != -1 && chan.second != -1) {
			if(dmap.IsDead(i, chan.first, 0) || dmap.IsDead(i, chan.second, 1)) break; //dead channel check
			coords = detectors[i].GetTrajectoryCoordinates(nucleus->theta, nucleus->phi);
			thetaIncident = std::acos(coords.Dot(detectors[i].GetNormTilted())/(coords.GetR()));
			eloss = deadlayer.getEnergyLossTotal(nucleus->Z, nucleus->A, nucleus->KE, M_PI - thetaIncident);
			if((nucleus->KE - eloss) <= ENERGY_THRESHOLD) break; //deadlayer check
			e_deposited = sabre_eloss.getEnergyLossTotal(nucleus->Z, nucleus->A, nucleus->KE - eloss, M_PI - thetaIncident);
			return std::make_pair(true, e_deposited);
		}
	}

	return std::make_pair(false,0.0);
}

void SabreEfficiency::RunDecay(const std::string& filename) {
	TFile* input = TFile::Open(filename.c_str(), "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	THashTable* table = new THashTable();

	Mask::NucData* eject = new Mask::NucData();
	Mask::NucData* resid = new Mask::NucData();


	tree->SetBranchAddress("ejectile", &eject);
	tree->SetBranchAddress("residual", &resid);

	double nevents = tree->GetEntries();

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	std::pair<bool, double> eject_result, resid_result;
	int detected_eject=0, detected_resid=0;
	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		eject_result = IsSabre(eject);
		resid_result = IsSabre(resid);

		if(eject_result.first) {
			detected_eject++;
			MyFill(table, "detected_eject_theta_ke","detected ejectile theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,eject->theta/DEG2RAD,300,0.0,30.0,eject->KE);
			MyFill(table, "detected_eject_theta_edep","detected ejectile theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,eject->theta/DEG2RAD,300,0.0,30.0,eject_result.second);
		}

		if(resid_result.first) {
			detected_resid++;
			MyFill(table, "detected_resid_theta_ke","detected residual theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,resid->theta/DEG2RAD,300,0.0,30.0,resid->KE);
			MyFill(table, "detected_resid_theta_edep","detected residual theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,resid->theta/DEG2RAD,300,0.0,30.0,resid_result.second);
		}
	}

	double ejecteff = ((double) detected_eject)/nevents;
	double resideff = ((double) detected_resid)/nevents;
	TParameter<double> eject_eff("Ejectile Efficiency", ejecteff);
	TParameter<double> resid_eff("Residual Efficiency", resideff);

	input->cd();
	eject_eff.Write();
	resid_eff.Write();
	table->Write();
	input->Close();

}

void SabreEfficiency::Run1Step(const std::string& filename) {
	std::cout<<"SabreEfficiency::Run1Step Not Implemented!"<<std::endl;
	return;
}

void SabreEfficiency::Run2Step(const std::string& filename) {

	TFile* input = TFile::Open(filename.c_str(), "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	Mask::NucData* break1 = new Mask::NucData();
	Mask::NucData* break2 = new Mask::NucData();

	THashTable* table = new THashTable();

	tree->SetBranchAddress("breakup1", &break1);
	tree->SetBranchAddress("breakup2", &break2);

	double nevents = tree->GetEntries();

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	double costheta_cm =0;

	std::pair<bool, double> break1_result, break2_result;
	int detected_b1=0, detected_b2=0;

	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		break1_result = IsSabre(break1);
		break2_result = IsSabre(break2);

		if(break1_result.first) {
			detected_b1++;
			costheta_cm = std::cos(break1->theta_cm);
			MyFill(table,"detected_break1_cm_theta","cos(#theta_{CM})",20,-1.0,1.0,costheta_cm);
			MyFill(table, "detected_break1_theta_ke","detected break1 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break1->theta/DEG2RAD,300,0.0,30.0,break1->KE);
			MyFill(table, "detected_break1_theta_edep","detected break1 theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,break1->theta/DEG2RAD,300,0.0,30.0,break1_result.second);
		}

		if(break2_result.first) {
			detected_b2++;
			MyFill(table, "detected_break2_theta_ke","detected break2 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break2->theta/DEG2RAD,300,0.0,30.0,break2->KE);
			MyFill(table, "detected_break2_theta_edep","detected break2 theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,break2->theta/DEG2RAD,300,0.0,30.0,break2_result.second);
		
		}
	}

	double b1eff = ((double) detected_b1)/nevents;
	double b2eff = ((double) detected_b2)/nevents;
	TParameter<double> break1_eff("Breakup1 Efficiency", b1eff);
	TParameter<double> break2_eff("Breakup2 Efficiency", b2eff);

	input->cd();
	table->Write();
	break1_eff.Write();
	break2_eff.Write();
	input->Close();

}

void SabreEfficiency::Run3Step(const std::string& filename) {
	TFile* input = TFile::Open(filename.c_str(), "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	THashTable* table = new THashTable();

	Mask::NucData* break1 = new Mask::NucData();
	Mask::NucData* break3 = new Mask::NucData();
	Mask::NucData* break4 = new Mask::NucData();


	tree->SetBranchAddress("breakup1", &break1);
	tree->SetBranchAddress("breakup3", &break3);
	tree->SetBranchAddress("breakup4", &break4);

	double nevents = tree->GetEntries();

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	std::pair<bool, double> break1_result, break3_result, break4_result;
	int detected_b1=0, detected_b3=0, detected_b4=0;
	int b1b3_count=0, b1b4_count=0, b3b4_count=0, b1b3b4_count=0;


	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		break1_result = IsSabre(break1);
		break3_result = IsSabre(break3);
		break4_result = IsSabre(break4);


		tree->GetEntry(i);

		if(break1_result.first) {
			detected_b1++;
			MyFill(table, "detected_break1_theta_ke","detected break1 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break1->theta/DEG2RAD,300,0.0,30.0,break1->KE);
			MyFill(table, "detected_break1_theta_edep","detected break1 theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,break1->theta/DEG2RAD,300,0.0,30.0,break1_result.second);
		}

		if(break3_result.first) {
			detected_b3++;
			MyFill(table, "detected_break3_theta_ke","detected break3 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break3->theta/DEG2RAD,300,0.0,30.0,break3->KE);
			MyFill(table, "detected_break3_theta_edep","detected break3 theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,break3->theta/DEG2RAD,300,0.0,30.0,break3_result.second);
		}

		if(break4_result.first) {
			detected_b1++;
			MyFill(table, "detected_break4_theta_ke","detected break4 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break4->theta/DEG2RAD,300,0.0,30.0,break4->KE);
			MyFill(table, "detected_break4_theta_edep","detected break4 theta vs edep;#theta (deg);E deposited(MeV)",180,0.0,180.0,break4->theta/DEG2RAD,300,0.0,30.0,break4_result.second);
		}

		if(break1_result.first && break3_result.first && break4_result.first) {
			b1b3b4_count++;
			b1b3_count++;
			b1b4_count++;
			b3b4_count++;
			MyFill(table,"b3b4_theta_theta","b3b4_theta_theta;#theta_{3};#theta_{4}",180,0.0,180.0,break3->theta/DEG2RAD,180,0,180,break4->theta/DEG2RAD);
			MyFill(table,"b3b4_ke_ke","b3b4_ke_ke;KE_{3};KE_{4}",150,0.0,15.0,break3->KE,150,0.0,15.0,break4->KE);
			MyFill(table,"b3b4_phi_phi","b3b4_phi_phi;#phi_{3};#phi_{4}",360,0.0,360.0,break3->phi/DEG2RAD,360,0.0,360.0,break4->phi/DEG2RAD);
		} else if(break1_result.first && break3_result.first) {
			b1b3_count++;
		} else if(break1_result.first && break4_result.first) {
			b1b4_count++;
		} else if(break3_result.first && break4_result.first) {
			b3b4_count++;
			MyFill(table,"b3b4_theta_theta","b3b4_theta_theta;#theta_{3};#theta_{4}",180,0.0,180.0,break3->theta/DEG2RAD,180,0,180,break4->theta/DEG2RAD);
			MyFill(table,"b3b4_ke_ke","b3b4_ke_ke;KE_{3};KE_{4}",150,0.0,15.0,break3->KE,150,0.0,15.0,break4->KE);
			MyFill(table,"b3b4_phi_phi","b3b4_phi_phi;#phi_{3};#phi_{4}",360,0.0,360.0,break3->phi/DEG2RAD,360,0.0,360.0,break4->phi/DEG2RAD);
		}
	}

	double b1eff = ((double) detected_b1)/nevents;
	double b3eff = ((double) detected_b3)/nevents;
	double b4eff = ((double) detected_b4)/nevents;
	double b1b3eff = b1b3_count/nevents;
	double b1b4eff = b1b4_count/nevents;
	double b3b4eff = b3b4_count/nevents;
	double b1b3b4eff = b1b3b4_count/nevents;
	TParameter<double> break1_eff("Breakup1 Efficiency", b1eff);
	TParameter<double> break3_eff("Breakup3 Efficiency", b3eff);
	TParameter<double> break4_eff("Breakup4 Efficiency", b4eff);
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
	table->Write();

	input->Close();
}
