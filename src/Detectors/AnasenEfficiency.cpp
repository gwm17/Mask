#include "AnasenEfficiency.h"
#include "Kinematics.h"
#include <TH2.h>
#include <TH1.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TParameter.h>

AnasenEfficiency::AnasenEfficiency() {
	for(int i=0; i<n_sx3_per_ring; i++) {
		m_Ring1.emplace_back(4, sx3_length, sx3_width, ring_phi[i], ring1_z, ring_rho[i]);
		m_Ring2.emplace_back(4, sx3_length, sx3_width, ring_phi[i], -1.0*ring1_z, ring_rho[i]);
	}
	for(int i=0; i<n_qqq; i++) {
		m_forwardQQQs.emplace_back(qqq_rinner, qqq_router, qqq_deltaphi, qqq_phi[i], qqq_z[i]);
		m_backwardQQQs.emplace_back(qqq_rinner, qqq_router, qqq_deltaphi, qqq_phi[i], (-1.0)*qqq_z[i]);
	}
}

AnasenEfficiency::~AnasenEfficiency() {}

void AnasenEfficiency::MyFill(THashTable* table, const std::string& name, const std::string& title, int bins, float min, float max, double val) {
	TH1F* h = (TH1F*) table->FindObject(name.c_str());
	if(h) {
		h->Fill(val);
	} else {
		h = new TH1F(name.c_str(), title.c_str(), bins, min, max);
		h->Fill(val);
		table->Add(h);
	}
}

void AnasenEfficiency::MyFill(THashTable* table, const std::string& name, const std::string& title, int binsx, float minx, float maxx, double valx, int binsy, float miny, float maxy, double valy) {
	TH2F* h = (TH2F*) table->FindObject(name.c_str());
	if(h) {
		h->Fill(valx, valy);
	} else {
		h = new TH2F(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
		h->Fill(valx, valy);
		table->Add(h);
	}
}

void AnasenEfficiency::DrawDetectorSystem(const std::string& filename) {
	TFile* file = TFile::Open(filename.c_str(), "RECREATE");

	std::vector<double> x, y, z;
	std::vector<double> cx, cy, cz;
	Mask::Vec3 coords;
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			for(int k=0; k<4; k++) {
				coords = m_Ring1[i].GetRotatedFrontStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_Ring1[i].GetRotatedBackStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			coords = m_Ring1[i].GetHitCoordinates(j, 0);
			cx.push_back(coords.GetX());
			cy.push_back(coords.GetY());
			cz.push_back(coords.GetZ());
		}
	}
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			for(int k=0; k<4; k++) {
				coords = m_Ring2[i].GetRotatedFrontStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_Ring2[i].GetRotatedBackStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			coords = m_Ring2[i].GetHitCoordinates(j, 0);
			cx.push_back(coords.GetX());
			cy.push_back(coords.GetY());
			cz.push_back(coords.GetZ());
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<4; k++) {
				coords = m_forwardQQQs[i].GetRingCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_forwardQQQs[i].GetWedgeCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			for(int k=0; k<16; k++) {
				coords = m_forwardQQQs[i].GetHitCoordinates(j, k);
				cx.push_back(coords.GetX());
				cy.push_back(coords.GetY());
				cz.push_back(coords.GetZ());
			}
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<4; k++) {
				coords = m_backwardQQQs[i].GetRingCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_backwardQQQs[i].GetWedgeCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			for(int k=0; k<16; k++) {
				coords = m_backwardQQQs[i].GetHitCoordinates(j, k);
				cx.push_back(coords.GetX());
				cy.push_back(coords.GetY());
				cz.push_back(coords.GetZ());
			}
		}
	}
	TGraph2D* graph = new TGraph2D(x.size(), &(x[0]), &(y[0]), &(z[0]));
	graph->SetName("CornerGraph");
	graph->SetMarkerStyle(2);
	graph->SetLineColor(1);

	TGraph2D* graph2 = new TGraph2D(cx.size(), &(cx[0]), &(cy[0]), &(cz[0]));
	graph2->SetName("CenterGraph");
	graph2->SetMarkerStyle(2);
	graph2->SetMarkerColor(4);

	TCanvas* canvas = new TCanvas();
	canvas->SetName("ANASEN Detector");
	graph->Draw("A|P");
	graph2->Draw("same|P");
	canvas->Write();
	graph->Write();
	graph2->Write();
	file->Close();
}

double AnasenEfficiency::RunConsistencyCheck() {
	std::vector<Mask::Vec3> r1_points;
	std::vector<Mask::Vec3> r2_points;
	std::vector<Mask::Vec3> fqqq_points;
	std::vector<Mask::Vec3> bqqq_points;
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			r1_points.push_back(m_Ring1[i].GetHitCoordinates(j, 0));
		}
	}
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			r2_points.push_back(m_Ring2[i].GetHitCoordinates(j, 0));
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<16; k++) {
				fqqq_points.push_back(m_forwardQQQs[i].GetHitCoordinates(j, k));
			}
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<16; k++) {
				bqqq_points.push_back(m_backwardQQQs[i].GetHitCoordinates(j, k));
			}
		}
	}

	int npoints = r1_points.size() + r2_points.size() + fqqq_points.size() + bqqq_points.size();
	int count = 0;
	Mask::Vec3 coords;
	for(auto& point : r1_points) {
		for(auto& sx3 : m_Ring1) {
			auto result = sx3.GetChannelRatio(point.GetTheta(), point.GetPhi());
			coords = sx3.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}
	for(auto& point : r2_points) {
		for(auto& sx3 : m_Ring2) {
			auto result = sx3.GetChannelRatio(point.GetTheta(), point.GetPhi());
			coords = sx3.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}
	for(auto& point : fqqq_points) {
		for(auto& qqq : m_forwardQQQs) {
			auto result = qqq.GetTrajectoryRingWedge(point.GetTheta(), point.GetPhi());
			coords = qqq.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}
	for(auto& point : bqqq_points) {
		for(auto& qqq : m_backwardQQQs) {
			auto result = qqq.GetTrajectoryRingWedge(point.GetTheta(), point.GetPhi());
			coords = qqq.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}

	double ratio = ((double)count)/((double)npoints);

	return ratio;

}

bool AnasenEfficiency::IsRing1(double theta, double phi) {
	for(auto& sx3 : m_Ring1) {
		auto result = sx3.GetChannelRatio(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}

	return false;
}

bool AnasenEfficiency::IsRing2(double theta, double phi) {
	for(auto& sx3 : m_Ring2) {
		auto result = sx3.GetChannelRatio(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}

	return false;
}

bool AnasenEfficiency::IsQQQ(double theta, double phi) {
	for(auto& qqq : m_forwardQQQs) {
		auto result = qqq.GetTrajectoryRingWedge(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}

	
	for(auto& qqq : m_backwardQQQs) {
		auto result = qqq.GetTrajectoryRingWedge(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}



	return false;
}

void AnasenEfficiency::CalculateEfficiency(const std::string& file) {
	std::cout<<"----------ANASEN Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<file<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

	switch(m_rxn_type) {
		case Mask::Kinematics::ONESTEP_DECAY:
		{
			RunDecay(file);
			break;
		}
		case Mask::Kinematics::TWOSTEP:
		{
			RunTwoStep(file);
			break;
		}
		case Mask::Kinematics::THREESTEP:
		{
			RunThreeStep(file);
			break;
		}
	}
	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}

void AnasenEfficiency::RunDecay(const std::string& filename) {
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

	int detected_eject=0;
	int detected_resid=0;
	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		tree->GetEntry(i);

		if(eject->KE >= threshold && (IsRing1(eject->theta, eject->phi) || IsRing2(eject->theta, eject->phi) || IsQQQ(eject->theta, eject->phi))) {
			detected_eject++;
			MyFill(table, "detected_eject_theta_ke","detected ejectile theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,eject->theta/deg2rad,300,0.0,30.0,eject->KE);
		}

		if(resid->KE > threshold && (IsRing1(resid->theta, resid->phi) || IsRing2(resid->theta, resid->phi) || IsQQQ(resid->theta, resid->phi))) {
			detected_resid++;
			MyFill(table, "detected_resid_theta_ke","detected residual theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,resid->theta/deg2rad,300,0.0,30.0,resid->KE);
		}

	}

	double ejecteff = ((double) detected_eject)/nevents;
	double resideff = ((double) detected_resid)/nevents;
	TParameter<double> eject_eff("Light Breakup Efficiency", ejecteff);
	TParameter<double> resid_eff("Heavy Breakup Efficiency", resideff);

	input->cd();
	eject_eff.Write();
	resid_eff.Write();
	table->Write();
	input->Close();
}

void AnasenEfficiency::RunTwoStep(const std::string& filename) {

	TFile* input = TFile::Open(filename.c_str(), "UPDATE");
	TTree* tree = (TTree*) input->Get("DataTree");

	Mask::NucData* eject = new Mask::NucData();
	Mask::NucData* break1 = new Mask::NucData();
	Mask::NucData* break2 = new Mask::NucData();

	THashTable* table = new THashTable();

	tree->SetBranchAddress("ejectile", &eject);
	tree->SetBranchAddress("breakup1", &break1);
	tree->SetBranchAddress("breakup2", &break2);

	double nevents = tree->GetEntries();

	//Progress tracking
	int percent5 = nevents*0.05;
	int count = 0;
	int npercent = 0;

	double costheta_cm =0;

	bool break1_flag, break2_flag, eject_flag;
	int b1_count=0, b2_count=0, eject_count=0;
	int b1b2_count=0, b1e_count=0, b2e_count=0, b1b2e_count=0;
	for(int i=0; i<tree->GetEntries(); i++) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		break1_flag = false;
		break2_flag = false;
		eject_flag = false;

		tree->GetEntry(i);

		if(eject->KE >= threshold && (IsRing1(eject->theta, eject->phi) || IsRing2(eject->theta, eject->phi) || IsQQQ(eject->theta, eject->phi))) {
			eject_count++;
			eject_flag = true;
			MyFill(table, "detected_eject_theta_ke","detected ejectile theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,eject->theta/deg2rad,300,0.0,30.0,eject->KE);
		}

		if(break1->KE >= threshold && (IsRing1(break1->theta, break1->phi) || IsRing2(break1->theta, break1->phi) || IsQQQ(break1->theta, break1->phi))) {
			b1_count++;
			break1_flag = true;
			costheta_cm = std::cos(break1->theta_cm);
			MyFill(table,"detected_break1_cm_theta","detected breakup1 CM theta;cos(#theta_{CM})",20,-1.0,1.0,costheta_cm);
			MyFill(table, "detected_break1_theta_ke","detected breakup1 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break1->theta/deg2rad,300,0.0,30.0,break1->KE);
		}

		if(break2->KE >= threshold && (IsRing1(break2->theta, break2->phi) || IsRing2(break2->theta, break2->phi) || IsQQQ(break2->theta, break2->phi))) {
			b2_count++;
			break2_flag = true;
			MyFill(table, "detected_break2_theta_ke","detected breakup2 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break2->theta/deg2rad,300,0.0,30.0,break2->KE);
		}

		if(break1_flag && break2_flag && eject_flag) {
			b1b2e_count++;
			b1b2_count++;
			b1e_count++;
			b2e_count++;
			MyFill(table,"b1b2_theta_theta","b1b2_theta_theta;#theta_{3};#theta_{4}",180,0.0,180.0,break1->theta/deg2rad,180,0.0,180.0,break2->theta/deg2rad);
			MyFill(table,"b1b2_ke_ke","b1b2_ke_ke;KE_{3};KE_{4}",300,0.0,30.0,break1->KE,300,0.0,30.0,break2->KE);
			MyFill(table,"b1b2_phi_phi","b1b2_phi_phi;#phi_{3};#phi_{4}",360,0.0,360.0,break1->phi/deg2rad,360,0.0,360.0,break2->phi/deg2rad);
		} else if(break1_flag && eject_flag) {
			b1e_count++;
		} else if(break2_flag && eject_flag) {
			b2e_count++;
		} else if(break1_flag && break2_flag) {
			b1b2_count++;
			MyFill(table,"b1b2_theta_theta","b1b2_theta_theta;#theta_{3};#theta_{4}",180,0.0,180.0,break1->theta/deg2rad,180,0.0,180.0,break2->theta/deg2rad);
			MyFill(table,"b1b2_ke_ke","b1b2_ke_ke;KE_{3};KE_{4}",300,0.0,30.0,break1->KE,300,0.0,30.0,break2->KE);
			MyFill(table,"b1b2_phi_phi","b1b2_phi_phi;#phi_{3};#phi_{4}",360,0.0,360.0,break1->phi/deg2rad,360,0.0,360.0,break2->phi/deg2rad);
		}
	}

	double eeff = ((double) eject_count)/nevents;
	double b1eff = ((double) b1_count)/nevents;
	double b2eff = ((double) b2_count)/nevents;
	double b1b2eff = b1b2_count/nevents;
	double b1eeff = b1e_count/nevents;
	double b2eeff = b2e_count/nevents;
	double b1b2eeff = b1b2e_count/nevents;
	TParameter<double> eject_eff("Ejectile Efficiency", eeff);
	TParameter<double> break1_eff("Breakup 1 Efficiency", b1eff);
	TParameter<double> break2_eff("Breakup 2 Efficiency", b2eff);
	TParameter<double> break1break2_eff("Breakup1 Coincident with Breakup2", b1b2eff);
	TParameter<double> break1eject_eff("Breakup1 Coincident with Ejectile", b1eeff);
	TParameter<double> break2eject_eff("Breakup2 Coincident with Ejectile", b2eeff);
	TParameter<double> break1break2eject_eff("All Three Particles Coincident", b1b2eeff);


	input->cd();
	table->Write();
	eject_eff.Write();
	break1_eff.Write();
	break2_eff.Write();
	break1break2_eff.Write();
	break1eject_eff.Write();
	break2eject_eff.Write();
	break1break2eject_eff.Write();

	input->Close();
}

void AnasenEfficiency::RunThreeStep(const std::string& filename) {
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

	bool break1_flag, break3_flag, break4_flag;
	int b1_count=0, b3_count=0, b4_count=0;
	int b1b3_count=0, b1b4_count=0, b3b4_count=0, b1b3b4_count=0;

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

		if(break1->KE >= threshold && (IsRing1(break1->theta, break1->phi) || IsRing2(break1->theta, break1->phi) || IsQQQ(break1->theta, break1->phi))) {
			b1_count++;
			break1_flag = true;
			MyFill(table, "detected_break1_theta_ke","detected breakup1 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break1->theta/deg2rad,300,0.0,30.0,break1->KE);
		}


		if(break3->KE >= threshold && (IsRing1(break3->theta, break3->phi) || IsRing2(break3->theta, break3->phi) || IsQQQ(break3->theta, break3->phi))) {
			b3_count++;
			break3_flag = true;
			MyFill(table, "detected_break3_theta_ke","detected breakup3 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break3->theta/deg2rad,300,0.0,30.0,break3->KE);
		}

		if(break4->KE >= threshold && (IsRing1(break4->theta, break4->phi) || IsRing2(break4->theta, break4->phi) || IsQQQ(break4->theta, break4->phi))) {
			b4_count++;
			break4_flag = true;
			MyFill(table, "detected_break4_theta_ke","detected breakup4 theta vs ke;#theta (deg);KE(MeV)",180,0.0,180.0,break4->theta/deg2rad,300,0.0,30.0,break4->KE);
		}

		if(break1_flag && break3_flag && break4_flag) {
			b1b3b4_count++;
			b1b3_count++;
			b1b4_count++;
			b3b4_count++;
			MyFill(table,"b3b4_theta_theta","b3b4_theta_theta;#theta_{3};#theta_{4}",180,0.0,180.0,break3->theta/deg2rad,180,0.0,180.0,break4->theta/deg2rad);
			MyFill(table,"b3b4_ke_ke","b3b4_ke_ke;KE_{3};KE_{4}",300,0.0,30.0,break3->KE,300,0.0,30.0,break4->KE);
			MyFill(table,"b3b4_phi_phi","b3b4_phi_phi;#phi_{3};#phi_{4}",360,0.0,360.0,break3->phi/deg2rad,360,0.0,360.0,break4->phi/deg2rad);
		} else if(break1_flag && break3_flag) {
			b1b3_count++;
		} else if(break1_flag && break4_flag) {
			b1b4_count++;
		} else if(break3_flag && break4_flag) {
			b3b4_count++;
			MyFill(table,"b3b4_theta_theta","b3b4_theta_theta;#theta_{3};#theta_{4}",180,0.0,180.0,break3->theta/deg2rad,180,0.0,180.0,break4->theta/deg2rad);
			MyFill(table,"b3b4_ke_ke","b3b4_ke_ke;KE_{3};KE_{4}",300,0.0,30.0,break3->KE,300,0.0,30.0,break4->KE);
			MyFill(table,"b3b4_phi_phi","b3b4_phi_phi;#phi_{3};#phi_{4}",360,0.0,360.0,break3->phi/deg2rad,360,0.0,360.0,break4->phi/deg2rad);
		}
	}

	double b1eff = ((double) b1_count)/nevents;
	double b3eff = ((double) b3_count)/nevents;
	double b4eff = ((double) b4_count)/nevents;
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