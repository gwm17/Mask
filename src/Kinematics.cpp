#include "Kinematics.h"
#include "MaskFile.h"
#include <fstream>
#include <iostream>

namespace Mask {

Kinematics::Kinematics() :
	sys(nullptr), save_tree_flag(false), do_plotter_flag(false), global_generator(new TRandom3(0))
{
	std::cout<<"----------GWM Kinematics Simulation----------"<<std::endl;
}

Kinematics::~Kinematics() {
	delete global_generator;
	if(sys) delete sys;
}

bool Kinematics::LoadConfig(const char* filename) {
	std::cout<<"Loading configuration in "<<filename<<"..."<<std::endl;

	std::ifstream input(filename);
	if(!input.is_open()) {
		std::cerr<<"Unable to load configuration in "<<filename<<", check that it exists"<<std::endl;
		return false;
	}

	std::string junk;
	getline(input, junk);
	input>>junk>>m_outfile_name;
	input>>junk>>junk;
	if(junk == "yes") save_tree_flag = true;
	input>>junk>>junk;
	if(junk == "yes") do_plotter_flag = true;

	std::vector<int> avec, zvec, svec;
	int z, a, s;
	getline(input, junk);
	getline(input, junk);
	input>>junk>>m_rxn_type;
	getline(input, junk);
	getline(input, junk);
	switch(m_rxn_type) {
		case 0:
		{
			sys = new DecaySystem();
			m_rxn_type = ONESTEP_DECAY;
			for(int i=0; i<2; i++) {
				input>>z>>a;
				avec.push_back(a);
				zvec.push_back(z);
			}
			break;
		}
		case 1:
		{
			sys = new OneStepSystem();
			m_rxn_type = ONESTEP_RXN;
			for(int i=0; i<3; i++) {
				input>>z>>a;
				avec.push_back(a);
				zvec.push_back(z);
			}
			break;
		}
		case 2:
		{
			sys = new TwoStepSystem();
			m_rxn_type = TWOSTEP;
			for(int i=0; i<4; i++) {
				input>>z>>a;
				avec.push_back(a);
				zvec.push_back(z);
			}
			break;
		}
		case 3:
		{
			sys = new ThreeStepSystem();
			m_rxn_type = THREESTEP;
			for(int i=0; i<5; i++) {
				input>>z>>a;
				avec.push_back(a);
				zvec.push_back(z);
			}
			break;
		}
		default:
			return false;
	}
	sys->SetNuclei(zvec, avec);

	int nlayers;
	double thickness;
	getline(input, junk);
	getline(input, junk);
	input>>junk>>junk;
	input>>junk>>nlayers;
	for(int i=0; i<nlayers; i++) {
		input>>junk>>junk>>thickness;
		getline(input, junk);
		getline(input, junk);
		avec.clear(); zvec.clear(); svec.clear();
		while(input>>z) {
			if(z == 0) break;
			input>>a>>s;
			zvec.push_back(z); avec.push_back(a); svec.push_back(s);
		}
		sys->AddTargetLayer(zvec, avec, svec, thickness);
		input>>junk;
	}
	std::cout<<"Reaction equation: "<<GetSystemName()<<std::endl;

	double par1, par2;
	std::string dfile1, dfile2;
	getline(input, junk);
	getline(input, junk);

	input>>junk>>m_nsamples;
	input>>junk>>par1>>junk>>par2;
	sys->SetBeamDistro(par1, par2);
	input>>junk>>par1;
	switch(m_rxn_type) {
		case ONESTEP_RXN :
		{
			dynamic_cast<OneStepSystem*>(sys)->SetReactionThetaType(par1);
			break;
		}
		case TWOSTEP:
		{
			dynamic_cast<TwoStepSystem*>(sys)->SetReactionThetaType(par1);
			break;
		}
		case THREESTEP:
		{
			dynamic_cast<ThreeStepSystem*>(sys)->SetReactionThetaType(par1);
			break;
		}
	}
	input>>junk>>par1>>junk>>par2;
	sys->SetTheta1Range(par1, par2);
	input>>junk>>par1>>junk>>par2;
	sys->SetPhi1Range(par1, par2);
	input>>junk>>par1>>junk>>par2;
	sys->SetExcitationDistro(par1, par2);
	input>>junk>>dfile1;
	input>>junk>>dfile2;
	switch(m_rxn_type) {
		case ONESTEP_DECAY:
		{
			DecaySystem* this_sys = dynamic_cast<DecaySystem*>(sys);
			this_sys->SetDecay1Distribution(dfile1);
			std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<std::endl;
			std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<std::endl;
			break;
		}
		case TWOSTEP:
		{
			TwoStepSystem* this_sys = dynamic_cast<TwoStepSystem*>(sys);
			this_sys->SetDecay1Distribution(dfile1);
			std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<std::endl;
			std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<std::endl;
			break;
		}
		case THREESTEP:
		{
			ThreeStepSystem* this_sys = dynamic_cast<ThreeStepSystem*>(sys);
			this_sys->SetDecay1Distribution(dfile1);
			this_sys->SetDecay2Distribution(dfile2);
			std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<" Decay2 angular momentum: "<<this_sys->GetDecay2AngularMomentum()<<std::endl;
			std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<" Decay2 total branching ratio: "<<this_sys->GetDecay2BranchingRatio()<<std::endl;
			break;
		}
	}
	sys->SetRandomGenerator(global_generator);

	std::cout<<"Number of samples: "<<GetNumberOfSamples()<<std::endl;

	return true;
}

bool Kinematics::SaveConfig(const char* filename) { return true; }

NucData Kinematics::ConvertNucleus(const Nucleus& nuc) {
	NucData datum;
	datum.E = nuc.GetE();
	datum.KE = nuc.GetKE();
	datum.p = nuc.GetP();
	datum.theta = nuc.GetTheta();
	datum.theta_cm = nuc.GetThetaCM();
	datum.phi = nuc.GetPhi();
	datum.Ex = nuc.GetExcitationEnergy();
	datum.Z = nuc.GetZ();
	datum.A = nuc.GetA();
	return datum;
}

void Kinematics::Run() {
	std::cout<<"Running simulation..."<<std::endl;
	switch(m_rxn_type) {
		case ONESTEP_DECAY:
		{
			RunOneStepDecay();
			break;
		}
		case ONESTEP_RXN:
		{
			RunOneStepRxn();
			break;
		}
		case TWOSTEP:
		{
			RunTwoStep();
			break;
		}
		case THREESTEP:
		{
			RunThreeStep();
			break;
		}
	}
	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}

void Kinematics::RunOneStepRxn() {
	OneStepSystem* this_sys = dynamic_cast<OneStepSystem*>(sys);
	if(this_sys == nullptr) {
		return;
	}

	
	TFile* output = TFile::Open(m_outfile_name.c_str(), "RECREATE");
	TTree* tree;
	NucData targ, proj, eject, residual;
	if(save_tree_flag) {
		tree = new TTree("DataTree","DataTree");
		tree->Branch("target", &targ);
		tree->Branch("projectile", &proj);
		tree->Branch("ejectile", &eject);
		tree->Branch("residual", &residual);
	}
	

	//For progress tracking
	int percent5 = 0.05*m_nsamples;
	int count = 0;
	int npercent = 0;

	std::vector<Nucleus> data;
	data.resize(4);
	for(int i=0; i<m_nsamples; i++) {
		if(++count == percent5) {//Show update every 5 percent
			npercent++;
			count = 0;
			std::cout<<"\rPercent complete: "<<npercent*5<<"%"<<std::flush;
		}

		this_sys->RunSystem();
		
		if(save_tree_flag) {
			targ = ConvertNucleus(this_sys->GetTarget());
			proj = ConvertNucleus(this_sys->GetProjectile());
			eject = ConvertNucleus(this_sys->GetEjectile());
			residual = ConvertNucleus(this_sys->GetResidual());
			tree->Fill();
		}
		if(do_plotter_flag) {
			plotman.FillData(this_sys->GetTarget());
			plotman.FillData(this_sys->GetProjectile());
			plotman.FillData(this_sys->GetEjectile());
			plotman.FillData(this_sys->GetResidual());
		}

	}

	
	output->cd();
	if(save_tree_flag) tree->Write(tree->GetName(), TObject::kOverwrite);
	if(do_plotter_flag) {
		plotman.GetTable()->Write();
		plotman.ClearTable();
	}
	output->Close();
	
}

void Kinematics::RunOneStepDecay() {
	DecaySystem* this_sys = dynamic_cast<DecaySystem*>(sys);
	if(this_sys == nullptr) {
		return;
	}

	TFile* output = TFile::Open(m_outfile_name.c_str(), "RECREATE");
	TTree* tree;
	NucData targ, eject, residual;
	if(save_tree_flag) {
		tree = new TTree("DataTree","DataTree");
		tree->Branch("target", &targ);
		tree->Branch("ejectile", &eject);
		tree->Branch("residual", &residual);
	}

	//For progress tracking
	int percent5 = 0.05*m_nsamples;
	int count = 0;
	int npercent = 0;

	for(int i=0; i<m_nsamples; i++) {
		if(++count == percent5) {//Show update every 5 percent
			npercent++;
			count = 0;
			std::cout<<"\rPercent complete: "<<npercent*5<<"%"<<std::flush;
		}

		this_sys->RunSystem();
		if(save_tree_flag) {
			targ = ConvertNucleus(this_sys->GetTarget());
			eject = ConvertNucleus(this_sys->GetEjectile());
			residual = ConvertNucleus(this_sys->GetResidual());
			tree->Fill();
		}
		if(do_plotter_flag) {
			plotman.FillData(this_sys->GetTarget());
			plotman.FillData(this_sys->GetEjectile());
			plotman.FillData(this_sys->GetResidual());
		}

	}

	output->cd();
	if(save_tree_flag) tree->Write(tree->GetName(), TObject::kOverwrite);
	if(do_plotter_flag) {
		plotman.GetTable()->Write();
		plotman.ClearTable();
	}
	output->Close();
}

void Kinematics::RunTwoStep() {

	TwoStepSystem* this_sys = dynamic_cast<TwoStepSystem*>(sys);
	if(this_sys == nullptr) {
		return;
	}

	
	TFile* output = TFile::Open(m_outfile_name.c_str(), "RECREATE");
	TTree* tree;
	NucData targ, proj, eject, residual, break1, break2;
	if(save_tree_flag) {
		tree = new TTree("DataTree","DataTree");
		tree->Branch("target", &targ);
		tree->Branch("projectile", &proj);
		tree->Branch("ejectile", &eject);
		tree->Branch("residual", &residual);
		tree->Branch("breakup1", &break1);
		tree->Branch("breakup2", &break2);
	}
	


	//For progress tracking
	int percent5 = 0.05*m_nsamples;
	int count = 0;
	int npercent = 0;

	for(int i=0; i<m_nsamples; i++) {
		if(++count == percent5) {//Show update every 5 percent
			npercent++;
			count = 0;
			std::cout<<"\rPercent complete: "<<npercent*5<<"%"<<std::flush;
		}

		this_sys->RunSystem();
		if(save_tree_flag) {
			targ = ConvertNucleus(this_sys->GetTarget());
			proj = ConvertNucleus(this_sys->GetProjectile());
			eject = ConvertNucleus(this_sys->GetEjectile());
			residual = ConvertNucleus(this_sys->GetResidual());
			break1 = ConvertNucleus(this_sys->GetBreakup1());
			break2 = ConvertNucleus(this_sys->GetBreakup2());
			tree->Fill();
		}
		if(do_plotter_flag) {
			plotman.FillData(this_sys->GetTarget());
			plotman.FillData(this_sys->GetProjectile());
			plotman.FillData(this_sys->GetEjectile());
			plotman.FillData(this_sys->GetResidual());
			plotman.FillData(this_sys->GetBreakup1());
			plotman.FillData(this_sys->GetBreakup2());
		}
		
	}

	
	output->cd();
	if(save_tree_flag) tree->Write(tree->GetName(), TObject::kOverwrite);
	if(do_plotter_flag) {
		plotman.GetTable()->Write();
		plotman.ClearTable();
	}
	output->Close();
	
}

void Kinematics::RunThreeStep() {

	ThreeStepSystem* this_sys = dynamic_cast<ThreeStepSystem*>(sys);
	if(this_sys == nullptr) {
		return;
	}

	TFile* output = TFile::Open(m_outfile_name.c_str(), "RECREATE");
	TTree* tree;
	NucData targ, proj, eject, residual, break1, break2, break3, break4;
	if(save_tree_flag) {
		tree = new TTree("DataTree","DataTree");
		tree->Branch("target", &targ);
		tree->Branch("projectile", &proj);
		tree->Branch("ejectile", &eject);
		tree->Branch("residual", &residual);
		tree->Branch("breakup1", &break1);
		tree->Branch("breakup2", &break2);
		tree->Branch("breakup3", &break3);
		tree->Branch("breakup4", &break4);
	}
	
	//For progress updating
	int percent5 = 0.05*m_nsamples;
	int count = 0;
	int npercent = 0;

	for(int i=0; i<m_nsamples; i++) {
		if(++count == percent5) {//Show update every 5 percent
			npercent++;
			count = 0;
			std::cout<<"\rPercent complete: "<<npercent*5<<"%"<<std::flush;
		}

		this_sys->RunSystem();
		if(save_tree_flag) {
			targ = ConvertNucleus(this_sys->GetTarget());
			proj = ConvertNucleus(this_sys->GetProjectile());
			eject = ConvertNucleus(this_sys->GetEjectile());
			residual = ConvertNucleus(this_sys->GetResidual());
			break1 = ConvertNucleus(this_sys->GetBreakup1());
			break2 = ConvertNucleus(this_sys->GetBreakup2());
			break3 = ConvertNucleus(this_sys->GetBreakup3());
			break4 = ConvertNucleus(this_sys->GetBreakup4());
			tree->Fill();
		}
		if(do_plotter_flag) {
			plotman.FillData(this_sys->GetTarget());
			plotman.FillData(this_sys->GetProjectile());
			plotman.FillData(this_sys->GetEjectile());
			plotman.FillData(this_sys->GetResidual());
			plotman.FillData(this_sys->GetBreakup1());
			plotman.FillData(this_sys->GetBreakup2());
			plotman.FillData(this_sys->GetBreakup3());
			plotman.FillData(this_sys->GetBreakup4());
		}
	}

	output->cd();
	if(save_tree_flag) tree->Write(tree->GetName(), TObject::kOverwrite);
	if(do_plotter_flag) {
		plotman.GetTable()->Write();
		plotman.ClearTable();
	}
	output->Close();
}

};
