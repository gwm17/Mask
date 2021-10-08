#include "MaskApp.h"
#include "MaskFile.h"
#include <fstream>
#include <iostream>

namespace Mask {

	MaskApp::MaskApp() :
		sys(nullptr)
	{
		std::cout<<"----------GWM Kinematics Simulation----------"<<std::endl;
	}
	
	MaskApp::~MaskApp() 
	{
		if(sys) delete sys;
	}
	
	bool MaskApp::LoadConfig(const std::string& filename) 
	{
		std::cout<<"Loading configuration in "<<filename<<"..."<<std::endl;
	
		std::ifstream input(filename);
		if(!input.is_open()) 
		{
			return false;
		}
	
		std::string junk;
		getline(input, junk);
		input>>junk>>m_outfile_name;
	
		std::vector<int> avec, zvec, svec;
		int z, a, s;
		getline(input, junk);
		getline(input, junk);
		input>>junk>>junk;
		m_rxn_type = GetRxnTypeFromString(junk);
		getline(input, junk);
		getline(input, junk);
		switch(m_rxn_type) 
		{
			case RxnType::PureDecay:
			{
				sys = new DecaySystem();
				m_rxn_type = RxnType::PureDecay;
				for(int i=0; i<2; i++) {
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				break;
			}
			case RxnType::OneStepRxn:
			{
				sys = new OneStepSystem();
				m_rxn_type = RxnType::OneStepRxn;
				for(int i=0; i<3; i++) {
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				break;
			}
			case RxnType::TwoStepRxn:
			{
				sys = new TwoStepSystem();
				m_rxn_type = RxnType::TwoStepRxn;
				for(int i=0; i<4; i++) {
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				break;
			}
			case RxnType::ThreeStepRxn:
			{
				sys = new ThreeStepSystem();
				m_rxn_type = RxnType::TwoStepRxn;
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
		input>>junk>>nlayers;
		for(int i=0; i<nlayers; i++) 
		{
			input>>junk>>junk>>thickness;
			avec.clear(); zvec.clear(); svec.clear();
			while(input>>junk) 
			{
				if(junk == "begin_elements")
				{
					input>>junk>>junk>>junk;
					continue;
				} 
				else if (junk == "end_elements")
					break;
				input>>z>>a>>s;
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
		switch(m_rxn_type) 
		{
			case RxnType::OneStepRxn :
			{
				dynamic_cast<OneStepSystem*>(sys)->SetReactionThetaType(par1);
				break;
			}
			case RxnType::TwoStepRxn :
			{
				dynamic_cast<TwoStepSystem*>(sys)->SetReactionThetaType(par1);
				break;
			}
			case RxnType::ThreeStepRxn :
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
		switch(m_rxn_type) 
		{
			case RxnType::OneStepRxn :
			{
				DecaySystem* this_sys = dynamic_cast<DecaySystem*>(sys);
				this_sys->SetDecay1Distribution(dfile1);
				std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<std::endl;
				std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<std::endl;
				break;
			}
			case RxnType::TwoStepRxn :
			{
				TwoStepSystem* this_sys = dynamic_cast<TwoStepSystem*>(sys);
				this_sys->SetDecay1Distribution(dfile1);
				std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<std::endl;
				std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<std::endl;
				break;
			}
			case RxnType::ThreeStepRxn :
			{
				ThreeStepSystem* this_sys = dynamic_cast<ThreeStepSystem*>(sys);
				this_sys->SetDecay1Distribution(dfile1);
				this_sys->SetDecay2Distribution(dfile2);
				std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<" Decay2 angular momentum: "<<this_sys->GetDecay2AngularMomentum()<<std::endl;
				std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<" Decay2 total branching ratio: "<<this_sys->GetDecay2BranchingRatio()<<std::endl;
				break;
			}
		}
	
		std::cout<<"Number of samples: "<<GetNumberOfSamples()<<std::endl;
	
		return true;
	}
	
	bool MaskApp::SaveConfig(const std::string& filename) { return true; }
	
	void MaskApp::Run() {
		std::cout<<"Running simulation..."<<std::endl;
		if(sys == nullptr) 
		{
			return;
		}
	
		MaskFile output(m_outfile_name, MaskFile::FileType::write);
		output.WriteHeader(m_rxn_type, m_nsamples);
	
		//For progress tracking
		uint32_t percent5 = 0.05*m_nsamples;
		uint32_t count = 0;
		uint32_t npercent = 0;
	
		for(uint32_t i=0; i<m_nsamples; i++) 
		{
			if(++count == percent5) 
			{
				npercent++;
				count = 0;
				std::cout<<"\rPercent complete: "<<npercent*5<<"%"<<std::flush;
			}
	
			sys->RunSystem();
			output.WriteData(sys->GetNuclei());
		}
	
		output.Close();
		
		std::cout<<std::endl;
		std::cout<<"Complete."<<std::endl;
		std::cout<<"---------------------------------------------"<<std::endl;
	}

}
