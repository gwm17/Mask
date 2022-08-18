#include "MaskApp.h"
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

namespace Mask {

	MaskApp::MaskApp() :
		m_system(nullptr)
	{
		std::cout<<"----------GWM Kinematics Simulation----------"<<std::endl;
	}
	
	MaskApp::~MaskApp() 
	{
		delete m_system;
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
		input>>junk>>m_outputName;
	
		std::vector<int> avec, zvec, svec;
		int z, a, s;
		getline(input, junk);
		getline(input, junk);
		input>>junk>>junk;
		m_rxnType = GetRxnTypeFromString(junk);
		getline(input, junk);
		getline(input, junk);
		switch(m_rxnType) 
		{
			case RxnType::PureDecay:
			{
				m_system = new DecaySystem();
				for(int i=0; i<2; i++)
				{
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				break;
			}
			case RxnType::OneStepRxn:
			{
				m_system = new OneStepSystem();
				for(int i=0; i<3; i++)
				{
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				std::cout<<"here"<<std::endl;
				break;
			}
			case RxnType::TwoStepRxn:
			{
				m_system = new TwoStepSystem();
				for(int i=0; i<4; i++)
				{
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				break;
			}
			case RxnType::ThreeStepRxn:
			{
				m_system = new ThreeStepSystem();
				for(int i=0; i<5; i++)
				{
					input>>z>>a;
					avec.push_back(a);
					zvec.push_back(z);
				}
				break;
			}
			default:
				return false;
		}
		m_system->SetNuclei(zvec, avec);
	
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
			m_system->AddTargetLayer(zvec, avec, svec, thickness);
			input>>junk;
		}
		std::cout<<"Reaction equation: "<<GetSystemName()<<std::endl;
	
		double par1, par2;
		std::string dfile1, dfile2;
		getline(input, junk);
		getline(input, junk);
	
		input>>junk>>m_nsamples;
		input>>junk>>par1>>junk>>par2;
		m_system->SetBeamDistro(par1, par2);
		input>>junk>>par1;
		switch(m_rxnType) 
		{
			case RxnType::PureDecay : break;
			case RxnType::None : break;
			case RxnType::OneStepRxn :
			{
				dynamic_cast<OneStepSystem*>(m_system)->SetReactionThetaType(par1);
				break;
			}
			case RxnType::TwoStepRxn :
			{
				dynamic_cast<TwoStepSystem*>(m_system)->SetReactionThetaType(par1);
				break;
			}
			case RxnType::ThreeStepRxn :
			{
				dynamic_cast<ThreeStepSystem*>(m_system)->SetReactionThetaType(par1);
				break;
			}
		}
		input>>junk>>par1>>junk>>par2;
		m_system->SetTheta1Range(par1, par2);
		input>>junk>>par1>>junk>>par2;
		m_system->SetPhi1Range(par1, par2);
		input>>junk>>par1>>junk>>par2;
		m_system->SetExcitationDistro(par1, par2);
		input>>junk>>dfile1;
		input>>junk>>dfile2;
		switch(m_rxnType) 
		{
			case RxnType::PureDecay :
			{
				DecaySystem* this_sys = dynamic_cast<DecaySystem*>(m_system);
				this_sys->SetDecay1Distribution(dfile1);
				std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<std::endl;
				std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<std::endl;
				break;
			}
			case RxnType::None : break;
			case RxnType::OneStepRxn : break;
			case RxnType::TwoStepRxn :
			{
				TwoStepSystem* this_sys = dynamic_cast<TwoStepSystem*>(m_system);
				this_sys->SetDecay1Distribution(dfile1);
				std::cout<<"Decay1 angular momentum: "<<this_sys->GetDecay1AngularMomentum()<<std::endl;
				std::cout<<"Decay1 total branching ratio: "<<this_sys->GetDecay1BranchingRatio()<<std::endl;
				break;
			}
			case RxnType::ThreeStepRxn :
			{
				ThreeStepSystem* this_sys = dynamic_cast<ThreeStepSystem*>(m_system);
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
		if(m_system == nullptr) 
		{
			return;
		}
	
		TFile* output = TFile::Open(m_outputName.c_str(), "RECREATE");
		TTree* tree = new TTree("SimTree", "SimTree");
		tree->Branch("nuclei", m_system->GetNuclei());
	
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
	
			m_system->RunSystem();
			tree->Fill();
		}
	
		tree->Write(tree->GetName(), TObject::kOverwrite);
		output->Close();
		
		std::cout<<std::endl;
		std::cout<<"Complete."<<std::endl;
		std::cout<<"---------------------------------------------"<<std::endl;
	}

}
