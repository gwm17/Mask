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

		while(input>>junk)
		{
			if(junk == "begin_nuclei(Z,A)")
				continue;
			else if (junk == "end_nuclei(Z,A)")
				break;
			else
			{
				z = std::stoi(junk);
				input>>a;
				zvec.push_back(z);
				avec.push_back(a);
			}
		}

		m_system = CreateSystem(zvec, avec);
		if(m_system == nullptr)
		{
			std::cerr<<"Failure to parse reaction system... configuration not loaded."<<std::endl;
			return false;
		}

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
		std::string thetaTypeString;
		getline(input, junk);
		getline(input, junk);
	
		input>>junk>>m_nsamples;
		input>>junk>>par1>>junk>>par2;
		m_system->SetBeamDistro(par1, par2);
		input>>junk>>thetaTypeString;
		m_system->SetReactionThetaType(StringToRxnThetaType(thetaTypeString));
		
		input>>junk>>par1>>junk>>par2;
		m_system->SetTheta1Range(par1, par2);
		input>>junk>>par1>>junk>>par2;
		m_system->SetPhi1Range(par1, par2);
		input>>junk>>par1>>junk>>par2;
		m_system->SetExcitationDistro(par1, par2);
		input>>junk>>dfile1;
		input>>junk>>dfile2;
		m_system->SetDecay1Distribution(dfile1);
		m_system->SetDecay2Distribution(dfile2);
	
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
