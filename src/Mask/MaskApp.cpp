#include "MaskApp.h"
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

namespace Mask {

	MaskApp::MaskApp() :
		m_system(nullptr)
	{
		std::cout<<"----------Monte Carlo Simulation of Kinematics----------"<<std::endl;
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
		std::getline(input, junk);
		input>>junk>>m_outputName;
		std::getline(input, junk);
		std::getline(input, junk);
		input>>junk>>m_nsamples;

		std::vector<StepParameters> params;
		int z, a;
		while(input>>junk)
		{
			if(junk == "begin_chain")
				continue;
			else if (junk == "end_chain")
				break;
			else if(junk == "begin_step")
			{
				StepParameters currentParams;
				input >> junk >> junk;
				currentParams.rxnType = StringToRxnType(junk);
				if(currentParams.rxnType == RxnType::Reaction)
				{
					input >> junk;
					for(int i=0; i<3; i++)
					{
						input >> z >> a;
						currentParams.Z.push_back(z);
						currentParams.A.push_back(a);
					}
					input >> junk;
					input >> junk >> currentParams.meanBeamEnergy;
					input >> junk >> currentParams.sigmaBeamEnergy;
					input >> junk >> junk;
					currentParams.thetaType = StringToRxnThetaType(junk);
					input >> junk >> currentParams.thetaMin;
					input >> junk >> currentParams.thetaMax;
					input >> junk >> currentParams.phiMin;
					input >> junk >> currentParams.phiMax;
					input >> junk >> currentParams.meanResidualEx;
					input >> junk >> currentParams.sigmaResidualEx;
					params.push_back(currentParams);
				}
				else if(currentParams.rxnType == RxnType::Decay)
				{
					input >> junk;
					for(int i=0; i<2; i++)
					{
						input >> z >> a;
						currentParams.Z.push_back(z);
						currentParams.A.push_back(a);
					}
					input >> junk;
					input >> junk >> currentParams.phiMin;
					input >> junk >> currentParams.phiMax;
					input >> junk >> currentParams.meanResidualEx;
					input >> junk >> currentParams.sigmaResidualEx;
					input >> junk >> currentParams.angularDistFile;
					params.push_back(currentParams);
				}
				else
				{
					std::cerr << "Invalid reaction information at MaskApp::LoadConfig!" << std::endl;
					return false;
				}
			}
		}

		m_system = CreateSystem(params);
		if(m_system == nullptr || !m_system->IsValid())
		{
			std::cerr<<"Param size: "<<params.size()<<std::endl;
			std::cerr<<"Failure to parse reaction system... configuration not loaded."<<std::endl;
			return false;
		}
		std::getline(input, junk);
		std::getline(input, junk);

		LayeredTarget target;
		int nlayers;
		double thickness;
		std::vector<int> avec, zvec, svec;
		int s;
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
			target.AddLayer(zvec, avec, svec, thickness);
			input>>junk;
		}
		m_system->SetLayeredTarget(target);

		std::cout<<"Outputing data to file: "<<m_outputName<<std::endl;
		std::cout<<"Reaction equation: "<<GetSystemName()<<std::endl;
		std::cout<<"Number of samples: "<<GetNumberOfSamples()<<std::endl;
	
		return true;
	}
	
	//Not implemented... yet!
	bool MaskApp::SaveConfig(const std::string& filename)
	{ 
		return true;
	}
	
	void MaskApp::Run()
	{
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
