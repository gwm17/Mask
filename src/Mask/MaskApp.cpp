#include "MaskApp.h"
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

namespace Mask {

	MaskApp::MaskApp() :
		m_system(nullptr), m_resources(nullptr)
	{
		std::cout<<"----------Monte Carlo Simulation of Kinematics----------"<<std::endl;
	}
	
	MaskApp::~MaskApp() 
	{
		delete m_system;
		for(std::size_t i=0; i<m_systemList.size(); i++)
			delete m_systemList[i];
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
		input>>junk>>m_nthreads;
		std::getline(input, junk);
		std::getline(input, junk);
		input>>junk>>m_nsamples;

		std::cout<<"Allocating resources... Asking for " << m_nthreads << " threads...";
		m_resources = std::make_unique<ThreadPool>(m_nthreads);
		std::cout<<" Complete."<<std::endl;

		std::cout<<"Outputing data to file: " << m_outputName <<std::endl;
		m_fileWriter.Open(m_outputName, "SimTree");

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
			std::cerr<<"Failure to parse reaction system... configuration not loaded."<<std::endl;
			return false;
		}

		for(uint32_t i=0; i<m_nthreads; i++)
		{
			m_systemList.push_back(CreateSystem(params));
			if(m_systemList.back() == nullptr || !m_systemList.back()->IsValid())
			{
				std::cerr<<"Failure to parse reaction system... configuration not loaded."<<std::endl;
				return false;
			}
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

		for(auto system : m_systemList)
		{
			system->SetLayeredTarget(target);
		}

		std::cout<<"Reaction equation: "<<m_system->GetSystemEquation()<<std::endl;
		std::cout<<"Number of samples: "<<m_nsamples<<std::endl;
	
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
		if(m_systemList.size() != m_nthreads)
		{
			return;
		}

		//Give our thread pool some tasks
		for(auto system : m_systemList)
			m_resources->PushJob({std::bind(&MaskApp::RunChunk, std::ref(*this), std::placeholders::_1), system});

		uint64_t count = 0;
		double percent = 0.05;
		uint64_t flushVal = m_nsamples*percent;
		uint64_t flushCount = 0;
		while(true)
		{
			if(count == flushVal)
			{
				count = 0;
				++flushCount;
				std::cout<<"\rPercent of data written to disk: "<<percent*flushCount*100<<"%"<<std::flush;
			}

			if(m_resources->IsFinished() && m_fileWriter.GetQueueSize() == 0)
				break;
			else if(m_fileWriter.Write())
				++count;
		}

		std::cout<<std::endl;
		std::cout<<"Complete."<<std::endl;
		std::cout<<"---------------------------------------------"<<std::endl;
	}
	
	void MaskApp::RunSingleThread()
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
		uint64_t percent5 = 0.05*m_nsamples;
		uint64_t count = 0;
		uint64_t npercent = 0;
	
		for(uint64_t i=0; i<m_nsamples; i++) 
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

	void MaskApp::RunChunk(ReactionSystem* system)
	{
		if(system == nullptr)
			return;

		uint64_t samples = m_nsamples / m_nthreads;

		for(uint64_t i=0; i<samples; i++)
		{
			system->RunSystem();
			m_fileWriter.PushData(*(system->GetNuclei()));
		}
	}

}
