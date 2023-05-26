#include "MaskApp.h"
#include "ConfigSerializer.h"
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "yaml-cpp/yaml.h"

namespace Mask {

	MaskApp::MaskApp() :
		m_resources(nullptr)
	{
		std::cout<<"----------Monte Carlo Simulation of Kinematics----------"<<std::endl;
	}
	
	MaskApp::~MaskApp() 
	{
		for(std::size_t i=0; i<m_systemList.size(); i++)
			delete m_systemList[i];
	}

	bool MaskApp::LoadConfig(const std::string& filename)
	{
		std::cout << "Loading configuration in " << filename << std::endl;
		
		if (!ConfigSerializer::DeserializeConfig(filename, m_params))
		{
			return false;
		}

		//Initialize the system
		//Reaction chain
		for(uint32_t i=0; i<m_params.nThreads; i++)
		{
			m_systemList.push_back(CreateSystem(m_params.chainParams));
			if(m_systemList.back() == nullptr || !m_systemList.back()->IsValid())
			{
				std::cerr<<"Failure to parse reaction system... configuration not loaded."<<std::endl;
				return false;
			}
		}
		//Link the target
		for(auto system : m_systemList)
		{
			system->SetLayeredTarget(m_params.target);
		}
		//Setup threading
		m_resources = std::make_unique<ThreadPool<ReactionSystem*, uint64_t>>(m_params.nThreads);
		//Little bit of integer division mangling to make sure we do the total number of samples
    	uint64_t quotient = m_params.nSamples / m_params.nThreads;
    	uint64_t remainder = m_params.nSamples % m_params.nThreads;
    	m_chunkSamples.push_back(quotient + remainder);
    	for(uint64_t i=1; i<m_params.nThreads; i++)
        	m_chunkSamples.push_back(quotient);
		m_fileWriter.Open(m_params.outputFileName, "SimTree");


		std::cout << "Reaction equation: " << m_systemList[0]->GetSystemEquation() << std::endl;
		std::cout << "Number of samples: " << m_params.nSamples << std::endl;
		std::cout << "Number of threads: " << m_params.nThreads << std::endl;
		std::cout << "Outputing data to file: " << m_params.outputFileName << std::endl;
		return true;
	}

	bool MaskApp::SaveConfig(const std::string& filename)
	{
		std::cout << "Writing configuration to " << filename << std::endl;
		return ConfigSerializer::SerializeConfig(filename, m_params);
	}
	
	void MaskApp::Run()
	{
		std::cout<<"Running simulation..."<<std::endl;
		if(m_systemList.size() != m_params.nThreads)
		{
			std::cerr << "System list not equal to number of threads" << std::endl;
			return;
		}

		//Give our thread pool some tasks
		for(std::size_t i=0; i<m_systemList.size(); i++)
		{
			//bind a lambda to the job, taking in a ReactionSystem, and then provide a reaction system as the tuple arguments.
			m_resources->PushJob({[this](ReactionSystem* system, uint64_t chunkSamples) 
				{
					if(system == nullptr)
						return;

					for(uint64_t i=0; i<chunkSamples; i++)
					{
						system->RunSystem();
						m_fileWriter.PushData(*(system->GetNuclei()));
					}
				}, 
			{m_systemList[i], m_chunkSamples[i]}});
		}

		uint64_t count = 0;
		double percent = 0.05;
		uint64_t flushVal = m_params.nSamples*percent;
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
	
}
