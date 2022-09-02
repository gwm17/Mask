#ifndef MASKAPP_H
#define MASKAPP_H

#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"
#include "ThreeStepSystem.h"
#include "RxnType.h"
#include "ThreadPool.h"
#include "FileWriter.h"

#include <memory>

namespace Mask {

	class MaskApp
	{
	public:
		MaskApp();
		~MaskApp();
		bool LoadConfig(const std::string& filename);
		bool SaveConfig(const std::string& filename);

		void Run();
		void RunSingleThread();

	private:
		void RunChunk(ReactionSystem* system);

		ReactionSystem* m_system;

		std::string m_outputName;
		RxnType m_rxnType;
		uint64_t m_nsamples;
		uint32_t m_nthreads;

		std::vector<ReactionSystem*> m_systemList; //One system for each thread
		FileWriter m_fileWriter;
		std::unique_ptr<ThreadPool> m_resources;
	
	};

}

#endif
