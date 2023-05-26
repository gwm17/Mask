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

	struct AppParameters
	{
		std::string outputFileName = "";
		uint64_t nSamples = 0;
		uint32_t nThreads = 1;
		std::vector<StepParameters> chainParams;
		LayeredTarget target;
	};

	class MaskApp
	{
	public:
		MaskApp();
		~MaskApp();
		bool LoadConfig(const std::string& filename);
		bool SaveConfig(const std::string& filename);

		void Run();

	private:
		AppParameters m_params;

		std::vector<ReactionSystem*> m_systemList; //One system for each thread
		std::vector<uint64_t> m_chunkSamples;
		FileWriter m_fileWriter;
		std::unique_ptr<ThreadPool<ReactionSystem*, uint64_t>> m_resources;
	
	};

}

#endif
