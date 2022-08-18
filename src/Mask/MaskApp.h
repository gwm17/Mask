#ifndef MASKAPP_H
#define MASKAPP_H

#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"
#include "ThreeStepSystem.h"
#include "RxnType.h"

namespace Mask {

	class MaskApp
	{
	public:
		MaskApp();
		~MaskApp();
		bool LoadConfig(const std::string& filename);
		bool SaveConfig(const std::string& filename);
		int GetNumberOfSamples() const { return m_nsamples; }
		const std::string GetSystemName() const { return m_system == nullptr ? "Invalid System" : m_system->GetSystemEquation(); }
		const std::string GetOutputName() const { return m_outputName; }
		const RxnType GetReactionType() const { return m_rxnType; }

		void Run();

	private:
		ReactionSystem* m_system;
		std::string m_outputName;
		RxnType m_rxnType;
		int m_nsamples;
	
	};

}

#endif
