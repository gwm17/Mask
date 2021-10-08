#ifndef MASKAPP_H
#define MASKAPP_H

#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"
#include "ThreeStepSystem.h"
#include "RxnType.h"

namespace Mask {

	class MaskApp {
	public:
		MaskApp();
		~MaskApp();
		bool LoadConfig(const std::string& filename);
		bool SaveConfig(const std::string& filename);
		inline int GetNumberOfSamples() const { return m_nsamples; };
		inline const std::string GetSystemName() const { return sys == nullptr ? "" : sys->GetSystemEquation(); };
		inline const std::string GetOutputName() const { return m_outfile_name; };
		inline const RxnType GetReactionType() const { return m_rxn_type; };
		void Run();

	private:
		void RunOneStepRxn();
		void RunOneStepDecay();
		void RunTwoStep();
		void RunThreeStep();
	
		ReactionSystem* sys;
	
		std::string m_outfile_name;
	
		RxnType m_rxn_type;
		uint32_t m_nsamples;
	
	};

}

#endif
