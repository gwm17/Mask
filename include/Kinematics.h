#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"
#include "ThreeStepSystem.h"

#include <random>

namespace Mask {

	class Kinematics {
	public:
		Kinematics();
		~Kinematics();
		bool LoadConfig(const std::string& filename);
		bool SaveConfig(const std::string& filename);
		inline int GetNumberOfSamples() { return m_nsamples; };
		inline const std::string GetSystemName() { return sys == nullptr ? "" : sys->GetSystemEquation(); };
		inline const std::string GetOutputName() { return m_outfile_name; };
		inline const int GetReactionType() { return m_rxn_type; };
		void Run();
	
		enum RxnType {
			ONESTEP_DECAY,
			ONESTEP_RXN,
			TWOSTEP,
			THREESTEP
		};
	
	private:
		void RunOneStepRxn();
		void RunOneStepDecay();
		void RunTwoStep();
		void RunThreeStep();
	
		ReactionSystem* sys;
	
		std::string m_outfile_name;
	
		int m_rxn_type, m_nsamples;
	
		std::mt19937* global_generator;
	};

}

#endif
