#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"
#include "ThreeStepSystem.h"
#include "Plotter.h"

#include <TFile.h>
#include <TTree.h>

namespace Mask {
//For tree
struct NucData {
	double KE = -1;
	double E = -1;
	double Ex = -1;
	double p = -1;
	double theta = -1;
	double theta_cm = -1;
	double phi = -1;
	int Z = -1;
	int A = -1;
};

class Kinematics {
public:
	Kinematics();
	~Kinematics();
	bool LoadConfig(const char* filename);
	bool SaveConfig(const char* filename);
	inline void SetTreeFlag() { save_tree_flag = true; };
	inline void SetPlotterFlag() { do_plotter_flag = true; };
	inline int GetNumberOfSamples() { return m_nsamples; };
	inline const char* GetSystemName() { return sys == nullptr ? "" : sys->GetSystemEquation().c_str(); };
	inline const char* GetOutputName() { return m_outfile_name.c_str(); };
	inline int GetReactionType() { return m_rxn_type; };
	void Run();

	enum RxnType {
		ONESTEP_RXN,
		ONESTEP_DECAY,
		TWOSTEP,
		THREESTEP
	};

private:
	void RunOneStepRxn();
	void RunOneStepDecay();
	void RunTwoStep();
	void RunThreeStep();

	ReactionSystem* sys;
	Plotter plotman;

	NucData ConvertNucleus(const Nucleus& nuc);

	bool save_tree_flag, do_plotter_flag;

	std::string m_outfile_name;

	int m_rxn_type, m_nsamples;

	TRandom3* global_generator;
};

};

#endif
