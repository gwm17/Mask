#ifndef SABREEFFICIENCY_H
#define SABREEFFICIENCY_H

#include "SabreDetector.h"

class SabreEfficiency {
public:
	SabreEfficiency();
	~SabreEfficiency();
	inline void SetReactionType(int t) { m_rxn_type = t; };
	void CalculateEfficiency(const char* file);

private:
	void Run2Step(const char*);
	void Run3Step(const char*);
	void RunDecay(const char*);

	int m_rxn_type;
	std::vector<SabreDetGeometry> detectors;


	//Sabre constants
	const double INNER_R = 0.0326;
    const double OUTER_R = 0.1351;
    const double TILT = 40.0;
    //const double DIST_2_TARG = 0.14549;
    const double DIST_2_TARG = 0.1245;
    const double PHI_COVERAGE = 54.4; //delta phi for each det
    const double PHI0 = 36.0; //center phi values for each det in array
    const double PHI1 = 108.0; //# is equal to detID in channel map
    const double PHI2 = 324.0;
    const double PHI3 = 252.0;
    const double PHI4 = 180.0;
    const double DEG2RAD = M_PI/180.0;

    const double ENERGY_THRESHOLD = 0.1; //in MeV

};

#endif