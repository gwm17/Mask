#ifndef SABREEFFICIENCY_H
#define SABREEFFICIENCY_H

#include "SabreDetector.h"
#include "Target.h"

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
	std::vector<SabreDetector> detectors;
    std::vector<double> ringxs, ringys, ringzs;
    std::vector<double> wedgexs, wedgeys, wedgezs;
	Target deadlayer;


	//Sabre constants
	const double INNER_R = 0.0326;
    const double OUTER_R = 0.1351;
    const double TILT = 40.0;
    //const double DIST_2_TARG = 0.14549;
    const double DIST_2_TARG = -0.1245;
    const double PHI_COVERAGE = 54.4; //delta phi for each det
    const double PHI0 = 234.0; //center phi values for each det in array
    const double PHI1 = 162.0; //# is equal to detID in channel map
    const double PHI2 = 306.0;
    const double PHI3 = 18.0;
    const double PHI4 = 90.0;
    const double DEG2RAD = M_PI/180.0;
    static constexpr double DEADLAYER_THIN = 50 * 1e-7 * 2.3296 * 1e6; // ug/cm^2 (50 nm thick * density)

    const double ENERGY_THRESHOLD = 0.2; //in MeV

};

#endif
